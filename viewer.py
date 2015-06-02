import wx
from wx.lib.agw.floatspin import FloatSpin
import h5py
import numpy as np
import scipy.linalg as la
from scipy import fftpack
from scipy.stats import norm
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy.optimize import nnls
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import fftconvolve
import multiprocessing 
import ctypes
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.pyplot import cm
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from functools import partial

########################################################################################################
#
#			 Fitting scans to exponents
#
########################################################################################################



########################################## Fitting routines ############################################

def expfit(y, deg):
	"""
	Using Prony's method to fit a sum of exponentials. Adaptation of the procedure published at 
	http://scipy-central.org/item/31/2/using-pronys-method-to-fit-a-sum-of-exponentials

	Fit the data in <y> with $\sum_{j=1}^{deg} c_j exp(a_j (x-x_initial))$
	
	this function is useful to get initial guess for <curve_fit>.
	"""
	# Making sure that y is float
	y = y.astype(dtype=np.float)
	#Build matrix
	A=la.hankel(y[:-deg],y[-deg-1:])
	a=-A[:,:deg]    
	b= A[:,deg]
	#Solve it
	s=la.lstsq(a,b)[0]
	#Solve polynomial
    	p=np.flipud(np.hstack((s,1)))
	u=np.roots(p)
    	#Calc exponential factors
    	a=np.log(u + 0j).real
	u = np.real(u)
	#Build power matrix
    	A=np.power((np.ones((len(y),1))*u[:,None].T),np.arange(len(y))[:,None]*np.ones((1,deg)))
    	#solve it and calculate amplitudes (we enforce amplitudes' positivity)
    	c=nnls(A,y)[0]
    	return a, c

def fitted_histogram (sfft_response_function, background_signal, respose_fun_offset, *params) :
	"""
	This function is used in <post_process_histogram> to curve fit 
	Arguments: <sfft_respose_function> -- the symmetric FFT of response function all other characterize the sum of exps. 
	"""
	# Calculating exponential sum
	N = sfft_response_function.size
	exp_sum = np.zeros(N, dtype=np.complex)
	exp_sum += background_signal
	x = np.arange(N)

	for a, c in zip( params[:len(params)/2], params[len(params)/2:] ) : exp_sum += c*np.exp(a*x)
	
	# Shitting the response function using the symmetric FFT shift theorem
	response_function = np.exp( -2j*np.pi * respose_fun_offset * fftpack.fftshift(fftpack.fftfreq(N)) ) 
	response_function *= sfft_response_function
	response_function = fftpack.ifftshift(fftpack.ifft( fftpack.ifftshift(response_function), overwrite_x=True)) 
	
	# Performing circular convolution of the response function with the sum of exponentials
	exp_sum = fftpack.fft(exp_sum, overwrite_x=True)
	exp_sum *= fftpack.fft(response_function, overwrite_x=True)

	return fftpack.ifft(exp_sum, overwrite_x=True).real

def post_process_histogram (histogram, response_function, deg, get_fvec = False) :
	'''
	This function fits a histogram to the convolution of response function with the sum of exp.
	'''
	# Check weather the histogram is absent
	if not isinstance(histogram, np.ndarray) :
		nan_array = np.array([ np.nan for k in range(get_fit_params_size(deg)) ])
		if get_fvec : return np.array([]), nan_array
		else : return nan_array

	# Using Prony method to get the initial guess for the exponential fit 
	
	# Preparing histogram for Prony fitting
	initial_point = np.argmax(histogram)
	try : final_point = np.nonzero(histogram > 2)[0][-1]
	except IndexError : final_point = initial_point

	histogram_cut = histogram[initial_point:final_point].astype(dtype=np.float)

	# Checking whether the cut histogram is not empty
	if histogram_cut.size <= 1 : 
		nan_array = np.array([ np.nan for k in range(get_fit_params_size(deg)) ])
		if get_fvec : return np.array([]), nan_array
		else : return nan_array

	# Get the initial guess
	exponent, ampl = expfit(gaussian_filter(histogram_cut, sigma=0.02*histogram_cut.size), deg)

	# Accounting for <initial_point> in the initial guess
	ampl *= np.exp(-exponent*initial_point)
	p0 = np.append(exponent, ampl)
	
	# Adjusting sizes of the histogram and response function
	if response_function.size > histogram.size : response_function = response_function[:histogram.size]
	else : histogram = histogram[:response_function.size]

	# Saving the symmetric FFT of taylored expression
	sfft_response_function = fftpack.fftshift( fftpack.fft( fftpack.fftshift(response_function) ) )

	try :
		# Perform accurate fit accounting for the response function
		params, pcov = curve_fit(fitted_histogram, sfft_response_function, histogram, 
								p0=np.append([0, 0],p0), maxfev=5000)

		# Perform the chi-square test to characterize the goodness of fit and append the results to <fit_params> 
		fitted_func = fitted_histogram(sfft_response_function, *params)
		# Note that the number of fitting parameters (i.e., number of arguments of <fitted_histogram>) is used as the delta degree of fredom
		fit_params = np.append(chisquare(histogram, fitted_func, ddof=(2*deg+2)), params)
	except RuntimeError, err :
		fit_params = np.array([ np.nan for k in range(get_fit_params_size(deg)) ])
		fitted_func = np.array([])
		print "RuntimeError : %s " % str(err)

	if get_fvec : return fitted_func, fit_params # return the fitted function and fitting parameters
	else : return fit_params # return only fitting parameters

########################### Utilities interpreting returns of <post_process_histogram> ###############

def get_chisq (fit_params) :
	"""
	chi-square characterizing the goodness of fit
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,0]
	else : return fit_params[0]

def get_p_val_chisq_test (fit_params) :
	"""
	probability value of the chi-square test (i.e., the probability that the fit is wrong)
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,1]
	else : return fit_params[1]

def get_background_signal (fit_params) :
	"""
	background noise (i.e., additive constant)
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,2]
	else : return fit_params[2]

def get_respose_fun_offset (fit_params) :
	"""
	response function offset
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,3]
	else : return fit_params[3]

def get_respose_fun_offset_ps (fit_params, Resolution) :
	"""
	Return the picosecond value of  the response function offset, obtained as a result of fitting procedure.
	"""
	return Resolution * get_respose_fun_offset(fit_params)
			
def get_exp_factors (fit_params) :
	"""
	time-life
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,4:(fit_params.shape[1]/2 +2)]
	else : return fit_params[4:(len(fit_params)/2 +2)]

def get_time_life_ns (fit_params, Resolution) :
	"""
	Convert exp factors, extracted during fitting into time life in nano seconds.
	"""
	return 	(1.e-3 * Resolution) / np.abs(get_exp_factors(fit_params))

def get_exp_ampl (fit_params) :
	"""
	pre-exponential amplitude
	"""
	if isinstance(fit_params, np.ndarray) and len(fit_params.shape) == 2 : 
		return fit_params[:,(fit_params.shape[1]/2 +2):]
	return fit_params[(len(fit_params)/2 +2):]

def get_fit_params_size (deg) :
	"""
	number of parameters to represent the fit with multiple (<deg>) exponents.
	"""
	return 2*deg+4

##########################################################################################################
def post_process_scans (file_name, response_function, histogram_position, deg, processed_scans_buffer, post_process_event) :
	"""
	This function will be run as a separate process. Post process entire 3D scan stored in file <file_name>.
	<processed_scans_buffer> saves the results of calculations.
	Communications are done through <post_process_event>
	<histogram_position> contains the coordinates of histograms to be processed.
	"""
	print "\nBegin post processing...\n"
	
	data_file = h5py.File (file_name, 'r')
	histogram_group = data_file["histograms"]

	# Create a numpy wrapper over <processed_scans_buffer>
	processed_scans = np.frombuffer(processed_scans_buffer.get_obj())
	processed_scans = processed_scans.reshape( (len(histogram_position), get_fit_params_size(deg))  )

	# Pool of workers
	pool_size = multiprocessing.cpu_count()-2
	pool = multiprocessing.Pool (processes=pool_size)
	
	# Create map function
	map_func = partial(post_process_histogram, response_function=response_function, deg=deg, get_fvec=False)

	# Schedule of loading the histograms 
	Loads = range(0, histogram_position.shape[0], 100*pool_size)	
	del Loads[0]
	if len(Loads) == 0 : Loads = [ histogram_position.shape[0] ]
	elif Loads[-1] != histogram_position.shape[0] : Loads.append (histogram_position.shape[0])
	
	L_previous = 0
	for L in Loads :
		# Check weather the user requested to stop post processing
		if not post_process_event.is_set() : break
		
		# Load histograms from specified positions
		histograms = [ load_histogram(histogram_group, p) for p in histogram_position[L_previous:L] ]
		
		# Running post_processing on multiple cores
		map_results = pool.map(map_func, histograms)

		# Saving results
		processed_scans_buffer.acquire()
		processed_scans[L_previous:L] = map_results
		processed_scans_buffer.release()

		# Cleaning lists
		del histograms, map_results
		
		L_previous = L

	# Signal that post processing is over
	post_process_event.clear ()
	data_file.close()
	
	print "\nPost processing completed\n"

############################################ Utilities ######################################################

def load_response_function (Resolution, path="") :
	"""
	Load response function that corresponds to histograms with <Resolution>. 
	<path> list of directories, where search the file should be searched for.
	"""	
	filename = "response_function_%dps.dat" % Resolution  
	
	# Trying to load the file from the list of paths
	if isinstance(path, basestring) : path = [path]
	
	for P in path :
		try : response_function = np.fromfile (P + filename)
		except IOError : pass
	
	try : response_function 
	except NameError :
		# The file has not been found so let's ask user for a directory
		openFileDialog = wx.FileDialog (None, "Chose file %s" % filename, "", "", \
			"DATA (*.data)|*.dat", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
		#"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

		# Check whether user canceled
		if openFileDialog.ShowModal() == wx.ID_CANCEL : 
				#print ("User did not chose file containing the response function")
				response_function = None
		else :
			response_function = np.fromfile (openFileDialog.GetPath())

	return response_function
	
def get_file_checksum (filename, path="") :
	"""
	Calculate a checksum of a file specified by <filename>. 
	<path> list of directories, where search the file should be searched for.
	"""
	# Trying to load the file from the list of paths
	if isinstance(path, basestring) : path = [path]
	
	for P in path :
		try : 
			file_to_analize = open (P + filename, 'rb'); break
		except IOError : pass
	
	try : file_to_analize 
	except NameError :
		# The file has not been found so let's ask user for a directory
		openFileDialog = wx.FileDialog (None, "Chose file %s" % filename, "", "", \
			"ALL (*.*)|*.*", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

		# Check whether user canceled
		if openFileDialog.ShowModal() == wx.ID_CANCEL : 
				raise ValueError ("User did not chose the file")
	
	# Finally calculating hash
	import hashlib
	m = hashlib.sha1()	
	for chunk in file_to_analize : m.update(chunk)
	file_to_analize.close()
	return m.hexdigest()
	
def array_to_materials (A, N=None) :
	"""
	Convert array <A> into a list of materials. 
	<N> a number of materials to be used, i.e., the level of quantization
	"""		
	# Normalize a copy of the array
	A = np.array(A, dtype=np.float, copy=True)
	A -= A.min(); A /= A.max()

	# Convert to the desired quantization
	if np.isscalar (N) : 
		N = np.abs(N)
		A *= N; np.round(A, out=A)
		N -= 1
		if N > 1 : A /= N 
		
	# find unique elements and their incidences
	U, U_indecies = np.unique (A, return_inverse=True)
				
	# Create a list of materials, which are dictionaries of properties 
	materials = [ { "diffuse_color" : tuple(color.tolist())[:3], "specular_intensity" : 0.3, "alpha" : 0.3	} 
					for color in cm.jet(U) ]
	
	# Create a dictionary of materials
	materials = dict( ( "material_%d" % num, material ) for num, material in enumerate(materials) )
	
	# Convert array A to materials 
	array_converted_to_materials = [ "material_%d" % num for num in U_indecies ]
	return materials, array_converted_to_materials

def load_histogram (histograms_group, position) :
	"""
	Load histogram in the specified triple <position> from <histograms_group>. If the histogram is absent then <np.nan> will be returned
	"""
	try : 
		histogram = histograms_group["histogram_%d_%d_%d" % tuple(position)][...]
		# Deleting zeros at the end
		mask = (histogram > 0 )
		if mask.sum() == 0 :
			# Istead of returing a zero size array, we will return a nan value
			raise KeyError
		cut_off = np.nonzero(mask)[0][-1]
		return histogram[:cut_off]
	except (KeyError, IOError) : 
		return np.nan
		
#####################################################################################################################
#
#				Viewer
#
#####################################################################################################################

class CViewer (wx.Frame) :
	"""
	Application for viewing the data taken by microscope
	"""
	
	def __init__ (self, parent=None, filename = None, title=None) :
		"""
		Constructor
		"""
		# Saving the current path
		import os, sys
		self.paths = [ "" ]
		self.paths.append( os.path.abspath(os.path.dirname(sys.argv[0])) + "\\" )
			
		# If file name is not specified, then ask user what file should be open 
		if filename == None :
			openFileDialog = wx.FileDialog (parent, "Chose HDF5 file for viewing", "", "", \
				"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)

			# Check whether user canceled
			if openFileDialog.ShowModal() == wx.ID_CANCEL : 
				raise ValueError ("User did not chose the file to view")
			else : filename = openFileDialog.GetPath()
		
		# Open file
		self.data_file = h5py.File (filename, 'r')
		
		# Adding new path
		self.paths.append( os.path.abspath('') + "\\" )
		
		# Loading groups
		self.parameters_group 	= self.data_file["parameters"]
		self.histograms_group	= self.data_file["histograms"]
		try : self.photos_group	= self.data_file["photos"]
		except KeyError : pass
		
		self.Resolution	= int(self.parameters_group["Resolution"][...])

		# Loading response function
		self.response_function = load_response_function(self.Resolution, self.paths)

		# Read some parameters
		try :
			# This is for backward compatibilities
			#self.piezo_calibrated_displacement_ax1 = int(self.parameters_group["piezo_calibrated_displacement_ax1"][...])
			#self.piezo_calibrated_displacement_ax2 = int(self.parameters_group["piezo_calibrated_displacement_ax2"][...])
			#self.piezo_calibrated_displacement_ax3 = int(self.parameters_group["piezo_calibrated_displacement_ax3"][...])
			self.piezo_calibrated_displacement_ax1 = int(self.parameters_group["unit_size"][...])
			self.piezo_calibrated_displacement_ax2 = self.piezo_calibrated_displacement_ax1
			self.piezo_calibrated_displacement_ax3 = int(self.parameters_group["piezo_calibrated_displacement_ax3"][...])
		except KeyError : 
			self.piezo_calibrated_displacement_ax1 = int(self.parameters_group["unit_size"][...])
			self.piezo_calibrated_displacement_ax2 = self.piezo_calibrated_displacement_ax1
			self.piezo_calibrated_displacement_ax3 = self.piezo_calibrated_displacement_ax1

		try :
			# If the HDF5 file contains predicated the rate of fluorescence, 
			# then load them otherwise calculate manually 
			self.voltages 			= self.parameters_group["fluoresence_rate_keys"][...]
			self.fluoresence_rate	= self.parameters_group["fluoresence_rate_values"][...]
		except KeyError :
			# Crating voltages by looking at the name of the histograms
			self.voltages = set( tuple(map(int, key.split('_')[-3:])) for key in self.histograms_group.iterkeys() )
			self.voltages = np.array(list(self.voltages))
			
			# Loading fluorescence rate 
			self.fluoresence_rate = np.array( map( lambda p : np.sum(load_histogram(self.histograms_group, p)), self.voltages), dtype=np.float )
			self.fluoresence_rate /= self.parameters_group["Tacq"][...]
			
			# Deleting all elements without histogram 
			idx_to_delete = np.nonzero( np.logical_not(np.isnan(self.fluoresence_rate)) )
			self.voltages = np.asarray( self.voltages[idx_to_delete], order='C')
			self.fluoresence_rate = np.asarray( self.fluoresence_rate[idx_to_delete], order='C')
		
		# Convert voltages to microns
		self.x_micron = 1.e-3 * self.piezo_calibrated_displacement_ax1 * self.voltages[:,0]
		self.y_micron = 1.e-3 * self.piezo_calibrated_displacement_ax2 * self.voltages[:,1]
		self.z_micron = 1.e-3 * self.piezo_calibrated_displacement_ax3 * self.voltages[:,2]
		
		# Create GUI
		if title == None : 
			import ntpath
			title = "Viewing file [%s]" % ntpath.basename(filename)

		wx.Frame.__init__ (self, parent, title=title)
		self.ConstructGUI ()
		self.Show ()
		self.Maximize (True)

	def __del__ (self) :
		# Stopping post processing
		try : self.post_process_event.clear()
		except AttributeError : pass

	def ConstructGUI (self) :
		"""
		Build GUI
		"""
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer (2, 4)

		################################ Display all measurements in the file ########################
		# Matplotlib canvas to display fluorescence rate 	
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		self.dpi = 80
		display_width, display_hight = wx.DisplaySize()

		self.fluoresence_rate_fig = Figure((0.49*display_width/self.dpi, 0.8*display_hight/self.dpi), dpi=self.dpi)
		self.fluoresence_rate_canvas = FigCanvas (panel, -1, self.fluoresence_rate_fig)
		
		# Process canvas clicking by displaying histogram
		self.fluoresence_rate_canvas.mpl_connect ('pick_event', partial(self.show_histogram, is_animation=False))

		# Plot fluorescence rate
		self.display_fluoresence_rate ()

		boxsizer.Add(self.fluoresence_rate_canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
	
		toolbar = NavigationToolbar (self.fluoresence_rate_canvas)
		boxsizer.Add(toolbar, 0, wx.EXPAND)
				
		sizer.Add(boxsizer, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)			
		
		################################### Viewing region #########################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)
	
		boxsizer.Add (wx.StaticText(panel, label="Viewing region: value 1 min"))
		self.value1min = FloatSpin (panel, increment=0.1, digits=3)
		boxsizer.Add (self.value1min)
		
		boxsizer.Add (wx.StaticText(panel, label=" value 1 max"))
		self.value1max = FloatSpin (panel, increment=0.1, digits=3)
		boxsizer.Add (self.value1max)

		sizer.Add (boxsizer, pos=(1, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)			

		####################################### Commands (left panel) ###########################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

		# Plot the time life distribution
		self.time_life_distribution_button = wx.Button (panel, label="Time life distribution")

		self.time_life_distribution_button.__get_data__ = partial (get_time_life_ns, Resolution=self.Resolution)
		self.time_life_distribution_button.__histogram_xlabel__ 	= "decay time [ns]"
		self.time_life_distribution_button.__histogram_title__ 		= "Distribution of %d decay time (mean = %.2e ; sigma = %.2e)"
		self.time_life_distribution_button.__histogram_exception_title__= "Distribution of time decay"
		self.time_life_distribution_button.__3d_plot_title__ 		= "%d decay time [ns]"
		self.time_life_distribution_button.__3d_plot_exception_title__	= "Time life [ns]"
		
		self.time_life_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.time_life_distribution_button.Bind (wx.EVT_LEFT_DCLICK, self.display_post_processing_data)
		boxsizer.Add (self.time_life_distribution_button, flag=wx.LEFT, border=5)

		# Distribution of amplitudes in front of exponents
		self.exp_ampl_distribution_button = wx.Button (panel, label="Exponential amplitude distribution")
		
		self.exp_ampl_distribution_button.__get_data__ 			= get_exp_ampl
		self.exp_ampl_distribution_button.__histogram_xlabel__ 		= "Exponential amplitude (counts)" 
		self.exp_ampl_distribution_button.__histogram_title__		= "Distribution of %d exp amplitude (mean = %.2e ; sigma = %.2e)"
		self.exp_ampl_distribution_button.__histogram_exception_title__	= "Distribution of exp ampl"
		self.exp_ampl_distribution_button.__3d_plot_title__		= "Amplitude in front of %d exponent (counts)"
		self.exp_ampl_distribution_button.__3d_plot_exception_title__	= "Amplitude in front of exponent (counts)"

		self.exp_ampl_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.exp_ampl_distribution_button.Bind (wx.EVT_LEFT_DCLICK, self.display_post_processing_data)
		boxsizer.Add (self.exp_ampl_distribution_button, flag=wx.LEFT, border=5)

		# Background noise distribution (i.e., additive constant)
		self.background_noise_distribution_button = wx.Button (panel, label="Background noise distribution")
		
		self.background_noise_distribution_button.__get_data__ = get_background_signal
		self.background_noise_distribution_button.__histogram_xlabel__ 		= "Additive noise [counts]"
		self.background_noise_distribution_button.__histogram_title__		= ""
		self.background_noise_distribution_button.__histogram_exception_title__ = "Distribution of additive noise"
		self.background_noise_distribution_button.__3d_plot_title__		= ""
		self.background_noise_distribution_button.__3d_plot_exception_title__	= "Background additive noise [counts]"

		self.background_noise_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.background_noise_distribution_button.Bind (wx.EVT_LEFT_DCLICK, self.display_post_processing_data)
		boxsizer.Add (self.background_noise_distribution_button, flag=wx.LEFT, border=5)
	
		sizer.Add (boxsizer, pos=(2, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)	

		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

		# PicoHarp time jitter (PicoHarp does not maintain the stable t=0 so we had to offset the response function during fitting)
		self.time_jitter_distribution_button = wx.Button (panel, label="Time jitter distribution")

		self.time_jitter_distribution_button.__get_data__ = partial(get_respose_fun_offset_ps, Resolution=self.Resolution)
			
		self.time_jitter_distribution_button.__histogram_xlabel__		= "Time offset [ps]"
		self.time_jitter_distribution_button.__histogram_title__		= ""
		self.time_jitter_distribution_button.__histogram_exception_title__	= "Distribution of time offset"
		self.time_jitter_distribution_button.__3d_plot_title__			= ""
		self.time_jitter_distribution_button.__3d_plot_exception_title__	= "Time jitter [ps]"

		self.time_jitter_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.time_jitter_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)

		boxsizer.Add (self.time_jitter_distribution_button, flag=wx.LEFT, border=5)

		# Chi-square distribution
		self.chisq_distribution_button = wx.Button (panel, label="Chi-square distribution") 

		self.chisq_distribution_button.__get_data__ 			= get_chisq
		self.chisq_distribution_button.__histogram_xlabel__		= "Chi-square"
		self.chisq_distribution_button.__histogram_title__		= ""
		self.chisq_distribution_button.__histogram_exception_title__	= "Distribution of chi-square"
		self.chisq_distribution_button.__3d_plot_title__		= ""
		self.chisq_distribution_button.__3d_plot_exception_title__	= "Chi-square"

		self.chisq_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.chisq_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)

		boxsizer.Add (self.chisq_distribution_button, flag=wx.LEFT, border=5)

		# P-values of chi-square test (i.e., the probability that the fit is wrong)
		self.pval_distribution_button = wx.Button (panel, label="P-value of chisquare test")

		self.pval_distribution_button.__get_data__ 			= get_p_val_chisq_test
		self.pval_distribution_button.__histogram_xlabel__		= "Probability"
		self.pval_distribution_button.__histogram_title__		= ""
		self.pval_distribution_button.__histogram_exception_title__	= "Distribution of P-value of chi-square test"
		self.pval_distribution_button.__3d_plot_title__			= ""
		self.pval_distribution_button.__3d_plot_exception_title__	= "P-value of chi-square test"

		self.pval_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
		self.pval_distribution_button.Bind (wx.EVT_LEFT_DOWN, self.display_post_processing_data)
	
		boxsizer.Add (self.pval_distribution_button, flag=wx.LEFT, border=5)
		
		# Draw fluorescence rate button
		self.fluoresence_rate_button = wx.Button (panel, label="fluorescence rate")
		self.fluoresence_rate_button.__get_data__ 		= lambda K : self.fluoresence_rate
		self.Bind (wx.EVT_BUTTON, self.display_fluoresence_rate, self.fluoresence_rate_button)
		boxsizer.Add(self.fluoresence_rate_button,  flag=wx.LEFT, border=5)

		sizer.Add (boxsizer, pos=(3, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)	

		############################# Setting up matplot lib #######################################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		self.histogram_fig = Figure((0.49*display_width/self.dpi, 0.8*display_hight/self.dpi), dpi=self.dpi)
		self.histogram_canvas = FigCanvas (panel, -1, self.histogram_fig)
		boxsizer.Add(self.histogram_canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
	
		toolbar = NavigationToolbar (self.histogram_canvas)
		boxsizer.Add(toolbar, 0, wx.EXPAND)
		
		sizer.Add(boxsizer, pos=(0, 1), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)

		################################### Settings ########################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

		# Number of exponent to be used in fitting
		boxsizer.Add(wx.StaticText(panel, label="Number of exponents to fit"), flag=wx.LEFT, border=5)  
		self.number_fitting_exponent = wx.SpinCtrl (panel, value="1", min=1, max=5)
		def fitting_deg_updated (event) :
			# Delete all fits
			try : del self.processed_scans_buffer, self.processed_scans
			except AttributeError : pass
			# redraw the current histogram
			self.show_histogram ()
		self.number_fitting_exponent.Bind (wx.EVT_SPINCTRL, fitting_deg_updated)
		boxsizer.Add(self.number_fitting_exponent, flag=wx.LEFT, border=5)

		self.to_process_histogram = wx.CheckBox (panel, label="Fit histogram to exponent")
		self.to_process_histogram.SetValue(True)
		self.Bind (wx.EVT_CHECKBOX, self.show_histogram, self.to_process_histogram)
		boxsizer.Add(self.to_process_histogram,  flag=wx.LEFT, border=5)

		# Disable exp fitting controls if response function has not been loaded 
		if self.response_function is None :
			self.number_fitting_exponent.SetValue(0)
			self.number_fitting_exponent.Disable()
			
			self.to_process_histogram.SetValue(False)
			self.to_process_histogram.Disable()
			
		self.to_use_logplot = wx.CheckBox (panel, label="Use log plot for histogram")
		self.to_use_logplot.SetValue(True)
		self.Bind (wx.EVT_CHECKBOX, self.show_histogram, self.to_use_logplot)
		
		boxsizer.Add(self.to_use_logplot,  flag=wx.LEFT, border=5)
		sizer.Add(boxsizer, pos=(1, 1), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)
		
		###################################### Commands (right panel) #########################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

		# Buttons to save fitted info that can be edited in Blender
		self.save_post_processed_button = wx.Button (panel, label="Save FLIM...") 
		self.Bind (wx.EVT_BUTTON, self.save_post_processed, self.save_post_processed_button)
		boxsizer.Add(self.save_post_processed_button, flag=wx.LEFT, border=5)
		
		# Buttons to export pickled fit data 
		self.load_post_processed_button = wx.Button (panel, label="Load FLIM...") 
		self.Bind (wx.EVT_BUTTON, self.load_post_processed, self.load_post_processed_button)
		boxsizer.Add(self.load_post_processed_button, flag=wx.LEFT, border=5)
		
		# Extract boundary 
		self.extract_boundary_button = wx.Button (panel, label="Extract boundary...")
		self.Bind (wx.EVT_BUTTON, self.extract_boundary, self.extract_boundary_button)
		boxsizer.Add(self.extract_boundary_button, flag=wx.LEFT, border=5)
		
		sizer.Add(boxsizer, pos=(2, 1), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)
		#################################### Animation ##########################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

		# Go to initial frame button
		self.initial_frame_button = wx.Button (panel, label="<<")
		def go_to_initial_frame (event) :
			self.current_frame.SetValue (0); self.show_histogram ()
		self.Bind (wx.EVT_BUTTON, go_to_initial_frame, self.initial_frame_button)
		boxsizer.Add(self.initial_frame_button, flag=wx.LEFT, border=5)

		# Go to the previous frame
		self.previous_frame_button = wx.Button (panel, label="<")
		def go_to_previous_frame (event) :
			current_value = self.current_frame.GetValue()
			if current_value > 0 : 
				self.current_frame.SetValue(current_value-1); self.show_histogram ()
		self.Bind (wx.EVT_BUTTON, go_to_previous_frame, self.previous_frame_button)
		boxsizer.Add(self.previous_frame_button, flag=wx.LEFT, border=5)
		
		# Variable storing current frame number
		self.current_frame = wx.SpinCtrl (panel, value="0", min=0, max=len(self.fluoresence_rate)-1)
		self.current_frame.Bind (wx.EVT_SPINCTRL, self.show_histogram)
		boxsizer.Add(self.current_frame, flag=wx.LEFT, border=5)

		# Animation button
		self.animation_button = wx.Button (panel)
		self.animation_button.__start_label__ = "Start animation"
		self.animation_button.__stop_label__ = "STOP animation"
		self.animation_button.SetLabel (self.animation_button.__start_label__)
		self.Bind (wx.EVT_BUTTON, self.animation, self.animation_button)
		boxsizer.Add(self.animation_button, flag=wx.LEFT, border=5)

		# Go to the next frame button
		self.next_frame_button = wx.Button (panel, label=">")
		def go_to_next_frame (event) :
			current_value = self.current_frame.GetValue()
			if current_value < len(self.fluoresence_rate)-1 : 
				self.current_frame.SetValue(current_value+1); self.show_histogram ()
		self.Bind (wx.EVT_BUTTON, go_to_next_frame, self.next_frame_button)
		boxsizer.Add(self.next_frame_button, flag=wx.LEFT, border=5)

		# Go to the last frame button
		self.final_frame_button = wx.Button (panel, label=">>")
		def go_to_last_frame (event) :
			self.current_frame.SetValue (len(self.fluoresence_rate)-1); self.show_histogram ()
		self.Bind (wx.EVT_BUTTON, go_to_last_frame, self.final_frame_button)
		boxsizer.Add(self.final_frame_button, flag=wx.LEFT, border=5)

		sizer.Add(boxsizer, pos=(3, 1), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)
		#########################################################################################			
	
		panel.SetSizer (sizer)		
		sizer.Fit (panel)
		
		################ Display the first frame ################################
		self.show_histogram ()

	def animation (self, event=None) :
		"""
		<self.animation_button> was clicked
		"""
		if self.animation_button.GetLabel() == self.animation_button.__start_label__ :
			# Initiate animation
			
			# set the timer
			timer_id = wx.NewId()
			self.animation_timer = wx.Timer (self, timer_id)
			self.animation_timer.Start (0, True)
			
			def do_animation (event) :
				current_value = self.current_frame.GetValue()
				if current_value < len(self.fluoresence_rate)-1 : 
					# Continue animation
					self.current_frame.SetValue(current_value+1); wx.CallAfter(self.show_histogram)
					# Setting next timer call 
					self.animation_timer.Start (0, True)
				else :
					# Stop animation
					self.animation()
				
			wx.EVT_TIMER (self, timer_id, do_animation)
			self.animation_button.SetLabel( self.animation_button.__stop_label__ )
		else :
			# Stop animation
			self.animation_timer.Stop()
			del self.animation_timer
			self.animation_button.SetLabel( self.animation_button.__start_label__ ) 

	def show_histogram (self, event=None, is_animation=True) :
		"""
		Image displaying total fluorescent scans is mouse clicked
		""" 
		if is_animation or event.mouseevent.button == 1 : 
			# Left button click
			try : del self.__histogram_axes__, self.__photo_axes__
			except AttributeError : pass
			self.histogram_fig.clear ()
			self.histogram_fig.set_facecolor('grey')
			self.__histogram_axes__ = self.histogram_fig.add_subplot (211)
			self.__photo_axes__ = self.histogram_fig.add_subplot (212)
		elif event.mouseevent.button == 3 :
			# Right button click
			pass
		# Ignore all other clicks
		else : return  

		# Loading histogram 
		if not is_animation : self.current_frame.SetValue(event.ind[0])
		position = tuple(self.voltages[self.current_frame.GetValue(), :])
		histogram = load_histogram( self.histograms_group, position)
		
		# Load photo
		try :
			photo = self.photos_group["photo_%d_%d_%d" % position][...]
			# plot photo
			if isinstance(photo, np.ndarray) :
				self.__photo_axes__.clear()
				self.__photo_axes__.imshow(photo, cmap=cm.gray)
				self.__photo_axes__.set_xlabel ("Sample's photo")
				self.__photo_axes__.get_xaxis().set_ticks([])
				self.__photo_axes__.get_yaxis().set_ticks([])
		except (AttributeError, KeyError) : pass

		if not isinstance(histogram, np.ndarray) :
			# There is nothing to plot
			self.histogram_canvas.draw(); return
	
		# time axis in nano sec
		pixels_to_ns = 1.e-3 * self.Resolution
		time_ns = pixels_to_ns * np.arange(histogram.size)

		# Perform histogram post processing, if requested
		if self.to_process_histogram.IsChecked() : 
			fitted_func, fit_params = post_process_histogram (histogram, 
					self.response_function, self.number_fitting_exponent.GetValue(), True)
			self.__histogram_axes__.set_title ("Processed histogram")	
		else :
			self.__histogram_axes__.set_title ("Raw histogram")	
	
		if self.to_use_logplot.IsChecked() :
			# use log plot
			if 0 in histogram : self.__histogram_axes__.set_ylim ([ 0.1, 10**np.ceil(np.log10(histogram.max())) ])
			PLOTTING_FUNCTION = self.__histogram_axes__.semilogy
		else : # use linear plot
			PLOTTING_FUNCTION = self.__histogram_axes__.plot

		if histogram.sum() == 0 :
			# there is no histogram
			self.__histogram_axes__.plot (time_ns, histogram)
		else :
			PLOTTING_FUNCTION (time_ns, histogram)
			# Draw exponential fit, if requested
			try : 
				# Adjusting time axis 
				time_ns = pixels_to_ns * np.arange(fitted_func.size)
				# Plotting the fit
				PLOTTING_FUNCTION (time_ns, fitted_func)
				self.__histogram_axes__.legend (["histogram", "exp fit"])
			except NameError : pass
			
		self.__histogram_axes__.set_xlabel ('time (ns)')
		self.__histogram_axes__.set_ylabel ('counts')
		self.histogram_canvas.draw()
	
	def display_post_processing_data (self, event) :
		"""
		Visualize data obtained from exponent fitting
		"""
		# Determine which button was clicked
		try : 
			button = event.GetEventObject()
			# Mouse double clicking stops post-processing
			if event.GetEventType() == wx.wxEVT_LEFT_DCLICK  : 
				try : 
					self.post_process_event.clear()
					return
				except AttributeError : pass
		except AttributeError : button = event
	
		# Check whether post processing data is available
		try : self.processed_scans
		except AttributeError :
			# the data is not available and needs to be obtained
			self.post_processing(button)
			return

		# Extract the data to be plotted
		self.processed_scans_buffer.acquire()
		data = button.__get_data__(self.processed_scans)
		# If the property is one parametric then reshape the array
		if len(data.shape) == 1 : data = data.reshape( (data.size,1) )
		# Make a copy of data so that the values can be altered 
		data = np.copy(data)
		self.processed_scans_buffer.release()
		
		# Draw the data histogram if post-processing is over
		try : self.post_processing_process	
		except AttributeError :
			# Post-processing is over
			self.gaussian_fits = []
			self.histogram_fig.clear (); self.histogram_fig.set_facecolor('grey')
			for cut in range(data.shape[1]) :
				data_cut = data[:,cut]
				
				# Removing nan values from histogram
				data_cut = data_cut[ np.nonzero(np.logical_not(np.isnan(data_cut))) ]
				
				# Calculate statistical fits
				self.gaussian_fits.append( norm.fit(data_cut) )
				mean, sigma = self.gaussian_fits[cut]
				
				axes = self.histogram_fig.add_subplot(data.shape[1], 1, cut+1, axisbg='grey')
				axes.hist (data_cut, 100)
				axes.set_xlabel (button.__histogram_xlabel__)
				try : axes.set_title ( button.__histogram_title__ % ((cut+1), mean, sigma) )
				except TypeError : axes.set_title(button.__histogram_exception_title__)
			
			self.histogram_canvas.draw()

			# Since the post processing is over, update the default values of the bounds
			# Do not update if the same button was previously clicked
			try : 
				if self.__previously_clicked_button__ is not button : raise AttributeError 
			except AttributeError :
				self.__previously_clicked_button__ = button
				mean, sigma = self.gaussian_fits[0]
				min_val = mean - 2*sigma; max_val = mean + 2*sigma
				if np.isnan(min_val) or np.isnan(max_val) :
					min_val = np.nanmin(data); max_val = np.nanmax(data); 
				self.value1min.SetValue (min_val); self.value1max.SetValue (max_val)
	
		# Plotting 3D plots 
		self.fluoresence_rate_fig.clear ()
		self.fluoresence_rate_fig.patch.set_facecolor('grey')

		for cut in range(data.shape[1]) :
			
			# Applying cut-off to the data
			data_cut = data[:, cut]
			try :
				mean, sigma = self.gaussian_fits[cut]
				try : 
					min_val = getattr(self, "value%dmin" % (cut+1)).GetValue() 
					max_val = getattr(self, "value%dmax" % (cut+1)).GetValue() 
				except AttributeError : 
					min_val = mean - 2*sigma; max_val = mean + 2*sigma
				data_cut[ np.nonzero(data_cut < min_val) ] = min_val
				data_cut[ np.nonzero(data_cut > max_val) ] = max_val
			except IndexError : pass
		
			axes = self.fluoresence_rate_fig.add_subplot(data.shape[1], 1, cut+1, projection='3d', axisbg='grey')
			axes.set_aspect ('equal')
			axes.view_init (90, 90)
			
			# Show relative values of the coordinates
			X = self.x_micron - self.x_micron.min()
			Y = self.y_micron - self.y_micron.min()
			
			sc = axes.scatter (X, Y, self.voltages[:,2], c=data_cut, cmap=cm.jet, picker=2, linewidths=0)
			
			axes.w_xaxis.set_pane_color((0.5,0.5,0.5))
			axes.w_yaxis.set_pane_color((0.5,0.5,0.5))
			axes.w_zaxis.set_pane_color((0.5,0.5,0.5))

			axes.set_xlabel ('ax1 (micron)')
			axes.set_ylabel ('ax2 (micron)')
			axes.set_zlabel ('ax3 (voltages)')
	
			try : axes.set_title ( button.__3d_plot_title__ % (cut+1) )
			except TypeError : axes.set_title (button.__3d_plot_exception_title__)

			try :
				if data.shape[1] > 1 : self.fluoresence_rate_fig.colorbar (sc, use_gridspec=True, orientation='vertical')
				else : self.fluoresence_rate_fig.colorbar (sc, use_gridspec=True, orientation='horizontal')
			except TypeError : pass

		self.fluoresence_rate_canvas.draw()	

	def post_processing (self, button) :
		"""
		Post process all scans
		"""
		# Begin post processing
		self.post_process_event = multiprocessing.Event()
		self.post_process_event.set()

		deg = self.number_fitting_exponent.GetValue()
			
		# Cleaning up if possible
		self.gaussian_fits = []
		try : del self.processed_scans_buffer, self.processed_scans
		except AttributeError : pass

		# Allocating memory
		buffer_shape = (len(self.voltages), get_fit_params_size(deg))
		self.processed_scans_buffer = multiprocessing.Array (ctypes.c_double, np.prod(buffer_shape) )

		# Creating numpy wrapper over the buffer
		self.processed_scans_buffer.acquire()
		self.processed_scans = np.frombuffer (self.processed_scans_buffer.get_obj())
		self.processed_scans[:] = np.nan
		self.processed_scans = self.processed_scans.reshape( buffer_shape  )
		self.processed_scans_buffer.release()

		self.post_processing_process = multiprocessing.Process(target=post_process_scans, args=(self.data_file.filename, \
				self.response_function, self.voltages, deg, self.processed_scans_buffer, self.post_process_event) ) 
		self.post_processing_process.start()
			
		# Set up timer to check whether post processing is completed
		def is_processing_over (event) :
			if not self.post_process_event.is_set() : 
				# Finish post processing
				self.post_processing_timer.Stop()
				self.post_process_event.clear()
				self.post_processing_process.join()
				# Clear
				del self.post_processing_timer, self.post_process_event, self.post_processing_process
			# Plot
			wx.CallAfter(self.display_post_processing_data, button)

		timer_id = wx.NewId()
		self.post_processing_timer = wx.Timer (self, timer_id)
		self.post_processing_timer.Start (5000)
		wx.EVT_TIMER (self, timer_id, is_processing_over)

	def display_fluoresence_rate (self, event=None) :
		"""
		Plot requested Ax3 cut
		"""
		self.fluoresence_rate_fig.clear()
		self.fluoresence_rate_fig.patch.set_facecolor('grey')
		
		axes = self.fluoresence_rate_fig.add_subplot(111, projection='3d', axisbg='grey')
		axes.set_aspect ('equal')
		axes.view_init (90, 90)

		axes.w_xaxis.set_pane_color((0.5,0.5,0.5))
		axes.w_yaxis.set_pane_color((0.5,0.5,0.5))
		axes.w_zaxis.set_pane_color((0.5,0.5,0.5))

		# Show relative values of the coordinates
		X = self.x_micron - self.x_micron.min()
		Y = self.y_micron - self.y_micron.min()
		Z = self.z_micron - self.z_micron.min()
		
		sc = axes.scatter (X, Y, Z, c=self.fluoresence_rate, cmap=cm.jet, picker=2, linewidths=0)		
		self.fluoresence_rate_fig.colorbar (sc, use_gridspec=True, orientation='horizontal')

		axes.set_xlabel ('ax1 (micron)')
		axes.set_ylabel ('ax2 (micron)')
		axes.set_zlabel ('ax3 (micron)')
		axes.set_title ("fluorescence rate (counts per millisecond)")

		self.fluoresence_rate_canvas.draw()			
		
	def extract_property (self, caption) :
		"""
		Open dialogue window for user to select property. The selected property is then returned
		as array, otherwise IOError exception is risen  
		"""
		######### Select which property will be extracted #####################
		property_name = [] 	# Names of the properties
		property_func = [] 	# Function that extract the values of properties
		for num in range(self.number_fitting_exponent.GetValue()) :
			property_name.append (self.time_life_distribution_button.GetLabel() + " %d" % (num+1) )
			property_func.append ( lambda FD : self.time_life_distribution_button.__get_data__(FD)[:,num])
			property_name.append (self.exp_ampl_distribution_button.GetLabel() + " %d" % (num+1) )
			property_func.append ( lambda FD : self.exp_ampl_distribution_button.__get_data__(FD)[:,num])
			
		ax_property = [self.fluoresence_rate_button, self.background_noise_distribution_button, 
			self.time_jitter_distribution_button, self.chisq_distribution_button, self.pval_distribution_button]
			
		property_name += map(lambda K : K.GetLabel(), ax_property)
		property_func += map(lambda K : K.__get_data__, ax_property)
		
		dlg = wx.SingleChoiceDialog (self,caption,"Select property", property_name)
		if dlg.ShowModal() <> wx.ID_OK : raise IOError 
		########### end selection #######################
		
		# Extract the specified property
		self.processed_scans_buffer.acquire()
		prop = np.copy( property_func[dlg.GetSelection()](self.processed_scans).flatten() )
		self.processed_scans_buffer.release()
		return prop
		
	def save_post_processed (self, event=None) :
		""" 
		<save_post_processed_button> was clicked to save fitted data. 
		The data will be pickled and can be subsequently edited in Blender.
		"""
		from itertools import izip
		
		# Saving image as a dictionary with key being coordinates 
		# and values representing parameters of that point
		points = {}
		
		# Check whether post processing data is available
		try : 
			self.processed_scans
			# fitted data is present

			try : prop = self.extract_property("Which property should be converted to color?")
			except IOError : return
			
			# Decide whether to use cut-off
			dlg = wx.MessageDialog(self, "Should the cut-off specified by Viewing Region (see the left panel) be applied to data?",
				'Cut-off', wx.YES | wx.NO | wx.YES_NO | wx.ICON_QUESTION )
			if dlg.ShowModal() == wx.ID_YES :
				min_val = self.value1min.GetValue(); max_val = self.value1max.GetValue()
				prop[ np.nonzero(prop < min_val) ] = min_val
				prop[ np.nonzero(prop > max_val) ] = max_val
			
			# Coveting selected property to material 
			materials, materials_of_points = array_to_materials(prop)
			
			# Forming the image data
			self.processed_scans_buffer.acquire()
			for V, TF, FD, MP in izip(self.voltages, self.fluoresence_rate, self.processed_scans, materials_of_points) :
				points[ tuple(V.tolist()) ] = { "fluoresence_rate" : np.asscalar(TF), "fitted_data" : tuple(FD.tolist()), "material" : MP } 
			self.processed_scans_buffer.release()
			
		except AttributeError :
			# the fitting data is not available so we can save only the fluorescence rate
			dlg = wx.MessageDialog(self, "There is no fitting data! Should we export the fluorescence rate data only?",
				'No fitted data', wx.OK | wx.YES_NO | wx.ICON_QUESTION )
			if dlg.ShowModal() == wx.ID_NO : return
			
			# Convert <self.fluoresence_rate> to colour 
			materials, materials_of_points = array_to_materials (self.fluoresence_rate)
			
			# Forming the image data
			for V, TF, MP in izip(self.voltages, self.fluoresence_rate, materials_of_points) :
				points[ tuple(V.tolist()) ] = { "fluoresence_rate" : np.asscalar(TF), "material" : MP }
		
		# Creating the dictionary to contain the post processed data
		fitted_data = {}
		
		# Saving image
		fitted_data["points"]  = points
		
		# Saving the dictionary of materials
		fitted_data["materials"] = materials
		
		# Saving number of exponents used for fitting
		fitted_data["number_exponent"] = self.number_fitting_exponent.GetValue()

		# Saving the checksum of underlying HDF5 file
		fitted_data["HDF5_check_sum"] = get_file_checksum(self.data_file.filename, self.paths)
		
		# Saving the name of the associated HDF5 file
		fitted_data["HDF5_filename"] = self.data_file.filename
		
		# Saving the conversion factors from units to microns
		fitted_data["unit_to_microns_ax1"] = 1.e-3 * self.piezo_calibrated_displacement_ax1
		fitted_data["unit_to_microns_ax2"] = 1.e-3 * self.piezo_calibrated_displacement_ax2
		fitted_data["unit_to_microns_ax3"] = 1.e-3 * self.piezo_calibrated_displacement_ax3
		
		# Time resolution in ps
		fitted_data["Resolution"] = self.Resolution
		
		# Save the centre of mass
		fitted_data["centre_of_mass"] = tuple( np.mean(self.voltages, axis=0).tolist() )
		
		# Save the pixel size
		tmp = np.diff(np.sort(self.voltages, axis=0),axis=0)
		max_value = np.iinfo(tmp.dtype).max
		tmp[ np.nonzero(tmp == 0) ] = max_value
		tmp = np.min(tmp,axis=0) # the pixel size is found
		# remove the outlier values in pixel size.
		tmp[ np.nonzero(tmp == max_value) ] = tmp.min()
		fitted_data["pixel_size"] = tuple( tmp.tolist() ) 
				
		# Choosing the file for dumping
		
		# Generate the default file name 
		import os
		default_filename = os.path.basename(self.data_file.filename)
		default_filename = os.path.splitext(default_filename)[0] + ".flim"
		
		openFileDialog = wx.FileDialog(self, "Choosing pickle file to save fits", "", default_filename,
                                       "pickle files (*.flim)|*.flim", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return	
		
		# dumping into the file
		import cPickle
		with open(openFileDialog.GetPath(), "wb") as out_file :
			cPickle.dump(fitted_data, out_file)
	
	def load_post_processed (self, event=None) :
		"""
		<load_post_processed_button> was clicked to load fitted data. 
		The data will be loaded from python pickle formate. 
		This method is closely related to the method <save_post_processed>.
		"""
		openFileDialog = wx.FileDialog(self, "Choosing pickle file to load fits", "", "",
                                       "pickle files (*.flim)|*.flim", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return	
		
		# loading data
		import cPickle
		with open(openFileDialog.GetPath(), "rb") as in_file :
			fitted_data = cPickle.load(in_file) 
		
		if fitted_data["HDF5_check_sum"] <> get_file_checksum(self.data_file.filename, self.paths) :
			dlg = wx.MessageDialog(None, "The fitted data may not correspond to measurements from loaded HDF5 data file. Should the fitted data be loaded?",
				'Data may not match', wx.YES_NO | wx.ICON_QUESTION )
			if dlg.ShowModal() == wx.ID_NO : return
		
		# Updating the number of exponents used for fitting
		self.number_fitting_exponent.SetValue( fitted_data["number_exponent"] )
		
		# Defying operators for fast extractions
		import operator
		get_fluoresence_rate 	= operator.itemgetter("fluoresence_rate")
		get_fitted_data			= operator.itemgetter("fitted_data")
		
		# Unrolling image
		try :
			extracted_voltages, extracted_fluoresence_rate, extracted_fitted_data \
				= zip( *( (key, get_fluoresence_rate(value), get_fitted_data(value)) for key,value in fitted_data["points"].iteritems()) )
		except KeyError :
			# No fitted data found, just load voltages and fluorescence rate
			extracted_voltages, extracted_fluoresence_rate \
				= zip( *( (key, get_fluoresence_rate(value)) for key,value in fitted_data["points"].iteritems()) )
			extracted_fitted_data = None
				
		self.initialize(extracted_voltages, extracted_fluoresence_rate, extracted_fitted_data)
		
	def initialize (self, new_voltages, new_fluoresence_rate, new_fitted_data) :
		"""
		Re-initialize arrays such as <self.voltages>, <self.fluoresence_rate> <self.processed_scans>
		"""
	
		# Converting list to arrays
		self.voltages 			= np.array(new_voltages) 
		self.fluoresence_rate = np.array(new_fluoresence_rate) 
		
		# Convert voltages to microns
		self.x_micron = 1.e-3 * self.piezo_calibrated_displacement_ax1 * self.voltages[:,0]
		self.y_micron = 1.e-3 * self.piezo_calibrated_displacement_ax2 * self.voltages[:,1]
		self.z_micron = 1.e-3 * self.piezo_calibrated_displacement_ax3 * self.voltages[:,2]
		
		if new_fitted_data is not None :
			new_fitted_data 			= np.array(new_fitted_data)
		
			########## Converting array <new_fitted_data> to <self.processed_scans> ##################
			try : del self.processed_scans_buffer, self.processed_scans
			except AttributeError : pass

			# Allocating memory 
			self.processed_scans_buffer = multiprocessing.Array (ctypes.c_double, new_fitted_data.size)

			# Creating numpy wrapper over the buffer
			self.processed_scans_buffer.acquire()
			self.processed_scans = np.frombuffer (self.processed_scans_buffer.get_obj())
			self.processed_scans = self.processed_scans.reshape( new_fitted_data.shape  )
			np.copyto (self.processed_scans, new_fitted_data)
			self.processed_scans_buffer.release()
		
		######### Re-plot with updated data  
		try : wx.CallAfter(self.display_post_processing_data, self.__previously_clicked_button__)
		except AttributeError : wx.CallAfter(self.display_fluoresence_rate)
		
	def extract_boundary (self, event) :
		"""
		Extract the boundary by tracing the image along rays coming out of the image's centre of mass.
		"""
		try : prop = self.extract_property ("Which property should be used to extract the boundary?\nValue 1 max of Viewing Region (see the left panel) is used as a typical value for the inside points.")
		except IOError : return
		
		# Calculate the Cartesian coordinates with respect to the centre of mass
		cartesian_voltages = self.voltages - np.mean(self.voltages, axis=0)
		
		# Going to the spherical coordinates (r, theta, phi)
		tmp 	= cartesian_voltages[:,0]**2 + cartesian_voltages[:,1]**2
		Radius	= np.sqrt(tmp + cartesian_voltages[:,2]**2)
		Theta 	= np.arctan2(cartesian_voltages[:,2], np.sqrt(tmp))
		Phi		= np.arctan2(cartesian_voltages[:,1], cartesian_voltages[:,0])
		del cartesian_voltages, tmp
		
		# Binning the values of the angular coordinates
		def round_binning (A, N_bins = 500) :
			""" <N_bins> - Number of bins """
			A -= A.min()
			if A.max() > 0 :
				A *= ( N_bins/A.max() )
				np.round(A, out=A)
		
		round_binning(Theta); round_binning(Phi)
		
		# Re-ranging the data in a form of rays starting from the origin  
		# i.e., representing the data as a dictionary with keys being the angular coordinates
		# and values are lists of tuple (the radial coordinate, the property, the position in  <self.voltages>)
		from itertools import izip
		Rays = {}
		for r, theta, phi, p, num in izip(Radius, Theta, Phi, prop, np.arange(self.voltages.shape[0])) :
			try : Rays[ (theta, phi) ].append( (r, p, num) )
			except KeyError : Rays[ (theta, phi) ] = [ (r, p, num) ]
		
		# Define the upper cut-off value. 
		# All the pixels with the value of the property above this threshold will be removed
		max_cut_off = self.value1max.GetValue()
		
		# Find list of indexes that are not the boundary
		index_to_be_removed = [] 
		for ray in Rays.values() : 
			# Sort the ray by radial coordinate
			ray.sort()
			# Um-zipping the ray 
			_, p, num = zip(*ray)
			try :
				# Find the last position where the value of property <p> is larger than cut-off
				max_cut_off_indx = np.nonzero(gaussian_filter(p, sigma=3) > max_cut_off)[0].max()
				# All points up to the found position should be marked for removal
				index_to_be_removed += num[:max_cut_off_indx]
			except ValueError : pass
		
		# Finally, get the indexes that represent the extracted boundary 
		boundary_indx = set(xrange(self.voltages.shape[0])) 
		boundary_indx.difference_update( index_to_be_removed )
		boundary_indx = np.array( sorted(boundary_indx) ) 
		
	
		# Extract the boundary
		new_voltages			= self.voltages[ boundary_indx, : ] 
		new_fluoresence_rate 	= self.fluoresence_rate[ boundary_indx ]
		new_processed_scans		= self.processed_scans[ boundary_indx, : ]
		
		self.initialize(new_voltages, new_fluoresence_rate, new_processed_scans)
	
#################################################################################################################

if __name__ == '__main__' :

	# Check weather the filename was given as a command line argument	
	import sys
	if len(sys.argv) == 2 : filename = sys.argv[1] 
	else : filename = None

	app = wx.App ()
	CViewer (filename=filename)
	app.MainLoop ()
