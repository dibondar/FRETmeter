import h5py
import wx
import multiprocessing

#####################################################################################################################
class PostProcessing :
	"""
	Parent class for all post processing tasks.
	"""

	def __init__ (self, source_filename=None, destination_filename=None) :
		"""
		This constructor opens a source file (to be processed) and a destination file
		(where the results of post-processing are to be saved). Moreover, all photos 
		and parameters will be copied from the source to destination
		"""
		if source_filename is None :
			# Open file dialogue to select the source file
			app = wx.App()
			openFileDialog = wx.FileDialog (None, "Choose the source file...", "", "",
										   "HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
			if openFileDialog.ShowModal() == wx.ID_CANCEL :
				raise RuntimeError("A source file must be selected")
			source_filename = openFileDialog.GetPath()
			del app
				
		if destination_filename is None :
			# Open file dialogue to select the destination file
			app = wx.App()
			openFileDialog = wx.FileDialog (None, "Open HDF5 file to load settings", "", "post_processed.hdf5",
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
			if openFileDialog.ShowModal() == wx.ID_CANCEL :
				raise RuntimeError("Destination file must be selected")
			destination_filename = openFileDialog.GetPath()
			del app
		
		# Opening source file
		self.SourceFile = h5py.File (source_filename, 'r')
		self.HistogramsSource = self.SourceFile["histograms"]
		self.ParametersSource = self.SourceFile["parameters"]
		
		# Opening destination file 
		self.DestinationFile = h5py.File (destination_filename, 'w', libver='latest')

		# Copying parameters
		self.DestinationFile.copy (self.ParametersSource, self.DestinationFile)
		self.ParametersDestination = self.DestinationFile["parameters"]

		try :
			# Copying photos, if they are present in the source
			self.PhotosSource = self.SourceFile["photos"]
			self.DestinationFile.copy (self.PhotosSource, self.DestinationFile)
			self.PhotosDestination = self.DestinationFile["photos"]
		except KeyError : pass
		
		# Creating histogram group
		self.HistogramsDestination = self.DestinationFile.create_group("histograms")

		# Get the list of revisions in the source file
		self.Revisions = set( int(key.split('_')[1]) for key in self.HistogramsSource.iterkeys () if "revision" in key )
		self.Revisions = sorted(self.Revisions)
	
	def SelectMultipleRevisions (self, title, selection=None) :
		"""
		Open a multiple choice dialogue box for selecting revisions and returns the list of revisions.

		if <selection> is not specified, then the dialogue will ask user to make selections.
		"""
		# Consistency check
		if len(self.Revisions) == 0 :
			raise RuntimeError("There are no revisions in the source file")
	
		if selection is None :
			# Ask use to make a choice
			app = wx.App()
			dlg = wx.MultiChoiceDialog(None, "Select multiple revisions", title, map(lambda n : "revision %d" %n, self.Revisions))
			if dlg.ShowModal() == wx.ID_OK : revisions_selected = map(lambda n : self.Revisions[n], dlg.GetSelections()) 
			else : raise RuntimeError ("User cancelled selections")
		elif selection == 'all' : 
			# Select all revisions
			revisions_selected = self.Revisions
		elif set(selection).issubset(self.Revisions) :
			# Provide <selection> is a proper subset of <self.Revisions>
			revisions_selected = selection
		else : raise ValueError("Invalid value of <selection>")
		
		return map( lambda n : "revision_%d_" % n, revisions_selected)
		
	def __del__ (self) :
		self.SourceFile.close()
		self.DestinationFile.close()
	
#####################################################################################################################
#
#	Implementation of utilities
#
#####################################################################################################################

def ExtractRevision (RunningEvent) :
	"""
	This function (run as a separate process) extracts revision histograms
	"""	
	try : self = PostProcessing()
	except RuntimeError, msg :
		RunningEvent.clear()
		print msg
		return
	
	app = wx.App()
	dlg = wx.SingleChoiceDialog (None, "Select revision", 
				"Revision to extract", map(lambda n : "revision %d" %n, self.Revisions))
	if dlg.ShowModal() == wx.ID_OK :
		RevisionLabel =  "revision_%d_" % self.Revisions[dlg.GetSelection()] 
	else :
		RunningEvent.clear(); return
	del dlg, app

	print "\nStart extracting %s..." % RevisionLabel

	for key, histogram in self.HistogramsSource.iteritems () :
		# check whether user requested to stop the process
		if not RunningEvent.is_set() : 
			print "Process was terminated"; break
		
		# Simply rename the histogram 
		if RevisionLabel in key :
			new_name = key[key.find("histogram"):]
			self.HistogramsDestination.copy (histogram, new_name)
	
	RunningEvent.clear()
	print "Finished extracting!"

#####################################################################################################################

def AddDirectlyRevisions (RunningEvent, source_filename=None, destination_filename=None, selection=None) :
	"""	
	This function, run as a separate process, merges selected revisions with histograms 
	by directly (without adjusting time zero offset) adding them up
	"""
	try : 
		self = PostProcessing(source_filename, destination_filename)
		SelectedRevisionLabels = self.SelectMultipleRevisions ("Revisions to be added directly to histograms", selection)
	except RuntimeError, msg :
		RunningEvent.clear()
		print msg
		return
	
	print "\nStart direct addition..." 

	# Retrieve the acquisition time
	Tacq = float(self.ParametersDestination["Tacq"][...])

	# Lists where the rate of fluorescence will be saved
	fluoresence_rate_values = [] # values of the rates
	fluoresence_rate_keys 	= [] # coordinates of the histograms for which rates are calculated
	
	for key, histogram in self.HistogramsSource.iteritems () :
		# check whether user requested to stop the process
		if not RunningEvent.is_set() : 
			print "Process was terminated"; break
		
		if key.find("histogram") == 0 :
			# We are processing the original histogram

			# Loading the original histogram
			orig_hist = histogram[...]

			# Number of scans
			# This variable accounts for the fact that different histograms may 
			# have different acquisition time.
			Nacq = 1
			for RLabel in SelectedRevisionLabels :
				# Adding to the original histogram 
				try : 
					orig_hist += self.HistogramsSource[RLabel + key][...]
					Nacq += 1
				except KeyError : pass
			
			# Saving the results of calculations into the destination file
			self.HistogramsDestination.create_dataset(key, data=orig_hist, compression='szip')
			
			# Saving the fluorescence rate for the current histogram
			# as a pair of histogram coordinates and the rate in counts per milliseconds
			fluoresence_rate_keys.append( tuple(map(int, key.split('_')[-3:])) )
			fluoresence_rate_values.append( orig_hist.sum() / (Nacq*Tacq) )
			
	# Saving the obtained fluorescence rate  
	# histogram's coordinates 
	self.ParametersDestination["fluoresence_rate_keys"] 	= fluoresence_rate_keys
	# the rate in counts per milliseconds
	self.ParametersDestination["fluoresence_rate_values"] 	= fluoresence_rate_values
	
	RunningEvent.clear()
	print "Finished adding!"


def FolderBatch_AddDirectlyRevisions (RunningEvent) :
	"""
	
	"""
	pass
	
#####################################################################################################################

def AddSmartlyRevisions (RunningEvent) :
	"""	
	This function, run as a separate process, merges selected revisions with histograms 
	by smartly (adjusting time zero offset) adding them up
	"""
	import scipy.fftpack as fftpack
	import numpy as np

	try : 
		self = PostProcessing()
		SelectedRevisionLabels = self.SelectMultipleRevisions ("Revisions to be added to histograms with time offset adjustments")
	except RuntimeError, msg :
		RunningEvent.clear()
		print msg
		return
	
	print "\nStart smart addition..." 
	
	# Retrieve the acquisition time
	Tacq = float(self.ParametersDestination["Tacq"][...])
	
	# Lists where the rate of fluorescence will be saved
	fluoresence_rate_values = [] # values of the rates
	fluoresence_rate_keys 	= [] # coordinates of the histograms for which rates are calculated

	for key in self.HistogramsSource.iterkeys () :
		# check whether user requested to stop the process
		if not RunningEvent.is_set() : 
			print "Process was terminated"; break
		
		if key.find("histogram") == 0 :
			# We found the original histogram
			
			# Loading the original histogram
			try : orig_hist = self.HistogramsSource[key][...]
			except IOError : 
				print "IOError in loading %s !!!" % key
				continue

			# Number of scans
			# This variable accounts for the fact that different histograms may 
			# have different acquisition time.
			Nacq = 1
				
			# Find out where non-zero part of it ends
			cut_off = np.nonzero(orig_hist > 0 )[0][-1]

			# Saving the FFT of the original histogram
			fft_orig_hist = fftpack.fft (orig_hist[:cut_off])

			# Creating the Gaussian for filtering
			gaussian = np.exp( -5.*fftpack.fftfreq(fft_orig_hist.size)**2 )

			for RLab in SelectedRevisionLabels :
				# Adding to the original histogram by finding proper time zero offset
				try :
					rev_hist = self.HistogramsSource[RLab + key][...]
					rev_hist = rev_hist[:cut_off]

					# Using the FFT circular convolution theorem
					# calculate the convolution function
					convol = fftpack.fft(rev_hist[::-1])
					convol *= fft_orig_hist
					convol *= gaussian # Smoothing the convolution by applying the Gaussian filter
					convol = fftpack.ifft(convol, overwrite_x=True)

					# Shifting the revision histogram such that the time zero matches
					# and add it to the original histogram
					orig_hist[:cut_off] += np.roll(rev_hist, np.argmax(convol))
					
					Nacq += 1
				except (KeyError, IOError) : pass
			
			# Saving the results of calculations into the destination file
			self.HistogramsDestination.create_dataset(key, data=orig_hist, compression='szip')
			
			# Saving the fluorescence rate for the current histogram
			# histogram's coordinates 
			fluoresence_rate_keys.append( tuple(map(int, key.split('_')[-3:])) )
			# the rate in counts per milliseconds
			fluoresence_rate_values.append( orig_hist.sum() / (Nacq*Tacq) )
			
	# Saving the obtained fluorescence rate  
	self.ParametersDestination["fluoresence_rate_keys"] 	= fluoresence_rate_keys
	self.ParametersDestination["fluoresence_rate_values"] 	= fluoresence_rate_values

	RunningEvent.clear()
	print "Finished adding!"

#####################################################################################################################

def ExtractAxisCut (RunningEvent, axis) :
	"""	
	This function, run as a separate process, extracts histograms with certain value of the coordinate at the given axis  
	(axis X, Y, Z corresponds to axis=-3,-2,-1). 
	"""
	try : self = PostProcessing()
	except RuntimeError, msg :
		RunningEvent.clear()
		print msg
		return

	# Extract all the values along the chosen axis
	axis_values = set( int(key.split('_')[axis]) for key in self.HistogramsSource.iterkeys() )
	axis_values = sorted(axis_values) 
	
	# Consistency check
	if len(axis_values) == 0 : raise RuntimeError("There are no revisions in the source file")

	# Select values to be extracted	
	app = wx.App()
	dlg = wx.MultiChoiceDialog(None, "Select multiple values", "Chose coordinate values to be extracted...", map(str, axis_values))
	if dlg.ShowModal() == wx.ID_OK : values_to_be_extracted = [ axis_values[n] for n in dlg.GetSelections() ]
	else :
		RunningEvent.clear(); raise RuntimeError ("User cancelled selections")
	del app

	print "\nStart extracting coordinate cut..."

	for key, histogram in self.HistogramsSource.iteritems () :
		# check whether user requested to stop the process
		if not RunningEvent.is_set() : 
			print "Process was terminated"; break
	
		# check whether this histogram belongs to the desired cut
		if int(key.split('_')[axis]) in values_to_be_extracted :
			self.HistogramsDestination.copy (histogram, key)		
		
	RunningEvent.clear()
	print "Finished extracting coordinate!"

#####################################################################################################################

#####################################################################################################################
#
#	GUI for post processing utilities
#
#####################################################################################################################

class PostProcessingUtilities(wx.Frame) :

	def __init__ (self) :
		wx.Frame.__init__ (self, None, title="Post processing utilities")
		self.ConstructGUI ()
		self.Show ()
		wx.EVT_CLOSE (self, self.OnClose)

	def OnClose (self, event) :
		# Stopping all processes
		for key in dir(self) :
			button = getattr (self, key)
			if isinstance(button, (wx.Button)) :
				try : button.__RunningEvent__.clear()
				except AttributeError : pass
		self.Destroy()

	def ConstructGUI (self) :
		"""
		Make GUI
		"""
		panel = wx.Panel(self)
		sizer = wx.BoxSizer (wx.VERTICAL)

		# utility for extracting rev histograms
		self.ExtractRevision_Button = wx.Button (panel)
		self.ExtractRevision_Button.__StartLabel__ 	= "Extract revision..."
		self.ExtractRevision_Button.__StopLabel__	= "STOP Extracting revision"
		self.ExtractRevision_Button.SetLabel (self.ExtractRevision_Button.__StartLabel__)
		self.ExtractRevision_Button.__exec__ = ExtractRevision
		self.ExtractRevision_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.ExtractRevision_Button, flag=wx.EXPAND, border=5)

		# utility for directly adding up histograms
		self.AddDirectlyRevision_Button = wx.Button (panel)
		self.AddDirectlyRevision_Button.__StartLabel__ 	= "Add directly revisions..."
		self.AddDirectlyRevision_Button.__StopLabel__	= "STOP adding revisions"
		self.AddDirectlyRevision_Button.SetLabel (self.AddDirectlyRevision_Button.__StartLabel__)
		self.AddDirectlyRevision_Button.__exec__ = AddDirectlyRevisions
		self.AddDirectlyRevision_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.AddDirectlyRevision_Button, flag=wx.EXPAND, border=5)

		# utility for smartly adding up histograms
		self.AddSmartlyRevision_Button = wx.Button (panel)
		self.AddSmartlyRevision_Button.__StartLabel__ 	= "Add smartly revisions..."
		self.AddSmartlyRevision_Button.__StopLabel__	= "STOP adding revisions"
		self.AddSmartlyRevision_Button.SetLabel (self.AddSmartlyRevision_Button.__StartLabel__)
		self.AddSmartlyRevision_Button.__exec__ = AddSmartlyRevisions
		self.AddSmartlyRevision_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.AddSmartlyRevision_Button, flag=wx.EXPAND, border=5)
		
		from functools import partial
		
		# utility to extract X cut
		self.ExtractXCut_Button = wx.Button (panel)
		self.ExtractXCut_Button.__StartLabel__  = "Extract X cut..."
		self.ExtractXCut_Button.__StopLabel__ 	= "STOP extracting X cut"
		self.ExtractXCut_Button.SetLabel (self.ExtractXCut_Button.__StartLabel__)
		self.ExtractXCut_Button.__exec__ = partial(ExtractAxisCut, axis=-3)
		self.ExtractXCut_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.ExtractXCut_Button, flag=wx.EXPAND, border=5)

		# utility to extract Y cut
		self.ExtractYCut_Button = wx.Button (panel)
		self.ExtractYCut_Button.__StartLabel__  = "Extract Y cut..."
		self.ExtractYCut_Button.__StopLabel__ 	= "STOP extracting Y cut"
		self.ExtractYCut_Button.SetLabel (self.ExtractYCut_Button.__StartLabel__)
		self.ExtractYCut_Button.__exec__ = partial(ExtractAxisCut, axis=-2)
		self.ExtractYCut_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.ExtractYCut_Button, flag=wx.EXPAND, border=5)

		# utility to extract Z cut
		self.ExtractZCut_Button = wx.Button (panel)
		self.ExtractZCut_Button.__StartLabel__  = "Extract Z cut..."
		self.ExtractZCut_Button.__StopLabel__ 	= "STOP extracting Z cut"
		self.ExtractZCut_Button.SetLabel (self.ExtractZCut_Button.__StartLabel__)
		self.ExtractZCut_Button.__exec__ = partial(ExtractAxisCut, axis=-1)
		self.ExtractZCut_Button.Bind (wx.EVT_LEFT_DOWN, self.OnUtility)

		sizer.Add (self.ExtractZCut_Button, flag=wx.EXPAND, border=5)

		panel.SetSizer (sizer)

	def OnUtility (self, event) :

		try : button = event.GetEventObject()
		except AttributeError :  button = event

		if button.GetLabel() == button.__StartLabel__ :
			# Starting utility
			button.__RunningEvent__ = multiprocessing.Event()
			button.__RunningEvent__.set()

			multiprocessing.Process( target=button.__exec__, args=(button.__RunningEvent__,) ).start()	
				
			# Starting timer to check when the utility is over
			TIMER_ID = wx.NewId()
			button.__timer__ = wx.Timer (button, TIMER_ID)
			button.__timer__.Start (3000) # check every 3 second

			def IsUtilityFinished (E) :
				if not button.__RunningEvent__.is_set() : self.OnUtility (button)

			wx.EVT_TIMER (button, TIMER_ID, IsUtilityFinished)

			# Changing button's label
			button.SetLabel( button.__StopLabel__ )

		elif button.GetLabel() == button.__StopLabel__ :
			# Stopping the utility and the timer
			button.__timer__.Stop()
			button.__RunningEvent__.clear()
			del button.__RunningEvent__, button.__timer__
			button.SetLabel( button.__StartLabel__ )

		else : raise NotImplementedError ("Button's label is cont recognized") 

#################################################################################################################

if __name__ == '__main__' :
	app = wx.App ()
	PostProcessingUtilities ()
	app.MainLoop ()
