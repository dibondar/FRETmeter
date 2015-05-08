import wx
import multiprocessing 
import itertools
import time
import operator
import numpy as np
import h5py
import ctypes
# Real time plotting
import visvis
from os.path import isfile
app = visvis.use('wx')

####################################################################################################

# Commands for intra process communication
CMD_EXIT = 0

CMD_MEASURE_HISTOGRAM 		= 1
CMD_INITIATE_PICO_HARP 		= 2
CMD_GET_COUNT_RATES_PICO_HARP	= 3
CMD_GET_RESOLUTION_PICO_HARP	= 4
CMD_SET_OFFSET_PICO_HARP	= 5
CMD_SET_CFDZERO_PICO_HARP	= 6
CMD_SET_CFDLEVEL_PICO_HARP	= 7
CMD_SET_DIVIDER_PICO_HARP	= 8
CMD_SET_RANGE_PICO_HARP		= 9
CMD_SET_TACQ_PICO_HARP		= 10
	
CMD_OPEN_SHUTTER = 12
CMD_CLOSE_SHUTTER = 13

CMD_TAKE_PHOTO			= 16
CMD_TAKE_PHOTO_ACCURATELY	= 17

# scanning
CMD_TEST_SCANNING 	= 18
CMD_SCAN		= 19
CMD_TAKE_SCANNING_VIDEO	= 20
CMD_ADAPTIVE_SCAN	= 21
CMD_BOUNDARY_SCAN	= 22
CMD_BLEACHING_STUDY	= 23

# it is important for debugging that the following two constants are not numeric 
RETURN_SUCCESS 	= []
RETURN_FAIL	= None

####################################### Hardware constants ##############################################

# Thorlabs camera parameters
CAMERA_IMG_WIDTH = 1280
CAMERA_IMG_HIGHT = 1024

######################################### Utilities ###########################################################

def camera_img_shape (camera_img_buffer) :
	"""
	Convert a multiprocess Array <camera_img_buffe> storing the image into a numpy array 
	"""
	return np.frombuffer(camera_img_buffer.get_obj(), dtype=np.uint8).reshape( (CAMERA_IMG_HIGHT, CAMERA_IMG_WIDTH) )

# PicoHarp histogram
HISTCHAN  =  65536 # number of histogram channel -- a constant from Phdefin.h

def histogram_shape (histogram_buffer) :
	"""
	Convert a multiprocess Array <histogram_buffer>, storing a Pico Harp histogram, into a numpy array
	"""
	return np.frombuffer (histogram_buffer.get_obj(), dtype=np.uint)
	

def get_scanning_polygon (ScanParameters) :
	"""
	This function generates a list of tuples (vol_ax1, vol_ax2, vol_ax3) that constitute a polygon to be used for scanning
	"""
	# Computational geometry
	from shapely.geometry import Point, Polygon
	
	# Conversion factor
	mm_per_V = 1e-6*ScanParameters["unit_size"] 
		
	# Converting points from mm to units
	points_ax1, points_ax2 = zip(*ScanParameters["points"])
	points_ax1 = np.asarray( np.round(np.array(points_ax1) / mm_per_V), dtype=np.int)
	points_ax2 = np.asarray( np.round(np.array(points_ax2) / mm_per_V), dtype=np.int)
		
	# Generating polygon out of specified points
	polygon = Polygon( zip(points_ax1, points_ax2) )
	del points_ax1, points_ax2
		
	# Find coordinate bounds
	min_ax1, min_ax2, max_ax1, max_ax2 = polygon.bounds
		
	# Generate intervals
	Ax1 = np.linspace(min_ax1, max_ax1, ScanParameters["number_steps_ax1"])
	Ax1 = np.asarray(np.round(Ax1), dtype=np.int)
	Ax2 = np.linspace(min_ax2, max_ax2, ScanParameters["number_steps_ax2"])
	Ax2 = np.asarray(np.round(Ax2), dtype=np.int)
	Ax1, Ax2 = np.meshgrid (Ax1, Ax2)
	Ax1 = Ax1.flatten(); Ax2 = Ax2.flatten()

	# Exclude points outside the polygon
	indx = np.nonzero( map(polygon.contains, map(Point, zip(Ax1, Ax2))) ) 
	Ax1 = Ax1[indx]; Ax2 = Ax2[indx]

	# Adding layers
	Ax1_copy = np.copy(Ax1); Ax2_copy = np.copy(Ax2) 
	Ax1 = []; Ax2 = []; Ax3 = []
	
	LayersAx3 = np.linspace (ScanParameters["piezo_min_volt_ax3"], 
						ScanParameters["piezo_max_volt_ax3"], ScanParameters["number_steps_ax3"])
	LayersAx3 = np.asarray(np.round(LayersAx3), dtype=np.int)
	for ax3_layer in LayersAx3 :
		Ax1 = np.append (Ax1, Ax1_copy)
		Ax2 = np.append (Ax2, Ax2_copy)
		Ax3 = np.append (Ax3, ax3_layer*np.ones_like(Ax1_copy))
	
	voltages = list(set( zip(Ax1, Ax2, Ax3) ))
	
	# Scanning ordering
	slow_varying_axis = ScanParameters["slow_varying_axis"]
	if slow_varying_axis == 3 :
		# Randomized ordering
		import random
		random.shuffle (voltages)
	else :
		# Order by slow varying axes
		voltages = sorted( voltages, key=operator.itemgetter(0) )
		voltages = sorted( voltages, key=operator.itemgetter(1) )
		voltages = sorted( voltages, key=operator.itemgetter(2) )
		voltages = sorted( voltages, key=operator.itemgetter(slow_varying_axis) )

	return voltages
		
####################################################################################################

def control_pico_harp (input_pico_harp_queue, output_pico_harp_queue, histogram_buffer) :
	"""
	This function runs as a separate process and controls the Pico Harp 3000.
	All communications are done through the input and output queues.
	<histogram_buffer> a buffer where the histogram is saved.
	"""
	from ctypes import c_int

	lib = ctypes.WinDLL ("C:/Windows/System32/phlib.dll")
	
	# String buffer
	str_p =  ctypes.create_string_buffer (500)

	# Extract the Pico Harp library version
	lib.PH_GetLibraryVersion (str_p)
	pico_harp_lib_version_dev = "2.3"
	if str_p.value <> pico_harp_lib_version_dev :
		print "Warning: The current code was developed for PicoHarp Library version %s. The version of the currently installed library is %s" \
				% (pico_harp_lib_version_dev, str_p.value)
	
	# Constants from Phdefin.h
	MODE_HIST = c_int(0)
	FLAG_OVERFLOW = 0x0040
	
	# Pico harp device ID to be utilized
	DEV_INDX = c_int(0)

	def check_for_warnings () :
		"""
		This function checks for wornings, and if there are any it prints them.
		"""
		warnings = lib.PH_GetWarnings (DEV_INDX)
		if warnings < 0 :
			print "\nError: PH_GetWarnings failed\n"; raise RuntimeError
		elif warnings > 0 :
			# Get the warning mesage
			if lib.PH_GetWarningsText (DEV_INDX, str_p, warnings) < 0 :
				print "\nError: PH_GetWarningsText failed\n"; raise RuntimeError
			print str_p.value


	###### Checking for the commands from other processes
	while True :
		# Getting the full command with arguments
		full_command = input_pico_harp_queue.get ()

		# Extracting just the name of the command
		try : command = full_command[0]
		except TypeError : command = full_command

		if   command == CMD_EXIT : break # exiting the process
		elif command == CMD_MEASURE_HISTOGRAM :
		##### Histogram measurement is requested. <RETURN_SUCCESS> or <RETURN_FAIL> will be retuned 
		##### and the histogram will be saved in <histogram_buffer> ####
			try :
				# Saving the acquisition time, if provided
				try : Tacq = full_command[1]
				except TypeError : pass

				if lib.PH_ClearHistMem (DEV_INDX,  c_int(0)) < 0 : # always use Block 0 if not Routing
					print "\nError: PH_ClearHistMem failed\n"; raise RuntimeError

				if lib.PH_StartMeas (DEV_INDX, Tacq) < 0 :
					print "\nError: PH_StartMeas failed\n"; raise RuntimeError

				# Measuring for <Tacq> ms
				ctcdone = 0
				while ctcdone == 0 : ctcdone = lib.PH_CTCStatus (DEV_INDX)

				if lib.PH_StopMeas (DEV_INDX) < 0 :
					print "\nError: PH_StopMeas failed\n"; raise RuntimeError
				
				# Retrieving the histogram
				histogram_buffer.acquire()
				ret = lib.PH_GetBlock (DEV_INDX, histogram_buffer.get_obj(), c_int(0)) 
				histogram_buffer.release()
				if ret < 0 :
					print "\nError: PH_GetBlock failed\n"; raise RuntimeError
	
				flags = lib.PH_GetFlags (DEV_INDX)
				if flags < 0 : 
					print "\nError: PH_GetFlags failed\n"; raise RuntimeError
				
				if flags & FLAG_OVERFLOW : print "\n\t!!!!!!Overflow!!!!!!!\n"

				# Returning the histogram through the output queue
				output_pico_harp_queue.put (RETURN_SUCCESS)

			except RuntimeError : output_pico_harp_queue.put (RETURN_FAIL)

		elif command == CMD_INITIATE_PICO_HARP :
		###### Initialization of Pico Harp is requested. <Resolution> will be returned ######
		
			# Extracting the arguments (i.e., settings for the measurement) from the full command
			ScanParameters 	= full_command[1]

			Offset		= c_int(ScanParameters["Offset"])
			CFDZeroX0	= c_int(ScanParameters["CFDZeroX0"])
			CFDLevel0	= c_int(ScanParameters["CFDLevel0"])
			CFDZeroX1 	= c_int(ScanParameters["CFDZeroX1"])
			CFDLevel1 	= c_int(ScanParameters["CFDLevel1"])
			SyncDiv 	= c_int(ScanParameters["SyncDiv"])
			Range 		= c_int(ScanParameters["Range"])
			Tacq 		= c_int(ScanParameters["Tacq"])
			
			# Closing Pico Harp, just in case if it was already open 
			lib.PH_CloseDevice (DEV_INDX)
	
			try :
				# Open the first Pico Harp
				if lib.PH_OpenDevice (DEV_INDX, str_p) < 0 :
					print "\nError: The first Pico Harp could not be opened\n"; raise RuntimeError

				# Initializing the Pico Harp device with device ID <DEV_INDEX> in the Histogram mode 
				if lib.PH_Initialize (DEV_INDX, MODE_HIST) < 0 :
					print "\nError: The first Pico Harp device could not be initiated in the Histogram mode\n"; raise RuntimeError
		
				if lib.PH_Calibrate (DEV_INDX) < 0 :
					print "\nError: Pico Harp calibration failed\n"; raise RuntimeError
		
				if lib.PH_SetSyncDiv (DEV_INDX, SyncDiv) < 0 :
					print "\nError: PH_SetSyncDiv failed\n"; raise RuntimeError

				if lib.PH_SetCFDLevel (DEV_INDX, c_int(0), CFDLevel0) < 0 :
					print "\nError: PH_SetCFDLevel 0 failed\n"; raise RuntimeError

				if lib.PH_SetCFDLevel (DEV_INDX, c_int(1), CFDLevel1) < 0 :
					print "\nError: PH_SetCFDLevel 1 failed\n"; raise RuntimeError

				if lib.PH_SetCFDZeroCross (DEV_INDX, c_int(0), CFDZeroX0) < 0 :
					print "\nError: PH_SetCFDZeroCross 0 failed\n"; raise RuntimeError

				if lib.PH_SetCFDZeroCross (DEV_INDX, c_int(1), CFDZeroX1) < 0 :
					print "\nError: PH_SetCFDZeroCross 1 failed\n"; raise RuntimeError

				if lib.PH_SetRange (DEV_INDX, Range) < 0 :
					print "\nError: PH_SetRange failed\n"; raise RuntimeError

				if lib.PH_SetOffset (DEV_INDX, Offset) < 0 :
					print "\nError: PH_SetOffset failed\n"; raise RuntimeError

				if lib.PH_SetStopOverflow (DEV_INDX, c_int(1), c_int(65535)) < 0 :
					print "\nError: PH_SetStopOverflow failed\n"; raise RuntimeError

				Resolution = lib.PH_GetResolution (DEV_INDX)
				if Resolution < 0 :
					print "\nError: PH_GetResolution failed\n"; raise RuntimeError

				# Note: after Init or SetSyncDiv you must allow 100 ms for valid new count rate readings
				time.sleep (0.2)
		
				Countrate0 = lib.PH_GetCountRate (DEV_INDX, c_int(0))
				Countrate1 = lib.PH_GetCountRate (DEV_INDX, c_int(1))

				print "Pico Harp Serial# %s is successfully initialized and ready for measurements" % str_p.value 
				print "Parameters:\nResolution=%1dps\nInput0 countrate = %1d/s\nInput1 countrate = %1d/s\n" % (Resolution, Countrate0, Countrate1)
				
				check_for_warnings ()

				### Initialization of Pico Harp is completed; <Resolution> is returned
				output_pico_harp_queue.put (Resolution)

			except RuntimeError : output_pico_harp_queue.put (RETURN_FAIL)

		elif command == CMD_GET_COUNT_RATES_PICO_HARP :
			# Get count rates
			Countrate0 = lib.PH_GetCountRate (DEV_INDX, c_int(0))
			Countrate1 = lib.PH_GetCountRate (DEV_INDX, c_int(1))
			if Countrate0 < 0 or Countrate1 < 0 : 
				print "Error in getting the Pico Harp count rate\n"
				output_pico_harp_queue.put (RETURN_FAIL)
			else : output_pico_harp_queue.put( (Countrate0, Countrate1) )

		elif command == CMD_GET_RESOLUTION_PICO_HARP :
			# Get time resolution of Pico Harp
			Resolution = lib.PH_GetResolution (DEV_INDX)
			if Resolution < 0 :
				print "Error in getting the Pico Harp time resolution\n"
				output_pico_harp_queue.put (RETURN_FAIL)
			else : output_pico_harp_queue.put(Resolution)

		# Interactive settings of PicoHarp settings (no returns)
		elif command == CMD_SET_OFFSET_PICO_HARP	: lib.PH_SetOffset (DEV_INDX, full_command[1])
		elif command == CMD_SET_CFDZERO_PICO_HARP	: lib.PH_SetCFDZeroCross (DEV_INDX, full_command[1], full_command[2])
		elif command == CMD_SET_CFDLEVEL_PICO_HARP	: lib.PH_SetCFDLevel (DEV_INDX, full_command[1], full_command[2])
		elif command == CMD_SET_DIVIDER_PICO_HARP	: lib.PH_SetSyncDiv (DEV_INDX,  full_command[1])
		elif command == CMD_SET_RANGE_PICO_HARP		: lib.PH_SetRange (DEV_INDX,  full_command[1])
		elif command == CMD_SET_TACQ_PICO_HARP		: Tacq = full_command[1]
		
		else : 
			print "PicoHarp process: Unrecognized command"; output_pico_harp_queue.put (RETURN_FAIL) 

	# Closing Pico Harp
	lib.PH_CloseDevice (DEV_INDX)
		
# The following function is a substitute 
def control_pico_harp_SUBSTITUTE (input_pico_harp_queue, output_pico_harp_queue, histogram_buffer) :
	"""
	!!!!!NOTE: THIS FUNCTION USES HYDRA HARP!!!!!
	
	This function runs as a separate process and controls the Pico Harp 3000.
	All communications are done through the input and output queues.
	<histogram_buffer> a buffer where the histogram is saved.
	"""
	from ctypes import c_int

	lib = ctypes.WinDLL ("C:/Windows/System32/hhlib.dll")
	
	# String buffer
	str_p =  ctypes.create_string_buffer (500)

	# Extract the Pico Harp library version
	lib.HH_GetLibraryVersion (str_p)
	pico_harp_lib_version_dev = "2.1"
	if str_p.value <> pico_harp_lib_version_dev :
		print "Warning: The current code was developed for PicoHarp Library version %s. The version of the currently installed library is %s" \
				% (pico_harp_lib_version_dev, str_p.value)
	
	# Constants from Phdefin.h
	MODE_HIST = c_int(0)
	FLAG_OVERFLOW = 0x0040
	
	# Pico harp device ID to be utilized
	DEV_INDX = c_int(0)

	def check_for_warnings () :
		"""
		This function checks for wornings, and if there are any it prints them.
		"""
		warnings = ctypes.c_int()
		if lib.HH_GetWarnings (DEV_INDX, ctypes.byref(warnings)) < 0 :
			print "\nError: HH_GetWarnings failed\n"; raise RuntimeError
		elif warnings.value > 0 :
			# Get the warning mesage
			if lib.HH_GetWarningsText (DEV_INDX, str_p, warnings) < 0 :
				print "\nError: HH_GetWarningsText failed\n"; raise RuntimeError
			print str_p.value


	###### Checking for the commands from other processes
	while True :
		# Getting the full command with arguments
		full_command = input_pico_harp_queue.get ()

		# Extracting just the name of the command
		try : command = full_command[0]
		except TypeError : command = full_command

		if   command == CMD_EXIT : break # exiting the process
		elif command == CMD_MEASURE_HISTOGRAM :
		##### Histogram measurement is requested. <RETURN_SUCCESS> or <RETURN_FAIL> will be retuned 
		##### and the histogram will be saved in <histogram_buffer> ####
			try :
				# Saving the acquisition time, if provided
				try : Tacq = full_command[1]
				except TypeError : pass

				if lib.HH_ClearHistMem (DEV_INDX) < 0 : # always use Block 0 if not Routing
					print "\nError: HH_ClearHistMem failed\n"; raise RuntimeError

				if lib.HH_StartMeas (DEV_INDX, Tacq) < 0 :
					print "\nError: HH_StartMeas failed\n"; raise RuntimeError

				# Measuring for <Tacq> ms
				ctcdone = ctypes.c_int(0)
				while ctcdone.value == 0 : 
					lib.HH_CTCStatus (DEV_INDX, ctypes.byref(ctcdone))

				if lib.HH_StopMeas (DEV_INDX) < 0 :
					print "\nError: HH_StopMeas failed\n"; raise RuntimeError
				
				# Retrieving the histogram
				histogram_buffer.acquire()
				ret = lib.HH_GetHistogram (DEV_INDX, histogram_buffer.get_obj(), c_int(0), c_int(0)) 
				histogram_buffer.release()
				if ret < 0 :
					print "\nError: HH_GetBlock failed\n"; raise RuntimeError
	
				flags = ctypes.c_int()
				if lib.HH_GetFlags (DEV_INDX, ctypes.byref(flags)) < 0 : 
					print "\nError: HH_GetFlags failed\n"; raise RuntimeError
				
				if flags.value & FLAG_OVERFLOW : print "\n\t!!!!!!Overflow!!!!!!!\n"

				# Returning the histogram through the output queue
				output_pico_harp_queue.put (RETURN_SUCCESS)

			except RuntimeError : output_pico_harp_queue.put (RETURN_FAIL)

		elif command == CMD_INITIATE_PICO_HARP :
		###### Initialization of Pico Harp is requested. <Resolution> will be returned ######
		
			# Extracting the arguments (i.e., settings for the measurement) from the full command
			ScanParameters 	= full_command[1]

			Offset		= c_int(ScanParameters["Offset"])
			CFDZeroX0	= c_int(ScanParameters["CFDZeroX0"])
			CFDLevel0	= c_int(ScanParameters["CFDLevel0"])
			CFDZeroX1 	= c_int(ScanParameters["CFDZeroX1"])
			CFDLevel1 	= c_int(ScanParameters["CFDLevel1"])
			SyncDiv 	= c_int(ScanParameters["SyncDiv"])
			Range 		= c_int(ScanParameters["Range"])
			Tacq 		= c_int(ScanParameters["Tacq"])
			
			# Closing Pico Harp, just in case if it was already open 
			lib.HH_CloseDevice (DEV_INDX)
	
			try :
				# Open the first Pico Harp
				if lib.HH_OpenDevice (DEV_INDX, str_p) < 0 :
					print "\nError: The first Pico Harp could not be opened\n"; raise RuntimeError

				# Initializing the Pico Harp device with device ID <DEV_INDEX> in the Histogram mode 
				if lib.HH_Initialize (DEV_INDX, MODE_HIST,0) < 0 :
					print "\nError: The first Pico Harp device could not be initiated in the Histogram mode\n"; raise RuntimeError
		
				if lib.HH_Calibrate (DEV_INDX) < 0 :
					print "\nError: Pico Harp calibration failed\n"; raise RuntimeError
		
				if lib.HH_SetSyncDiv (DEV_INDX, SyncDiv) < 0 :
					print "\nError: HH_SetSyncDiv failed\n"; raise RuntimeError

				if lib.HH_SetSyncCFD (DEV_INDX, CFDLevel0, CFDZeroX0) < 0 :
					print "\nError: HH_SetSyncCFD failed\n"; raise RuntimeError

				if lib.HH_SetInputCFD (DEV_INDX, c_int(0), CFDLevel1, CFDZeroX1) < 0 :
					print "\nError: HH_SetInputCFD (1 channel) failed\n"; raise RuntimeError

				if lib.HH_SetBinning (DEV_INDX, Range) < 0 :
					print "\nError: HH_SetBinning failed\n"; raise RuntimeError

				if lib.HH_SetOffset (DEV_INDX, Offset) < 0 :
					print "\nError: HH_SetOffset failed\n"; raise RuntimeError

				if lib.HH_SetStopOverflow (DEV_INDX, c_int(1), c_int(65535)) < 0 :
					print "\nError: HH_SetStopOverflow failed\n"; raise RuntimeError

				Resolution = ctypes.c_double()
				if lib.HH_GetResolution (DEV_INDX, ctypes.byref(Resolution)) < 0 :
					print "\nError: HH_GetResolution failed\n"; raise RuntimeError
				Resolution = Resolution.value
					
				# Note: after Init or SetSyncDiv you must allow 400 ms for valid new count rate readings
				time.sleep (0.4)
		
				Countrate0 = ctypes.c_int()
				lib.HH_GetSyncRate (DEV_INDX, ctypes.byref(Countrate0))
				Countrate0 = Countrate0.value
				
				Countrate1 = ctypes.c_int()
				lib.HH_GetCountRate (DEV_INDX, c_int(0), ctypes.byref(Countrate1))
				Countrate1 = Countrate1.value
				
				print "Pico Harp Serial# %s is successfully initialized and ready for measurements" % str_p.value 
				print "Parameters:\nResolution=%1dps\nInput0 countrate = %1d/s\nInput1 countrate = %1d/s\n" % (Resolution, Countrate0, Countrate1)
				
				check_for_warnings ()

				### Initialization of Pico Harp is completed; <Resolution> is returned
				output_pico_harp_queue.put (Resolution)

			except RuntimeError : output_pico_harp_queue.put (RETURN_FAIL)

		elif command == CMD_GET_COUNT_RATES_PICO_HARP :
			# Get count rates
			Countrate0 = ctypes.c_int()
			lib.HH_GetSyncRate (DEV_INDX, ctypes.byref(Countrate0))
			Countrate0 = Countrate0.value
				
			Countrate1 = ctypes.c_int()
			lib.HH_GetCountRate (DEV_INDX, c_int(0), ctypes.byref(Countrate1))
			Countrate1 = Countrate1.value
				
			if Countrate0 < 0 or Countrate1 < 0 : 
				print "Error in getting the Pico Harp count rate\n"
				output_pico_harp_queue.put (RETURN_FAIL)
			else : output_pico_harp_queue.put( (Countrate0, Countrate1) )

		elif command == CMD_GET_RESOLUTION_PICO_HARP :
			# Get time resolution of Pico Harp
			Resolution = ctypes.c_double()
			if lib.HH_GetResolution (DEV_INDX, ctypes.byref(Resolution)) < 0 :
				print "\nError: HH_GetResolution failed\n"
				output_pico_harp_queue.put (RETURN_FAIL)				
			else : output_pico_harp_queue.put(Resolution.value)

		# Interactive settings of PicoHarp settings (no returns)
		elif command == CMD_SET_OFFSET_PICO_HARP	: lib.HH_SetOffset (DEV_INDX, full_command[1])
		elif command == CMD_SET_CFDZERO_PICO_HARP	: pass #lib.HH_SetCFDZeroCross (DEV_INDX, full_command[1], full_command[2])
		elif command == CMD_SET_CFDLEVEL_PICO_HARP	: pass #lib.HH_SetCFDLevel (DEV_INDX, full_command[1], full_command[2])
		elif command == CMD_SET_DIVIDER_PICO_HARP	: lib.HH_SetSyncDiv (DEV_INDX,  full_command[1])
		elif command == CMD_SET_RANGE_PICO_HARP		: lib.HH_SetBinning (DEV_INDX,  full_command[1])
		elif command == CMD_SET_TACQ_PICO_HARP		: Tacq = full_command[1]
		
		else : 
			print "PicoHarp process: Unrecognized command"; output_pico_harp_queue.put (RETURN_FAIL) 

	# Closing Pico Harp
	lib.HH_CloseDevice (DEV_INDX)
		

####################################################################################################

def control_camera (input_camera_queue, output_camera_queue, camera_img_buffer) :
	"""
	This function runs a separate process to control Thorlabs camera.
	All communications are done through queues.
	<camera_img_buffer> buffer where the image is stored.
	"""
	# Create a numpy proxy for <camera_img_buffer>
	camera_img_buffer.acquire ()
	img_buffer = camera_img_shape (camera_img_buffer)
	camera_img_buffer.release ()
	
	################## Initialize the camera ###################
	import ctypes.wintypes
	
	# some constants
	IS_SUCCESS 			= 0
	IS_CAPTURE_RUNNING 	= 140
	
	# Load driver
	lib = ctypes.cdll.LoadLibrary("uc480")
	
	num_cams = ctypes.c_int()
	if lib.is_GetNumberOfCameras( ctypes.byref(num_cams) ) != IS_SUCCESS :
		raise RuntimeError("Error in is_GetNumberOfCameras")

	if num_cams.value > 1 :
		print "More than one camera is detect. The first one will be used."
	
	# Camera handle
	ph_cam = ctypes.wintypes.HDC(0)

	# Window handle
	h_wnd = ctypes.wintypes.HWND(0)

	# Initialize camera
	if lib.is_InitCamera( ctypes.byref(ph_cam), h_wnd ) != IS_SUCCESS :
		raise RuntimeError("Error in is_InitCamera")
	
	# ID of the memory for image storage
	pid = ctypes.c_int()
	# pointer to buffer
	pc_img_mem = img_buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_char))

	if lib.is_SetAllocatedImageMem (ph_cam,  CAMERA_IMG_WIDTH, CAMERA_IMG_HIGHT, 
		8*img_buffer.dtype.itemsize, pc_img_mem, ctypes.byref(pid) 
	) != IS_SUCCESS :
		raise RuntimeError("Error in is_SetAllocatedImageMem")

	# Set colour
	IS_CM_MONO8 = 6
	if lib.is_SetColorMode (ph_cam,  IS_CM_MONO8) != IS_SUCCESS :
		raise RuntimeError("Error in is_SetColorMode")
		
	# Set active memory
	if lib.is_SetImageMem(ph_cam, pc_img_mem, pid) != IS_SUCCESS :
		raise RuntimeError("Error in is_SetImageMem")	
		
	##################  Camera initiated ###################

	# Draw circule of 10 micorns radius 
	sqr = (np.arange(img_buffer.shape[0])[:,np.newaxis] - img_buffer.shape[0]/2.)**2 \
			+ (np.arange(img_buffer.shape[1])[np.newaxis,:] - img_buffer.shape[1]/2.)**2
	radius = 8.*img_buffer.shape[0]/190. # Adjust as needed
	circular_cut = np.nonzero( ( (radius-1)**2 < sqr)&(sqr < (radius+1)**2 )  )
	del sqr, radius

	def AquireImg () :
		"""
		Acquire image
		"""
		camera_img_buffer.acquire ()
		
		if lib.is_FreezeVideo(ph_cam, 0x0001) in [IS_SUCCESS, IS_CAPTURE_RUNNING] :
			# Drawing the "gun" sight
			img_buffer[img_buffer.shape[0]/2, :] = 0
			img_buffer[:, img_buffer.shape[1]/2] = 0
			# Drawing the circule of 10 micron radius
			img_buffer[ circular_cut ] = 0	
			
			result = RETURN_SUCCESS
		else :
			result = RETURN_FAIL
			print("Error in is_FreezeVideo")
			
		camera_img_buffer.release ()
		return result

	# record time when last photo was taken
	previous_photo_taken_time = time.clock()
	
	# Main loop: Monitoring input queue
	while True :
		command = input_camera_queue.get ()
		
		if command == CMD_EXIT : break # Exiting the process
		elif command == CMD_TAKE_PHOTO :
			# Taking picture if it not already in the buffer
			if time.clock() - previous_photo_taken_time > 0.05:
				# The photo in the buffer is more than 50ms old, take new one
				output_camera_queue.put ( AquireImg() )
				previous_photo_taken_time = time.clock()
			else : # no need to take new picture
				output_camera_queue.put (RETURN_SUCCESS)

		elif command == CMD_TAKE_PHOTO_ACCURATELY :
			# Taking picture 
			output_camera_queue.put ( AquireImg() )
			previous_photo_taken_time = time.clock()
		else : 
			print "Camera error: Unrecognised request"
			output_camera_queue.put (RETURN_FAIL)

	################### Close camera ##################
	if lib.is_ExitCamera(ph_cam) != IS_SUCCESS :
		raise RuntimeError("Error in is_ExitCamera")
	
####################################################################################################

def moveto (new_position, write_serial_port_queue, read_serial_port_queue, \
		unit_size, input_shutter_queue, output_shutter_queue, close_shutter) :
	"""
	Move to a <new_position>, specified as a triple of voltages. 
	<unit_size> must be in nano meters 
	"""
	# t0 = time.time()
	
	# Closing the shutter, if requested
	if close_shutter :
		input_shutter_queue.put (CMD_CLOSE_SHUTTER) 
		if output_shutter_queue.get() != RETURN_SUCCESS : print "\nError in closing the shutter!\n"
		"""	
		# Temporary decrease acceleration of the moving stage
		# Moving stage commands controlling acceleration	
		moving_stage_commands = [ "1AU", "2AU", "1AC", "2AC" ]
		
		# Saving values of the parameters controlling acceleration of the moving stage 
		# Note that at the end these values will be restored
		values_moving_stage_commands = []
		for c in moving_stage_commands :
			write_serial_port_queue.put( (True, c+"?\r") )
			values_moving_stage_commands.append( read_serial_port_queue.get()[:-2] )
	
		# Setting new values to make motion smooth
		command = '0.01;'.join(moving_stage_commands) + "0.01\r"
		write_serial_port_queue.put( (False, command) )
		"""
	
	# Converting voltages into positions 
	new_position = tuple( 1e-6*unit_size*p for p in new_position )

	# Moving 
	write_serial_port_queue.put( (True, "1PA%.6e;2PA%.6e;3PA%.6e;1WS0;2WS0;3WS0;2MD?\r"% new_position ) )

	result =  read_serial_port_queue.get()
	try : 
		if not int(result) : print "Error: Moving stage is still in motion!"
	except ValueError :
		print "Warning: Return from moving stage \"%s\" could not be converted to int" % result 
	
	# Open the shutter, if previously closed
	if close_shutter :
		"""
		# Restoring parameters of the moving stage
		command = ';'.join( map(lambda x : ''.join(x), zip(moving_stage_commands,values_moving_stage_commands) ) ) + "\r"
		write_serial_port_queue.put( (False, command) )
		"""
		# Open the shutter
		input_shutter_queue.put (CMD_OPEN_SHUTTER) 
		if output_shutter_queue.get() != RETURN_SUCCESS : print "\nError in openning the shutter!\n"

		
	# print time.time() - t0

class BaseScanning :
		"""
		This is the base class for any scanning algorithm.
		"""
		def __init__(self, *args) :
			"""
			Constructor
			"""
			# Save the all the arguments
			args_names = ["scaning_event", "pause_scannig_event", 
					"input_pico_harp_queue", "output_pico_harp_queue", "histogram_buffer", 
					"write_serial_port_queue", "read_serial_port_queue", 
					"input_camera_queue", "output_camera_queue", "camera_img_buffer", 
					"input_shutter_queue", "output_shutter_queue", "ScanParameters", "command"] 
			for key, value in zip(args_names, args) : setattr(self, key, value)
			
			# Arguments for function which controls the moving stages 
			self.move_to_args = (self.write_serial_port_queue, self.read_serial_port_queue, 
				self.ScanParameters["unit_size"], self.input_shutter_queue, self.output_shutter_queue)
				
			# Save HDF5 file name
			self.filename = self.ScanParameters["filename"]
						
			# Saving acquisition time
			self.Tacq = self.ScanParameters["Tacq"]
			self.short_scanning_time = 100 # milliseconds (this cannot be faster due to the hardware response time)

			# Ignore points during scanning that have less counts than the following threshold
			self.histogram_counts_threshold = self.ScanParameters["histogram_counts_threshold"]

			# The algorithm is designed to achive the following number of counts
			self.requested_histogram_counts = self.ScanParameters["requested_histogram_counts"]

			# How many times a histogram has to be retaken
			self.revisits_range = range(self.ScanParameters["revisit_number"])

			# Initializing arrays for saving histograms
			self.histogram 				= histogram_shape(self.histogram_buffer)
			self.collected_histogram 	= np.zeros_like(self.histogram)
			
		def pause_scanning (self) :
			"""
			Axillary function. Pausing scanning and saving collected data if self.scaning_event.is_set() 
			"""
			open_shutter = True
			if self.pause_scannig_event.is_set() :
				save_data = True
				while self.pause_scannig_event.is_set() and self.scaning_event.is_set() :
					if save_data :
						# Close shutter
						self.input_shutter_queue.put(CMD_CLOSE_SHUTTER); self.output_shutter_queue.get ()
						# Closing the file
						self.file_scan.close(); del self.file_scan
						# Re-opening file
						self.file_scan = h5py.File (self.filename, 'a', libver='latest')
						# Deleting and opening histogram group
						try : 
							del self.histograms_grp
							self.histograms_grp =  self.file_scan["histograms"]
						except (AttributeError, KeyError) : pass
						# Deleting and opening photo group
						try :
							del self.photos_grp 
							self.photos_grp =  self.file_scan["photos"]		
							# Non need to open the shutter because user works with photos
							open_shutter = False
						except (AttributeError, KeyError) : pass
						save_data = False
						
				# Open the shutter
				if open_shutter :
					self.input_shutter_queue.put(CMD_OPEN_SHUTTER); self.output_shutter_queue.get ()
		
		def __del__ (self) :
			try : self.file_scan.close()
			except AttributeError : pass
			# closing the shutter
			self.input_shutter_queue.put (CMD_CLOSE_SHUTTER); self.output_shutter_queue.get ()
			self.scaning_event.clear()
		
#####################################################################################################
#
#	Study bleaching
#
#####################################################################################################

class BleachingStudy (BaseScanning) : 
	"""
	Record the time histogram continuously from the same point
	"""
	def __init__(self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_BLEACHING_STUDY)

		if isfile (self.filename) :
			# the file is already exists ask user what to do
			app = wx.App()
			choices = ["overwrite", "exit"]
			choiceDialog = wx.SingleChoiceDialog(None, "The file %s already exists. What to do?" % self.filename, \
				"Adaptive scanning...", choices)
			if choiceDialog.ShowModal () == wx.ID_CANCEL or choiceDialog.GetSelection() == 1 :
				# Exiting the scanning process
				self.scaning_event.clear()
				return
				
			# Cleaning
			del app, choiceDialog
		
		print "Measurements of bleaching started..."

		# Opening HDF5 file 
		self.file_scan = h5py.File (self.filename, 'a', libver='latest')
	
		# Saving scanning parameters
		try : del self.file_scan["parameters"]
		except KeyError : pass
		parameters_grp = self.file_scan.create_group("parameters") 
		for key, value in self.ScanParameters.items () : parameters_grp[key] = value
	
		# Close the parameter grope
		del parameters_grp
		
		try : del self.file_scan["histograms"]
		except KeyError : pass
		self.histograms_grp =  self.file_scan.create_group("histograms")
	
		# open the shutter
		self.input_shutter_queue.put (CMD_OPEN_SHUTTER); self.output_shutter_queue.get ()

		# list where time will be saved
		frame_number_to_time = []

		# Number of frames saved to the file
		frame_number = 0
	
		# Saving initial time
		initial_time = time.time()

		# Time when the total count was previously printed
		time_counts_printed = initial_time

		# Recording histograms
		while self.scaning_event.is_set() :
		
			# Measure histogram
			self.input_pico_harp_queue.put( CMD_MEASURE_HISTOGRAM )
			if self.output_pico_harp_queue.get() == RETURN_FAIL :  
				self.scaning_event.clear (); break

			current_time =  time.time()
			frame_number_to_time.append( current_time-initial_time )
			frame_number += 1
		
			# Saving the histogram
			self.histogram_buffer.acquire()
			self.histograms_grp.create_dataset("histogram_%d" % frame_number, data=self.histogram)
			self.histogram_buffer.release()
		
			if current_time - time_counts_printed > 1 : 
				# Printing counts with 1 sec interval
				time_counts_printed = current_time
				print "total count %d" % self.histogram.sum()

			# Check whether user wants to pause
			self.pause_scanning ()

		self.file_scan["parameters"]["frame_number_to_time"] = frame_number_to_time
		print "Measurements of bleaching finished!"

#####################################################################################################
#
#	Different scanning strategies
#
#####################################################################################################

class TestScannig (BaseScanning) :
	"""
	Just move the sample without saving anything 
	"""
	def __init__ (self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_TEST_SCANNING)
		
		# Keep the shutter always closed
		self.move_to_args += (False,)

		# Constant counting number of measurements completed 
		completed = 0
	
		# Constant for measuring time it takes to complete scan	
		initial_time = time.clock()

		# Generate the scanning geometry
		voltages = get_scanning_polygon(self.ScanParameters)
		
		# Repeating measurements until scanning is over or user interrupted
		for position in voltages :

			# Pausing scanning if scaning_event.is_set()
			while self.pause_scannig_event.is_set() and self.scaning_event.is_set() : pass
			
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break
	
			# Going to new position
			moveto(position, *self.move_to_args)
				
			# Print progress info		
			completed += 1
			percentage_completed = 100.*completed / len(voltages)
			seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
			# convert to hours:min:sec
			m, s = divmod(seconds_left, 60)
			h, m = divmod(m, 60)
			print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s )

#####################################################################################################################

class VideScanPath (BaseScanning) :
	"""
	Record video of scanning path.
	"""
	def __init__ (self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_TAKE_SCANNING_VIDEO)
		
		# Opening HDF5 file 
		self.file_scan = h5py.File (self.filename, 'a', libver='latest')

		# Keep the shutter always closed
		self.move_to_args += (False,)
		
		# Saving scanning parameters
		try : del self.file_scan["parameters"]
		except KeyError : pass
		parameters_grp = self.file_scan.create_group("parameters") 
		for key, value in self.ScanParameters.items () : parameters_grp[key] = value

		try : del self.file_scan["photos"]
		except KeyError : pass
		self.photos_grp =  self.file_scan.create_group("photos")

		# Create wrapper for camera image
		camera_img = camera_img_shape(self.camera_img_buffer)

		# If the <histogram> group found in the file, choose <voltages> matching the position of histograms
		try :	
			voltages = map( lambda key : tuple(map(int, key.split('_')[-3:])), self.file_scan["histograms"].iterkeys() )
			voltages = list(set(voltages))

			# the following two lines are from the method <get_scanning_polygon>
			# Order by slow varying axes
			voltages = sorted( voltages, key=operator.itemgetter(0) )
			voltages = sorted( voltages, key=operator.itemgetter(1) )
			voltages = sorted( voltages, key=operator.itemgetter(2) )
			# Ignore random ordering
			try : voltages = sorted( voltages, key=operator.itemgetter(self.ScanParameters["slow_varying_axis"]) )
			except IndexError : pass

			# Saving new <voltages> in the file
			try : del parameters_grp["voltages"]
			except KeyError : pass
			parameters_grp["voltages"] = voltages
			
			print "\nTaking pictures corresponding to the histograms saved in the file.\n"
		except KeyError :
			# Otherwise generate the scanning geometry
			voltages = get_scanning_polygon(self.ScanParameters)

		# Close the parameter grope
		del parameters_grp
	
		# Constant counting number of measurements completed 
		completed = 0
	
		# Constant for measuring time it takes to complete scan	
		initial_time = time.clock()

		# Repeating measurements until scanning is over or user interrupted
		for position in voltages :

			# Pausing scanning and saving collected data if scaning_event.is_set() 
			self.pause_scanning()
			
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break
	
			# Going to new position
			moveto(position, *self.move_to_args)
		
			# Saving photo of current position
			self.input_camera_queue.put (CMD_TAKE_PHOTO_ACCURATELY)
			if self.output_camera_queue.get () == RETURN_FAIL :
				self.scaning_event.clear (); break
			# Saving photo
			self.camera_img_buffer.acquire ()
			self.photos_grp.create_dataset("photo_%d_%d_%d" % position, data=camera_img, compression='szip')	
			self.camera_img_buffer.release ()
		
			# Print progress info		
			completed += 1
			percentage_completed = 100.*completed / len(voltages)
			seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
			# convert to hours:min:sec
			m, s = divmod(seconds_left, 60)
			h, m = divmod(m, 60)
			print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s ) 

#####################################################################################################################

class ScanHistogramsBulk (BaseScanning) :
	"""
	Collecting multiple histograms (the exact number equals to ScanParameters["revisit_number"]+1) at positions specified by <voltages>.
	Multiple histograms are saved before we move to the next point. A point is never returned to. 
	There is a pause between acquiring revision histograms. The duration of the pause equals to the acquisition time.
	"""	
	def __init__ (self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_SCAN)

		# Check whether the file exits
		if not isfile (self.filename) : is_file_to_be_overwritten = True
		else :	# the file is present ask user what to do
			app = wx.App()
			choices = ["resume", "overwrite", "exit"]
			choiceDialog = wx.SingleChoiceDialog(None, "The file %s already exists. What to do?" % self.filename, \
				"Bulk scanning...", choices)
			
			if choiceDialog.ShowModal () == wx.ID_CANCEL or choiceDialog.GetSelection() == 2 :
				# Exiting the scanning process
				self.scaning_event.clear()
				return
			
			is_file_to_be_overwritten = ( choiceDialog.GetSelection() == 1 )
		
			# Cleaning
			del app, choiceDialog		
	
		# Open the HDF5 file
		self.file_scan = h5py.File (self.filename, 'a', libver='latest')

		# Generate the scanning geometry
		voltages = get_scanning_polygon(self.ScanParameters)
	
		if is_file_to_be_overwritten :
			# The content of the HDF5  file will be overwritten 
		
			# Saving scanning parameters
			try : del self.file_scan["parameters"]
			except KeyError : pass
			parameters_grp = self.file_scan.create_group("parameters") 
			for key, value in self.ScanParameters.items () : parameters_grp[key] = value	
			parameters_grp["voltages"] = voltages

			# Close the parameter group
			del parameters_grp

			try : del self.file_scan["histograms"]
			except KeyError : pass
			self.histograms_grp =  self.file_scan.create_group("histograms") 
		else :
			# Resume scanning where left
		
			self.histograms_grp =  self.file_scan["histograms"]

			# Extract the coordinates of the already scanned histograms
			already_scanned = map( lambda key : tuple(map(int, key.split('_')[-3:])), self.histograms_grp.iterkeys() )
		
			if len(already_scanned) > 0 :
				# the following two lines are from the method <get_scanning_polygon>
				# Order by slow varying axes
				already_scanned = sorted( already_scanned, key=operator.itemgetter(0) )
				already_scanned = sorted( already_scanned, key=operator.itemgetter(1) )
				already_scanned = sorted( already_scanned, key=operator.itemgetter(2) )
				# Ignore random ordering
				try : already_scanned = sorted( already_scanned, key=operator.itemgetter(self.ScanParameters["slow_varying_axis"]) )
				except IndexError : pass

				# Ignore all the positions prior to the position of the last scanned histogram
				voltages = voltages[voltages.index(already_scanned[-1]):]
				del voltages[0]

				# clean
				del already_scanned
	
		# Determine the distance between neighboring points
		min_ax1 = min(voltages, key=operator.itemgetter(0))[0]; max_ax1 = max(voltages, key=operator.itemgetter(0))[0] 
		min_ax2 = min(voltages, key=operator.itemgetter(1))[1]; max_ax2 = max(voltages, key=operator.itemgetter(1))[1]
	
		ax1_distance_between_neighbors = ( max_ax1 - min_ax1 )/self.ScanParameters["number_steps_ax1"]
		ax2_distance_between_neighbors = ( max_ax2 - min_ax2 )/self.ScanParameters["number_steps_ax2"]

		# Going to the initial position
		moveto(voltages[0], *self.move_to_args, close_shutter=True)
		previous_position = voltages[0]
  
		# Constant counting number of measurements completed 
		completed = 0
	
		# Count number of points for which <requested_histogram_counts> is reached before finishing all revisions
		bright_points = 0

		# Constant for measuring time it takes to complete scan	
		initial_time = time.clock()

		# Performing measurements on the grid to find out which points should be revisited
		for position in voltages :
			# how many points are checked
			completed += 1
		
			# pause scanning, if requested
			self.pause_scanning()
			
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break
	
			# check whether we move to a neighboring point
			if abs(position[0]-previous_position[0]) < 2*ax1_distance_between_neighbors and \
				abs(position[1]-previous_position[1]) < 2*ax2_distance_between_neighbors :
					# collect histogram for short scanning time 
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
					# Going to the neighbor without closing the shutter
					moveto(position, *self.move_to_args, close_shutter=False)
			else :
				# Going to the new position by closing the shutter
				moveto(position, *self.move_to_args, close_shutter=True)
				# collect histogram for short scanning time 
				self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
	
			# updating the previous position
			previous_position = position

			# waiting till histogram acquisition ends
			if self.output_pico_harp_queue.get() == RETURN_FAIL :  
				self.scaning_event.clear (); break
			
			# Calculate how many counts would be collected if we continue acquisition
			self.histogram_buffer.acquire()
			predicted_total_counts = self.histogram.sum() * float(self.Tacq)/float(self.short_scanning_time)
			self.histogram_buffer.release()
		
			# Decide whether to continue acquisition
			if predicted_total_counts < self.histogram_counts_threshold : continue
			else :  # Total counts are higher than the threshold
				if self.short_scanning_time >= self.Tacq :
					# No need for further scanning, save the histogram
					self.histogram_buffer.acquire()
					self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % position, data=self.histogram, compression='szip')
					self.histogram_buffer.release()
				else :	# Continue scanning
					# First, save already acquired histogram 
					self.histogram_buffer.acquire()
					np.copyto (self.collected_histogram, self.histogram) 
					self.histogram_buffer.release()
					# Continue measuring the histogram till the end of <Tacq>
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq-self.short_scanning_time) )
					if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break
					# Saving the histogram
					self.histogram_buffer.acquire()
					self.collected_histogram += self.histogram
					self.histogram_buffer.release()
					# Updating the count number
					predicted_total_counts = self.collected_histogram.sum()
					if predicted_total_counts > self.histogram_counts_threshold :
						self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % position, data=self.collected_histogram, compression='szip')
					else : continue

				# Since the fluorescence yield from a given point is high, we continue scanning
				for revision in self.revisits_range :
					# Determine whether we need to retake the histogram
					if predicted_total_counts >= self.requested_histogram_counts :
						print "\tBright point. Ignoring further revisions"
						bright_points += 1
						break # Requested number of counts is reached

					# close the shutter
					self.input_shutter_queue.put (CMD_CLOSE_SHUTTER); self.output_shutter_queue.get ()
					# Wait for fluorescence to recover
					time.sleep (self.Tacq*1e-3)
					# Open the shutter
					self.input_shutter_queue.put (CMD_OPEN_SHUTTER); self.output_shutter_queue.get ()  
					# Acquiring the histogram
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq) )
					if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break
					# Adding the histogram to already collected data
					self.histogram_buffer.acquire()
					self.histograms_grp.create_dataset( "revision_%d_histogram_%d_%d_%d" % ( (revision,)+position ), \
						data = self.histogram, compression='szip')
					# updating the total count info
					predicted_total_counts += self.histogram.sum()
					self.histogram_buffer.release()

			# Print progress info		
			percentage_completed = 100.*completed / len(voltages)
			seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
			# convert to hours:min:sec
			m, s = divmod(seconds_left, 60)
			h, m = divmod(m, 60)
			print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s ) 
	
		print "Number of bright points: %d" % bright_points

#####################################################################################################################

class ScanHistogramBoundary (BaseScanning) :
	"""
	Collecting histograms from the boundary. Boundary thickness is specified by "distance_between_adaptive_points"
	"""
	def __init__ (self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_BOUNDARY_SCAN)
	
		if isfile (self.filename) :
			# the file is already exists ask user what to do
			app = wx.App()
			choices = ["overwrite", "exit"]
			choiceDialog = wx.SingleChoiceDialog(None, "The file %s already exists. What to do?" % self.filename, \
				"Boundary scanning...", choices)
			if choiceDialog.ShowModal () == wx.ID_CANCEL or choiceDialog.GetSelection() == 1 :
				# Exiting the scanning process
				self.scaning_event.clear()
				return
			
			# Cleaning
			del app, choiceDialog
	
		# Opening HDF5 file 
		self.file_scan = h5py.File (self.filename, 'a', libver='latest')

		# Saving scanning parameters
		try : del self.file_scan["parameters"]
		except KeyError : pass
		parameters_grp = self.file_scan.create_group("parameters") 
		for key, value in self.ScanParameters.items () : parameters_grp[key] = value	
		
		# Close the parameter grope
		del parameters_grp

		try : del self.file_scan["histograms"]
		except KeyError : pass
		self.histograms_grp =  self.file_scan.create_group("histograms") 
		
		# Computational geometry
		from shapely.geometry import Point, Polygon
		
		# Conversion factor
		mm_per_V = 1e-6*self.ScanParameters["unit_size"] 
		
		# Converting points from mm to units
		points_ax1, points_ax2 = zip(*self.ScanParameters["points"])
		points_ax1 = np.asarray( np.round(np.array(points_ax1) / mm_per_V), dtype=np.int)
		points_ax2 = np.asarray( np.round(np.array(points_ax2) / mm_per_V), dtype=np.int)
		
		# Generating polygon out of specified points
		polygon = Polygon( zip(points_ax1, points_ax2) )
		del points_ax1, points_ax2
		
		# Find coordinate bounds
		min_ax1, min_ax2, max_ax1, max_ax2 = polygon.bounds
		
		# Generate coordinate ranges 
		if self.ScanParameters["slow_varying_axis"] in [0, 2, 3] :
			slow_varying_axis_range = np.linspace( min_ax1, max_ax1, self.ScanParameters["number_steps_ax1"])
			fast_varying_axis_range = np.linspace( min_ax2, max_ax2, self.ScanParameters["number_steps_ax2"])
			self.ORDER2 = lambda slow, fast : (slow, fast)
			self.ORDER3 = lambda slow, fast, ax3 : (slow, fast, ax3)
		else :
			fast_varying_axis_range = np.linspace( min_ax1, max_ax1, self.ScanParameters["number_steps_ax1"])
			slow_varying_axis_range = np.linspace( min_ax2, max_ax2, self.ScanParameters["number_steps_ax2"])
			self.ORDER2 = lambda slow, fast : (fast, slow)
			self.ORDER3 = lambda slow, fast, ax3 : (fast, slow, ax3)
			
		slow_varying_axis_range = np.asarray( np.round(slow_varying_axis_range), dtype=np.int)
		fast_varying_axis_range = np.asarray( np.round(fast_varying_axis_range), dtype=np.int)
		Ax3 = np.linspace(self.ScanParameters["piezo_min_volt_ax3"], 
							self.ScanParameters["piezo_max_volt_ax3"], self.ScanParameters["number_steps_ax3"])
		Ax3 = np.asarray(np.round(Ax3), dtype=np.int)
		
		# This parameter defines the thickness of the boundary
		self.radius_square = self.ScanParameters["distance_between_adaptive_points"]**2 		
		
		# Count number of points for which <requested_histogram_counts> is reached before finishing all revisions
		self.bright_points = 0
		
		# This variable is used for progress info
		completed = 0
		# Constant for measuring time it takes to complete scan	
		initial_time = time.clock()
		
		# Iterating over Z-axis
		for ax3 in Ax3 :
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break
			
			# Iterating over slow varying axis
			for slow in slow_varying_axis_range :
				# Check whether scanning is requested to stop
				if not self.scaning_event.is_set() : break
				
				# Find the fast varying axis range by exclude points outside the polygon for 
				fast_range = map(polygon.contains, map(Point, itertools.product(*self.ORDER2([slow],fast_varying_axis_range)) ) )
				fast_range = fast_varying_axis_range[ np.nonzero(fast_range) ]
				
				# This list is to ensure that the same point is not scanned twice
				self.points_visited = set()
				
				# First, scan in the forward direction
				self.scan_fast_varying_axis(slow, fast_range, ax3)
				
				# Then, scan in the backward direction (if the whole line has not yet been scanned)
				if len(self.points_visited) < len(fast_range) :
					self.scan_fast_varying_axis(slow, fast_range[::-1], ax3)
	
				# Print progress info
				completed += 1
				percentage_completed = 100.*completed / (slow_varying_axis_range.size*Ax3.size)
				seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
				# convert to hours:min:sec
				m, s = divmod(seconds_left, 60)
				h, m = divmod(m, 60)
				print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s ) 
			
			# This command ensures the appropriate thickness along the slow varying axis 
			# for further details see the implementation of <self.scan_fast_varying_axis>
			try : del self.slow_axis_init_value
			except AttributeError : pass
				
		print "Number of bright points: %d" % self.bright_points
		
		
	def scan_fast_varying_axis (self, slow, fast_range, ax3) :
		"""
		Record histograms by moving along the fast varying axis
		"""
		# New scan, thus the shutter should be closed
		close_shutter = True
		# New scan, thus the scanning is not continuous
		continuous_scan = False

		# Iterating over fast varying axis
		for position in itertools.product( *self.ORDER3([slow], fast_range, [ax3]) ) :
			
			# Check whether this point has already been scanned
			if position in self.points_visited :
				close_shutter = True; continue
			else : self.points_visited.add(position)
		
			# pause scanning, if requested
			self.pause_scanning()
			
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break
				
			if close_shutter :
				# Going to the new position by closing the shutter
				moveto(position, *self.move_to_args, close_shutter=True)
				# collect histogram for short scanning time 
				self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
				# Next time no need to close shutter 
				close_shutter = False
			else :
				# collect histogram for short scanning time 
				self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
				# Going to the neighbor without closing the shutter
				moveto(position, *self.move_to_args, close_shutter=False)
			
			# waiting till histogram acquisition ends 
			if self.output_pico_harp_queue.get() == RETURN_FAIL :  
				self.scaning_event.clear (); break
			
			# Calculate how many counts would be collected if we continue acquisition
			self.histogram_buffer.acquire()
			predicted_total_counts = self.histogram.sum() * float(self.Tacq)/float(self.short_scanning_time)
			self.histogram_buffer.release()
		
			# Decide whether to continue acquisition
			if predicted_total_counts < self.histogram_counts_threshold : 
				continuous_scan = False; continue
			else :  # Total counts are higher than the threshold
				if self.short_scanning_time >= self.Tacq :
					# No need for further scanning, save the histogram
					self.histogram_buffer.acquire()
					self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % position, data=self.histogram, compression='szip')
					self.histogram_buffer.release()
				else :	# Continue scanning
					# First, save already acquired histogram 
					self.histogram_buffer.acquire()
					np.copyto (self.collected_histogram, self.histogram)
					self.histogram_buffer.release()
					# Continue measuring the histogram till the end of <Tacq>
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq-self.short_scanning_time) )
					if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break
					# Saving the histogram
					self.histogram_buffer.acquire()
					self.collected_histogram += self.histogram
					self.histogram_buffer.release()
					# Updating the count number
					predicted_total_counts = self.collected_histogram.sum()
					if predicted_total_counts > self.histogram_counts_threshold :
						self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % position, data=self.collected_histogram, compression='szip')
					else : 
						continuous_scan = False; continue

				# Since the fluorescence yield from a given point is high, we continue scanning
				for revision in self.revisits_range :
					# Determine whether we need to retake the histogram
					if predicted_total_counts >= self.requested_histogram_counts :
						print "\tBright point. Ignoring further revisions"
						self.bright_points += 1
						break # Requested number of counts is reached

					# close the shutter
					self.input_shutter_queue.put (CMD_CLOSE_SHUTTER); self.output_shutter_queue.get ()
					# Wait for fluorescence to recover
					time.sleep (self.Tacq*1e-3)
					# Open the shutter
					self.input_shutter_queue.put (CMD_OPEN_SHUTTER); self.output_shutter_queue.get ()  
					# Acquiring the histogram
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq) )
					if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break
					# Adding the histogram to already collected data
					self.histogram_buffer.acquire()
					self.histograms_grp.create_dataset( "revision_%d_histogram_%d_%d_%d" % ( (revision,)+position ), \
							data = self.histogram, compression='szip')
					# updating the total count info
					predicted_total_counts += self.histogram.sum()
					self.histogram_buffer.release()
			
				# Analyse whether to continue scan
				if continuous_scan :
					# Check whether the boundary is sufficiently thick
					# checking thickness of the fast varying axis
					if ( ( boundary_beginning-np.array(position) )**2 ).sum() > self.radius_square :
						try :	
							# checking thickness of the slow varying axis
							if (slow - self.slow_axis_init_value)**2 > self.radius_square : break
						except AttributeError :
							# A fist point is discovered
							self.slow_axis_init_value = slow
				else :
					# This is the beginning of continuous scanning, save the current point
					boundary_beginning = np.array(position)
					continuous_scan = True
		
#####################################################################################################################

class ScanHistogramsAdaptively (BaseScanning) :
	"""
	Collecting histograms intelligently, i.e, the position of collecting the histograms is determined based on 
	the initial guess specified at <voltages>. Note that z-components of <voltages> will be set to <piezo_min_volt_ax3>.  
	"""
 	def __init__ (self, *args) :
		# Initialize
		BaseScanning.__init__ (self, *args)
		# consistency check
		assert (self.command == CMD_ADAPTIVE_SCAN)
			
		if isfile (self.filename) :
			# the file is already exists ask user what to do
			app = wx.App()
			choices = ["overwrite", "exit"]
			choiceDialog = wx.SingleChoiceDialog(None, "The file %s already exists. What to do?" % self.filename, \
				"Adaptive scanning...", choices)
			if choiceDialog.ShowModal () == wx.ID_CANCEL or choiceDialog.GetSelection() == 1 :
				# Exiting the scanning process
				self.scaning_event.clear()
				return
			
			# Cleaning
			del app, choiceDialog

		# Opening HDF5 file 
		self.file_scan = h5py.File (self.filename, 'a', libver='latest')

		# Saving scanning parameters
		try : del self.file_scan["parameters"]
		except KeyError : pass
		parameters_grp = self.file_scan.create_group("parameters") 
		for key, value in self.ScanParameters.items () : parameters_grp[key] = value	
		
		# Close the parameter grope
		del parameters_grp

		try : del self.file_scan["histograms"]
		except KeyError : pass
		self.histograms_grp =  self.file_scan.create_group("histograms") 
		
		# Generate the scanning geometry
		voltages = get_scanning_polygon(self.ScanParameters)
		
		# Getting rid off z-components in the initial region <voltages>
		piezo_min_volt_ax3 = self.ScanParameters["piezo_min_volt_ax3"]
		voltages = filter( lambda p : p[-1] == piezo_min_volt_ax3, voltages)
		Ax1, Ax2, _ = zip(*voltages)
		voltages = zip(Ax1, Ax2)
		
		# Copy initial coordinates for the adaptive construction
		Ax1_span = np.array(Ax1, copy=True); Ax2_span = np.array(Ax2, copy=True) 
		del Ax1, Ax2

		# Determine the distance between neighboring points
		ax1_distance_between_neighbors = ( Ax1_span.max() - Ax1_span.min() )/self.ScanParameters["number_steps_ax1"]
		ax2_distance_between_neighbors = ( Ax2_span.max() - Ax2_span.min() )/self.ScanParameters["number_steps_ax2"]

		# Points to be used in the Basian analysis
		points_to_revisit = []

		# Count number of points for which <requested_histogram_counts> is reached before finishing all revisions
		bright_points = 0

		# Initializing parameters for adaptive searching
		radius_square = int(self.ScanParameters["distance_between_adaptive_points"])**2 		# distance between adaptevly generated points
		weight_factor =	float(self.ScanParameters["adaptive_scan_point_weight"]) / 100.

		# Repeating measurements until scanning is over or user interrupted
		Ax3 = np.linspace(piezo_min_volt_ax3,self.ScanParameters["piezo_max_volt_ax3"],self.ScanParameters["number_steps_ax3"])
		Ax3 = np.asarray(np.round(Ax3), dtype=np.int)
		for volt_ax3 in Ax3 : 
			
			# Check whether scanning is requested to stop
			if not self.scaning_event.is_set() : break

			# There are no points to scan
			if len(voltages) == 0 : break

			# Going to the initial point (i.e., origin) to avoid the backslash
			moveto( (Ax1_span[0], Ax2_span[0], volt_ax3), *self.move_to_args, close_shutter=True)

			# Save the current position
			previous_position_ax1 = Ax1_span[0]; previous_position_ax2 = Ax2_span[0]

			print "Number of points to visit ", len(voltages)

			print "\n New Ax3 position = %d \n" % volt_ax3

			# Constant counting number of measurements completed 
			completed = 0
		
			# Constant for measuring time it takes to complete the Ax1-Ax2 scan	
			initial_time = time.clock()
		
			################################################################################################
			# Iterating over Ax1 and Ax2
			################################################################################################
			for volt_ax1, volt_ax2 in voltages :

				# how many points are checked
				completed += 1

				# pause scanning, if requested
				self.pause_scanning()

				# Check whether scanning is requested to stop
				if not self.scaning_event.is_set() : break

				# check whether we move to a neighboring point
				if abs(volt_ax1-previous_position_ax1) < 2*ax1_distance_between_neighbors and \
					abs(volt_ax2-previous_position_ax2) < 2*ax2_distance_between_neighbors :
						# collect histogram for short scanning time 
						self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
						# Going to the neighbor without closing the shutter
						moveto( (volt_ax1, volt_ax2, volt_ax3), *self.move_to_args, close_shutter=False)			
				else :
					# Going to the new position by closing the shutter
					moveto( (volt_ax1, volt_ax2, volt_ax3), *self.move_to_args, close_shutter=True)
					# collect histogram for short scanning time 
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.short_scanning_time) )
		
				# updating the previous position
				previous_position_ax1 = volt_ax1; previous_position_ax2 = volt_ax2

				# waiting till histogram acquisition ends 
				if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break

				# Calculate how many counts would be get if we continue acquisition
				self.histogram_buffer.acquire()
				predicted_total_counts = self.histogram.sum() * float(self.Tacq)/float(self.short_scanning_time)
				self.histogram_buffer.release()

				# Decide whether to continue acquisition
				if predicted_total_counts < self.histogram_counts_threshold : continue
				else :	# Total counts are higher than the threshold
					if self.short_scanning_time >= self.Tacq :
						# No need for further scanning, save the histogram
						self.histogram_buffer.acquire()
						self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % (volt_ax1, volt_ax2, volt_ax3), data=self.histogram, compression='szip')
						self.histogram_buffer.release()
					else :	# Continue scanning
						# First, save already acquired histogram 
						self.histogram_buffer.acquire()
						np.copyto (self.collected_histogram, self.histogram) 
						self.histogram_buffer.release()
						# Continue measuring the histogram till the end of <Tacq>
						self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq-self.short_scanning_time) )
						if self.output_pico_harp_queue.get() == RETURN_FAIL :  
							self.scaning_event.clear (); break
						# Saving the histogram
						self.histogram_buffer.acquire()
						self.collected_histogram += self.histogram
						self.histogram_buffer.release()
						# Updating the count number
						predicted_total_counts = self.collected_histogram.sum()
						if predicted_total_counts > self.histogram_counts_threshold :
							# Save histogram only if it has enough counts
							self.histograms_grp.create_dataset( "histogram_%d_%d_%d" % (volt_ax1, volt_ax2, volt_ax3), \
								data = self.collected_histogram, compression='szip')
						else : continue

				################################################################################################
				#  Revisiting Ax1 and Ax2 points were signal is non negligible
				################################################################################################
				for revision in self.revisits_range :
					# Determine whether we need to retake the histogram
					if predicted_total_counts >= self.requested_histogram_counts :
						print "\tBright point. Ignoring further revisions"
						bright_points += 1
						break # Desaired number of counts is reached

					# close the shutter
					self.input_shutter_queue.put (CMD_CLOSE_SHUTTER); self.output_shutter_queue.get ()
					# Wait for fluorescence to recover
					time.sleep (self.Tacq*1e-3)
					# Open the shutter
					self.input_shutter_queue.put (CMD_OPEN_SHUTTER); self.output_shutter_queue.get ()  
					# Acquiring the histogram
					self.input_pico_harp_queue.put( (CMD_MEASURE_HISTOGRAM, self.Tacq) )
					if self.output_pico_harp_queue.get() == RETURN_FAIL :  
						self.scaning_event.clear (); break
					# Adding the histogram to already collected data
					self.histogram_buffer.acquire()
					self.histograms_grp.create_dataset( "revision_%d_histogram_%d_%d_%d" % (revision, volt_ax1, volt_ax2, volt_ax3), \
							data = self.histogram, compression='szip')
					# updating the total count info
					predicted_total_counts += self.histogram.sum()
					self.histogram_buffer.release()

				# Saving this point for the Baisian analysis
				points_to_revisit.append( (volt_ax1, volt_ax2) )

				# Print progress info		
				percentage_completed = 100.*completed / len(voltages)
				seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
				# convert to hours:min:sec
				m, s = divmod(seconds_left, 60)
				h, m = divmod(m, 60)
				print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s ) 

			print "Number of bright points: %d" % bright_points

			################################################################################################
			# Naive Basian Inference to recalculate new Ax1 and Ax2 for the next value of <volt_ax3>
			################################################################################################

			# Calculate weights for each observed point 
			weights = np.zeros (Ax1_span.shape, dtype=np.float)
			for volt_ax1, volt_ax2 in points_to_revisit :
				weights += np.exp( -( (Ax1_span-volt_ax1)**2 + (Ax2_span-volt_ax2)**2 )//radius_square )
			weights *= weight_factor
			
			# Making a decision which points to scan at next iteration
			indx = np.nonzero(weights > 1)
			voltages = zip(Ax1_span[indx], Ax2_span[indx])

#####################################################################################################################

def scan_sample (*args):
	"""
	This function runs as a separate process scanning the sample.
	Communication are done through the boolean flag <scaning_event> and <pause_scanning_event>.
	Note: Pico Harp is assumed to be initialized! 
	<command> specifies what to do.
	"""
	command = args[-1]
	if  command == CMD_SCAN 				: ScanHistogramsBulk (*args)
	elif command == CMD_ADAPTIVE_SCAN 		: ScanHistogramsAdaptively (*args)
	elif command == CMD_BOUNDARY_SCAN 		: ScanHistogramBoundary (*args)
	elif command == CMD_TAKE_SCANNING_VIDEO : VideScanPath (*args)
	elif command == CMD_BLEACHING_STUDY 	: BleachingStudy (*args)
	elif command == CMD_TEST_SCANNING 		: TestScannig (*args)		
	else : raise ValueError ("Unrecognized value of <command>")

	print "Stop scanning"

####################################################################################################
#
#	End scanning 
#
####################################################################################################

def control_shutter (input_shutter_queue, output_shutter_queue) : 
	"""
	This function, run as a separate process, controls the Thorlabls shutter. All communications are done through queues.
	"""
	# Creating ActiveX control used to control the shutter
	import wx.lib.activex, threading

	class Shutter(wx.lib.activex.ActiveXCtrl) :

		def __init__ (self) :
			self.frame = wx.Frame (None, title="Shutter")
			panel = wx.Panel (self.frame)
			wx.lib.activex.ActiveXCtrl.__init__ (self, panel, 'MGMOTOR.MGMotorCtrl.1', size=self.frame.GetSizeTuple(), name='Shutter')
			# Initializing the device by specifying its serial number 
			self.ctrl.HWSerialNum = 85844994
			self.ctrl.StartCtrl()
			#self.open(return_success=False)
			# Building simple gui
			sizer = wx.BoxSizer ()
			sizer.Add (self, flag=wx.EXPAND)
			panel.SetSizer (sizer)
			self.frame.Show()
			#self.frame.Hide()
			# starting thread for checking commands sent from other processes
			threading.Thread(target=self.run_to_check_input_queue).start()

		def close (self, return_success=True) : 
			self.ctrl.SC_Disable(0)  
			if return_success : output_shutter_queue.put (RETURN_SUCCESS)

		def open (self, return_success=True) : 
			self.ctrl.SC_Enable(0)
			if return_success : output_shutter_queue.put (RETURN_SUCCESS)

		def run_to_check_input_queue (self) :
			"""
			This function run on a separate thread checks the queue.
			"""
			while True :
				command = input_shutter_queue.get() 
		
				if command == CMD_CLOSE_SHUTTER : 	wx.CallAfter(self.close)
				elif command == CMD_OPEN_SHUTTER : 	wx.CallAfter(self.open)
				elif command == CMD_EXIT :		
									wx.CallAfter(self.frame.Close); break
				else : 
					output_shutter_queue.put (RETURN_FAIL)
					print "\nUnrecognized commands sent to the shutter!\n"
	# Starting the application
	wx_app = wx.App(False)
	shutter = Shutter()
	wx_app.MainLoop()

####################################################################################################

def control_moving_stage (write_serial_port_queue, read_serial_port_queue) :
	"""
	This is function runs as a separate process and controls the NewPort 
	moving stage through the serial port. All communications are done through the read and write queues.
	"""
	import serial
	
	# Opening serial port to talk to NewPort moving stage
	serial_port = serial.Serial (port='COM5', baudrate=19200,  bytesize=8, parity=serial.PARITY_NONE, stopbits=1, timeout=0.5)
	
	while True :
		full_command = write_serial_port_queue.get ()
		if full_command == CMD_EXIT : break
	
		# Sending a command through the COM port
		is_read, command = full_command
		serial_port.write (command)
		if is_read : read_serial_port_queue.put (serial_port.readline())
	
	serial_port.close ()

####################################################################################################

def joystick_control_moving_stage (use_joystick_event, write_serial_port_queue, read_serial_port_queue) :
	"""
	Control the moving stages with joystick
	"""
	import pygame
	
	# Initializing joystick
	pygame.init()
	pygame.joystick.init()
	
	if pygame.joystick.get_count() == 0 or not pygame.joystick.get_init() :
		print "No joystick is detected!"
		return

	if pygame.joystick.get_count() <> 1 :
		print "Multiple joysticks are detected: the first one will be utilized"

	Joystick = pygame.joystick.Joystick(0)
	Joystick.init ()
	print 'Joystick "%s" was initialized' % Joystick.get_name() 

	# Linking joystick to the moving stage
	while use_joystick_event.is_set() :

		# reading off the position of the joystick 
		# Note that the labels are transposed to align the joystick with the moving stage
		possition0 = Joystick.get_axis(0)
		possition1 = -Joystick.get_axis(1)

		if abs(possition0) > 1e-2 or abs(possition1) > 1e-2 :  
			# User moved joystick
			
			# Check whether the moving stage does not move
			write_serial_port_queue.put( (True, "1MD?\r") )
			ax1 = int(read_serial_port_queue.get())
			write_serial_port_queue.put( (True, "2MD?\r") )
			ax2 = int(read_serial_port_queue.get())
			
			if ax1*ax2 :
				# The rudder control is used to scale the step size
				scale = 0.05*(1. - Joystick.get_axis(2))
				possition0 *= scale
				possition1 *= -scale

				# The moving stage is ready to be moved
				command = "1PR%e;2PR%e\r" % (possition0, possition1)
				write_serial_port_queue.put ( (False, command) ) 
		
		# Clearing the event stack
		pygame.event.clear ()	

	# Releasing resources
	Joystick.quit ()
	pygame.joystick.quit()
	pygame.quit ()

####################################################################################################

class CMicroscopeConrol (wx.Frame) :
	"""
	Application for controlling the microscope 
	"""

	def __init__ (self, parent, title=None) :
		"""
		Constructor
		"""
		# Create process controlling camera
		self.input_camera_queue 	= multiprocessing.Queue () 
		self.output_camera_queue	= multiprocessing.Queue ()
		# Buffer for storing image
		self.camera_img_buffer = multiprocessing.Array (ctypes.c_uint8, CAMERA_IMG_WIDTH*CAMERA_IMG_HIGHT) 

		self.control_camera = multiprocessing.Process (target=control_camera, \
						args=(self.input_camera_queue, self.output_camera_queue, self.camera_img_buffer) )
		self.control_camera.start ()

		# Create process controlling the NewPort Moving stage
		self.write_serial_port_queue 	= multiprocessing.Queue ()  
		self.read_serial_port_queue	= multiprocessing.Queue ()
	
		self.control_moving_stage = multiprocessing.Process (target=control_moving_stage, \
								args=(self.write_serial_port_queue, self.read_serial_port_queue))
		self.control_moving_stage.start ()

		# Create process controlling Pico Harp 3000
		self.input_pico_harp_queue	= multiprocessing.Queue ()
		self.output_pico_harp_queue	= multiprocessing.Queue ()
		# Buffer for storing histogram
		self.histogram_buffer = multiprocessing.Array (ctypes.c_uint, HISTCHAN)

		self.control_pico_harp = multiprocessing.Process (target=control_pico_harp, \
								args=(self.input_pico_harp_queue, self.output_pico_harp_queue, self.histogram_buffer))
		self.control_pico_harp.start ()

		# Creating process controlling shutter
		self.input_shutter_queue 	= multiprocessing.Queue ()
		self.output_shutter_queue	= multiprocessing.Queue ()
		self.control_shutter = multiprocessing.Process (target=control_shutter, args=(self.input_shutter_queue, self.output_shutter_queue))
		self.control_shutter.start()

		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title=title, size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Show ()
		wx.EVT_CLOSE (self, self.on_close)

	def on_close (self, event):
		"""
		Windows is about to be closed. Stop all timers.
		"""
		self.StopAllJobs ()
		self.Destroy ()

	def __del__ (self) :
		# Closing joystick control
		try :
			self.__use_joystick_event__.clear()
			self.joystick_process.join ()
		except AttributeError : pass

		# Closing shutter
		self.input_shutter_queue.put (CMD_EXIT)
		self.control_shutter.join ()
		self.input_shutter_queue.close ()
		self.output_shutter_queue.close ()

		# Closing camera
		self.input_camera_queue.put (CMD_EXIT)
		self.control_camera.join ()
		self.output_camera_queue.close ()
		self.input_camera_queue.close ()

		# Closing the NewPort moving stage
		self.write_serial_port_queue.put (CMD_EXIT)
		self.control_moving_stage.join()
		self.write_serial_port_queue.close ()
		self.read_serial_port_queue.close ()
	
		# Closing Pico Harp 3000
		self.input_pico_harp_queue.put (CMD_EXIT)
		self.control_pico_harp.join()
		self.input_pico_harp_queue.close ()
		self.output_pico_harp_queue.close ()
	
	def ConstructGUI (self) :
		"""
		Build GUI
		"""
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		###################################### Panel Pico Harp properties ##########################
		sb = wx.StaticBox(panel, label="Pico Harp settings")
		boxsizer = wx.StaticBoxSizer(sb, wx.VERTICAL)
		
		boxsizer.Add (wx.StaticText(panel, label="Offset (ns)"), flag=wx.LEFT|wx.TOP, border=5)
		self.PicoHarp_Offset = wx.SpinCtrl (panel, value="0", min=0, max=1e6)
		boxsizer.Add (self.PicoHarp_Offset, flag=wx.EXPAND, border=5)
		def update_offset_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_OFFSET_PICO_HARP, self.PicoHarp_Offset.GetValue()) )
		self.PicoHarp_Offset.Bind (wx.EVT_SPINCTRL, update_offset_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="CFDZeroX0 (mV)"), flag=wx.LEFT, border=5)
		self.PicoHarp_CFDZeroX0 = wx.SpinCtrl (panel, value="10", min=0, max=1e3)
		boxsizer.Add (self.PicoHarp_CFDZeroX0, flag=wx.EXPAND, border=5)
		def update_CFDZeroX0_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_CFDZERO_PICO_HARP, 0, self.PicoHarp_CFDZeroX0.GetValue()) )
		self.PicoHarp_CFDZeroX0.Bind (wx.EVT_SPINCTRL, update_CFDZeroX0_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="CFDLevel0 (mV)"), flag=wx.LEFT, border=5)
		self.PicoHarp_CFDLevel0 = wx.SpinCtrl (panel, value="80", min=0, max=1e3)
		boxsizer.Add (self.PicoHarp_CFDLevel0, flag=wx.EXPAND, border=5)
		def update_CFDLevel0_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_CFDLEVEL_PICO_HARP, 0, self.PicoHarp_CFDLevel0.GetValue()) )
		self.PicoHarp_CFDLevel0.Bind (wx.EVT_SPINCTRL, update_CFDLevel0_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="CFDZeroX1 (mV)"), flag=wx.LEFT, border=5)
		self.PicoHarp_CFDZeroX1 = wx.SpinCtrl (panel, value="10", min=0, max=1e3)
		boxsizer.Add (self.PicoHarp_CFDZeroX1, flag=wx.EXPAND, border=5)
		def update_CFDZeroX1_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_CFDZERO_PICO_HARP, 1, self.PicoHarp_CFDZeroX1.GetValue()) )
		self.PicoHarp_CFDZeroX1.Bind (wx.EVT_SPINCTRL, update_CFDZeroX1_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="CFDLevel1 (mV)"), flag=wx.LEFT, border=5)
		self.PicoHarp_CFDLevel1 = wx.SpinCtrl (panel, value="70", min=0, max=1e3)
		boxsizer.Add (self.PicoHarp_CFDLevel1, flag=wx.EXPAND, border=5)
		def update_CFDLevel1_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_CFDLEVEL_PICO_HARP, 1, self.PicoHarp_CFDLevel1.GetValue()) )
		self.PicoHarp_CFDLevel1.Bind (wx.EVT_SPINCTRL, update_CFDLevel1_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="Divider"), flag=wx.LEFT, border=5)
		self.PicoHarp_SyncDiv = wx.SpinCtrl (panel, value="8", min=1, max=1e3)
		boxsizer.Add (self.PicoHarp_SyncDiv, flag=wx.EXPAND, border=5)
		def update_SyncDiv_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_DIVIDER_PICO_HARP, self.PicoHarp_SyncDiv.GetValue()) )
		self.PicoHarp_SyncDiv.Bind (wx.EVT_SPINCTRL, update_SyncDiv_pico_harp)
					
		boxsizer.Add (wx.StaticText(panel, label="Range"), flag=wx.LEFT, border=5)
		self.PicoHarp_Range = wx.SpinCtrl (panel, value="2", min=0, max=1e6)
		boxsizer.Add (self.PicoHarp_Range, flag=wx.EXPAND, border=5)
		def update_Range_pico_harp (event) :
			self.input_pico_harp_queue.put( (CMD_SET_RANGE_PICO_HARP, self.PicoHarp_Range.GetValue()) )
		self.PicoHarp_Range.Bind (wx.EVT_SPINCTRL, update_Range_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="Acquisition time (ms)"), flag=wx.LEFT, border=5)
		self.PicoHarp_Tacq = wx.SpinCtrl (panel, value="1000", min=100, max=1e6)
		boxsizer.Add (self.PicoHarp_Tacq, flag=wx.EXPAND, border=5)
		def update_Tacq_pico_harp (event) :
			T = self.PicoHarp_Tacq.GetValue()
			self.input_pico_harp_queue.put( (CMD_SET_TACQ_PICO_HARP, T) )
			# Modifying the timer used to measure single histogram
			try :
				self.histogram_timer.Stop()
				self.histogram_timer.Start(T)
			except AttributeError : pass
		self.PicoHarp_Tacq.Bind (wx.EVT_SPINCTRL, update_Tacq_pico_harp)

		boxsizer.Add (wx.StaticText(panel, label="Time resolution (ps)"), flag=wx.LEFT, border=5)
		self.PicoHarp_Resolution = wx.SpinCtrl (panel, value="0", min=0, max=1e6)
		boxsizer.Add (self.PicoHarp_Resolution, flag=wx.EXPAND, border=5)
		self.PicoHarp_Resolution.Disable ()

		sizer.Add(boxsizer, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)
		
		###################### Panel Scan settings ################################
		sb = wx.StaticBox(panel, label="Scanning settings")
		boxsizer = wx.StaticBoxSizer(sb, wx.VERTICAL)
		
		# File name for saving results of scanning 
		boxsizer.Add (wx.StaticText(panel, label="File name to save scan"), flag=wx.LEFT|wx.TOP, border=5)
		self.scan_filename = wx.TextCtrl (panel, value="result.hdf5")
		boxsizer.Add (self.scan_filename,  flag=wx.EXPAND, border=5)
		
		# Info about scanning
		boxsizer.Add (wx.StaticText(panel, label="Info"), flag=wx.LEFT|wx.TOP, border=5)
		self.scan_info = wx.TextCtrl (panel, style=wx.TE_MULTILINE)
		boxsizer.Add (self.scan_info,  flag=wx.EXPAND, border=5)
		
		# How many times a given point will be revisited during scanning
		boxsizer.Add (wx.StaticText(panel, label="\nNumber of re-visits"), flag=wx.LEFT, border=5)
		self.revisit_number = wx.SpinCtrl (panel, value="3", min=0, max=200)
		boxsizer.Add (self.revisit_number,  flag=wx.EXPAND, border=5)
		
		# Weight of individual point during adaptive scan  
		boxsizer.Add (wx.StaticText(panel, label="Adaptive weight (%)"), flag=wx.LEFT, border=5)
		self.adaptive_scan_point_weight = wx.SpinCtrl (panel, value="80", min=0, max=200)
		boxsizer.Add (self.adaptive_scan_point_weight,  flag=wx.EXPAND, border=5)
	
		# Spacing between added points during adaptive scan
		boxsizer.Add (wx.StaticText(panel, label="Adaptive radius (units)"), flag=wx.LEFT, border=5)
		self.distance_between_adaptive_points = wx.SpinCtrl (panel, value="100", min=0, max=1e6)
		boxsizer.Add (self.distance_between_adaptive_points,  flag=wx.EXPAND, border=5)

		# Ignore histograms that have counts below following threshold
		boxsizer.Add (wx.StaticText(panel, label="Ignore points with counts under"), flag=wx.LEFT, border=5)
		self.histogram_counts_threshold = wx.SpinCtrl (panel, value="15000", min=0, max=1e7)
		boxsizer.Add (self.histogram_counts_threshold,  flag=wx.EXPAND, border=5)

		# Parameters used for adaptive scanning: Desired number of counts in each histogram
		boxsizer.Add (wx.StaticText(panel, label="Requested histogram counts"), flag=wx.LEFT, border=5)
		self.requested_histogram_counts = wx.SpinCtrl (panel, value="15000", min=0, max=1e7)
		boxsizer.Add (self.requested_histogram_counts,  flag=wx.EXPAND, border=5)

		# Number of steps
		boxsizer.Add (wx.StaticText(panel, label="\nNumber of steps in Ax1"), flag=wx.LEFT, border=5)
		self.number_steps_ax1 = wx.SpinCtrl (panel, value="100", min=1, max=1e6)
		boxsizer.Add (self.number_steps_ax1, flag=wx.EXPAND, border=5)

		boxsizer.Add (wx.StaticText(panel, label="Number of steps in Ax2"), flag=wx.LEFT, border=5)
		self.number_steps_ax2 = wx.SpinCtrl (panel, value="100", min=1, max=1e6)
		self.number_steps_ax2.SetRange (1, 1e6)
		boxsizer.Add (self.number_steps_ax2, flag=wx.EXPAND, border=5)

		boxsizer.Add (wx.StaticText(panel, label="Number of steps in Ax3"), flag=wx.LEFT, border=5)
		self.number_steps_ax3 = wx.SpinCtrl (panel, value="1", min=1, max=1e6)
		self.number_steps_ax3.SetRange (1, 1e6)
		boxsizer.Add (self.number_steps_ax3, flag=wx.EXPAND, border=5)
		
		# Voltage characteristics of piezo
		boxsizer.Add (wx.StaticText(panel, label="\nMin Ax3 (units)"), flag=wx.LEFT, border=5)
		self.piezo_min_volt_ax3 = wx.SpinCtrl (panel, value="0", min=-1e6, max=1e6)
		boxsizer.Add (self.piezo_min_volt_ax3, flag=wx.EXPAND, border=5)
		
		boxsizer.Add (wx.StaticText(panel, label="Max Ax3 (units)"), flag=wx.LEFT, border=5)
		self.piezo_max_volt_ax3 = wx.SpinCtrl (panel, value="100", min=-1e6, max=1e6)
		boxsizer.Add (self.piezo_max_volt_ax3, flag=wx.EXPAND, border=5)

		# Unit of distance used for moving stages
		boxsizer.Add (wx.StaticText(panel, label="Unit size in nm"), flag=wx.LEFT, border=5)
		self.unit_size = wx.SpinCtrl (panel, value="100", min=0, max=1e6)
		boxsizer.Add (self.unit_size, flag=wx.EXPAND, border=5)

		# Fast varying axes
		sb_ = wx.StaticBox(panel, label="\nSlow varying axis")
		boxsizer_ = wx.StaticBoxSizer(sb_, wx.HORIZONTAL)
		
		self.ax1_slow = wx.RadioButton(panel, label="Ax1", style=wx.RB_GROUP)
		boxsizer_.Add (self.ax1_slow, flag=wx.EXPAND, border=5)
		self.ax2_slow = wx.RadioButton(panel, label="Ax2")
		boxsizer_.Add (self.ax2_slow, flag=wx.EXPAND, border=5)
		self.ax3_slow = wx.RadioButton(panel, label="Ax3")
		boxsizer_.Add (self.ax3_slow, flag=wx.EXPAND, border=5)
		self.random_ax_slow = wx.RadioButton(panel, label="rand")
		boxsizer_.Add (self.random_ax_slow, flag=wx.EXPAND, border=5)

		boxsizer.Add (boxsizer_, flag=wx.EXPAND, border=10)
		
		sizer.Add(boxsizer, pos=(0, 1), span=(2, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		###################### Panel scanning geometry ###############################
		sb = wx.StaticBox(panel, label="Scanning geometry")
		boxsizer = wx.StaticBoxSizer(sb, wx.VERTICAL)

		# Buttons for specifying the boundary of scanning geometry

		self.add_point_button = wx.Button (panel, label="Add point")
		self.Bind (wx.EVT_BUTTON, self.add_point, self.add_point_button)
		boxsizer.Add (self.add_point_button, flag=wx.EXPAND, border=5)
		# List where all points are stored
		self.points = []

		self.delete_points_button = wx.Button (panel, label="Delete all points")
		
		def delete_points (event) :
			self.points = [] # Remove all points
			# Delete all associated photos
			for key in dir(self) :
				if "photo_" in key : delattr(self, key)
			# Changing the color
			self.add_point_button.SetBackgroundColour('')
			# Change the label
			self.go_to_points_ABCDE_button.SetLabel("Go to point_0")

		self.Bind (wx.EVT_BUTTON, delete_points, self.delete_points_button)
		boxsizer.Add (self.delete_points_button, flag=wx.EXPAND, border=5)

		# Button for drawing scanning geometry
		self.draw_scanning_geometry_button = wx.Button (panel, label="Draw scanning geometry")
		self.Bind (wx.EVT_BUTTON, self.draw_scanning_geometry, self.draw_scanning_geometry_button)
		boxsizer.Add(self.draw_scanning_geometry_button, flag=wx.EXPAND|wx.TOP, border=5)
		
		# Button for testing scanning geometry
		self.test_scanning_geometry_button = wx.Button (panel)
		self.test_scanning_geometry_button.__start_label__ = "Test scanning geometry"
		self.test_scanning_geometry_button.__pause_label__ = "PAUSE testing geometry"
		self.test_scanning_geometry_button.__resume_label__ = "RESUME testing geometry"
		self.test_scanning_geometry_button.__stop_label__ = "STOP testing geometry"
		self.test_scanning_geometry_button.SetLabel (self.test_scanning_geometry_button.__start_label__)
		self.test_scanning_geometry_button.command = CMD_TEST_SCANNING
		self.test_scanning_geometry_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.test_scanning_geometry_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.test_scanning_geometry_button, flag=wx.EXPAND|wx.TOP, border=5)

		sizer.Add(boxsizer, pos=(2, 1), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		###################### Panel Operations #######################################
		sb = wx.StaticBox(panel, label="Operations")
		boxsizer = wx.StaticBoxSizer(sb, wx.VERTICAL)

		# Perform bulk scanning (button)
		self.scanning_button = wx.Button (panel)
		self.scanning_button.__start_label__ = "Start bulk scanning"
		self.scanning_button.__pause_label__ = "PAUSE bulk scanning"
		self.scanning_button.__resume_label__ = "RESUME bulk scanning"
		self.scanning_button.__stop_label__ = "STOP bulk scanning"
		self.scanning_button.command = CMD_SCAN
		self.scanning_button.SetLabel (self.scanning_button.__start_label__)
		self.scanning_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.scanning_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.scanning_button, flag=wx.EXPAND, border=5)
	
		# Perform adaptive scanning (button)
		self.adaptive_scanning_button = wx.Button (panel)
		self.adaptive_scanning_button.__start_label__ = "Start adaptive scanning"
		self.adaptive_scanning_button.__pause_label__ = "PAUSE adaptive scanning"
		self.adaptive_scanning_button.__resume_label__ = "RESUME adaptive scanning"
		self.adaptive_scanning_button.__stop_label__ = "STOP adaptive scanning"
		self.adaptive_scanning_button.command = CMD_ADAPTIVE_SCAN
		self.adaptive_scanning_button.SetLabel (self.adaptive_scanning_button.__start_label__)
		self.adaptive_scanning_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.adaptive_scanning_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.adaptive_scanning_button, flag=wx.EXPAND, border=5)
		
		# Perform scanning of boundary (button)
		self.boundary_scanning_button = wx.Button (panel)
		self.boundary_scanning_button.__start_label__ = "Start boundary scanning"
		self.boundary_scanning_button.__pause_label__ = "PAUSE boundary scanning"
		self.boundary_scanning_button.__resume_label__ = "RESUME boundary scanning"
		self.boundary_scanning_button.__stop_label__ = "STOP boundary scanning"
		self.boundary_scanning_button.command = CMD_BOUNDARY_SCAN
		self.boundary_scanning_button.SetLabel (self.boundary_scanning_button.__start_label__)
		self.boundary_scanning_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.boundary_scanning_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.boundary_scanning_button, flag=wx.EXPAND, border=5)

		# Button for recording video of scanning path
		self.scanning_video_button = wx.Button (panel) 
		self.scanning_video_button.__start_label__ = "Video of scanning path"
		self.scanning_video_button.__pause_label__ = "PAUSE Video"
		self.scanning_video_button.__resume_label__ = "RESUME Video"
		self.scanning_video_button.__stop_label__ = "STOP Video"
		self.scanning_video_button.SetLabel (self.scanning_video_button.__start_label__)
		self.scanning_video_button.command = CMD_TAKE_SCANNING_VIDEO
		self.scanning_video_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.scanning_video_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.scanning_video_button, flag=wx.EXPAND, border=5)

		# Button for recording bleaching
		self.bleaching_study_button = wx.Button (panel) 
		self.bleaching_study_button.__start_label__ = "Bleaching study"
		self.bleaching_study_button.__pause_label__ = "PAUSE bleaching study"
		self.bleaching_study_button.__resume_label__ = "RESUME bleaching study"
		self.bleaching_study_button.__stop_label__ = "STOP bleaching study"
		self.bleaching_study_button.SetLabel (self.bleaching_study_button.__start_label__)
		self.bleaching_study_button.command = CMD_BLEACHING_STUDY
		self.bleaching_study_button.Bind (wx.EVT_LEFT_DOWN, self.do_scanning)
		self.bleaching_study_button.Bind (wx.EVT_LEFT_DCLICK, self.do_scanning)
		boxsizer.Add(self.bleaching_study_button, flag=wx.EXPAND, border=5)

		# Button for getting a single histogram
		self.__start_measuring_single_histogram_label__ 	= "Measure single histogram"
		self.__stop_measuring_single_histogram_label__		= "STOP measuring histogram"

		self.measure_single_histogram_button = wx.Button (panel, label=self.__start_measuring_single_histogram_label__)
		self.Bind (wx.EVT_BUTTON, self.measure_single_histogram, id=self.measure_single_histogram_button.GetId())
		boxsizer.Add(self.measure_single_histogram_button, flag=wx.EXPAND, border=5)
		
		# Button for opening the Thorlabs camera
		self.__start_camera_label__ = "Start Thorlabs camera"
		self.__stop_camera_label__ = "STOP Thorlabs camera"
		
		self.use_ccd_camera_button = wx.Button (panel, label=self.__start_camera_label__)
		self.Bind (wx.EVT_BUTTON, self.use_ccd_camera, id=self.use_ccd_camera_button.GetId())
		boxsizer.Add(self.use_ccd_camera_button, flag=wx.EXPAND, border=5)
	
		# Button for controlling the moving stages with the joystick
		self.__use_joystick_start_label__ = "Start joystick"
		self.__use_joystick_stop_label__ = "STOP joystick"

		self.use_joystick_button = wx.Button (panel, label=self.__use_joystick_start_label__)
		self.Bind (wx.EVT_BUTTON, self.use_joystick, id=self.use_joystick_button.GetId())
		boxsizer.Add(self.use_joystick_button, flag=wx.EXPAND, border=5)

		
		# Button for saving settings
		self.save_settings_button = wx.Button (panel, label="Save settings...")
		self.Bind (wx.EVT_BUTTON, self.save_settings, id=self.save_settings_button.GetId())
		boxsizer.Add(self.save_settings_button, flag=wx.EXPAND|wx.TOP, border=5)

		# Button for loading settings
		self.load_settings_button = wx.Button (panel, label="Load settings...")
		self.Bind (wx.EVT_BUTTON, self.load_settings, id=self.load_settings_button.GetId())
		boxsizer.Add(self.load_settings_button, flag=wx.EXPAND, border=5)
		
		sizer.Add(boxsizer, pos=(1, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)
		
		######################## Goto panel ###################################
		sb = wx.StaticBox(panel, label="Go to")
		boxsizer = wx.StaticBoxSizer(sb, wx.VERTICAL)
		
		# Button to determined current position
		self.get_current_location_button = wx.Button (panel, label="Update current location")
		self.Bind (wx.EVT_BUTTON, self.get_current_location, self.get_current_location_button)
		boxsizer.Add(self.get_current_location_button, flag=wx.EXPAND, border=5)

		# Coordinates where to move
		boxsizer.Add (wx.StaticText(panel, label="position Ax1 (unit)"), flag=wx.LEFT, border=5)
		self.go_to_position_ax1 = wx.SpinCtrl (panel, value="0", min=-1e4, max=1e4)
		boxsizer.Add (self.go_to_position_ax1, flag=wx.EXPAND, border=5)
		
		boxsizer.Add (wx.StaticText(panel, label="position Ax2 (unit)"), flag=wx.LEFT, border=5)
		self.go_to_position_ax2 = wx.SpinCtrl (panel, value="0", min=-1e4, max=1e4)
		boxsizer.Add (self.go_to_position_ax2, flag=wx.EXPAND, border=5)

		boxsizer.Add (wx.StaticText(panel, label="position Ax3 (unit)"), flag=wx.LEFT, border=5)
		self.go_to_position_ax3 = wx.SpinCtrl (panel, value="0", min=-1e4, max=1e4)
		boxsizer.Add (self.go_to_position_ax3, flag=wx.EXPAND, border=5)

		# Button to start motion
		self.go_to_button = wx.Button (panel, label="Move")
		self.Bind (wx.EVT_BUTTON, self.go_to, self.go_to_button)
		boxsizer.Add(self.go_to_button, flag=wx.EXPAND, border=5)

		# Button to move to points determining the scanning geometry
		self.go_to_points_ABCDE_button = wx.Button (panel, label="Go to point_0")
		self.Bind (wx.EVT_BUTTON, self.go_to_points_ABCDE, self.go_to_points_ABCDE_button)
		boxsizer.Add(self.go_to_points_ABCDE_button, flag=wx.EXPAND, border=5)

		sizer.Add(boxsizer, pos=(2, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		########################### End of constructing panel ######################################
		panel.SetSizer (sizer)
		
		############################# Setting visvis #######################################
		Figure = app.GetFigureClass()
		self.fig = Figure(self)
		
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)
		boxsizer.Add(panel, 1, wx.EXPAND)
		boxsizer.Add(self.fig._widget, 2, wx.EXPAND)

		self.SetSizer (boxsizer)
		self.SetAutoLayout(True)
		self.Layout() 	
		
	def use_joystick (self, event=None, start_camera=True) :
		"""
		<self.use_joystick_button> was clicked 
		"""
		if self.use_joystick_button.GetLabel() == self.__use_joystick_start_label__ :
			self.StopAllJobs ()	
			# Start using camera
			if start_camera : self.use_ccd_camera()
			
			# Begging using the joystick
			self.__use_joystick_event__ = multiprocessing.Event()
			self.__use_joystick_event__.set ()

			self.joystick_process = multiprocessing.Process(target=joystick_control_moving_stage, \
					args=(self.__use_joystick_event__, self.write_serial_port_queue, self.read_serial_port_queue))
			self.joystick_process.start () 
				
			self.use_joystick_button.SetLabel (self.__use_joystick_stop_label__)

		elif self.use_joystick_button.GetLabel() == self.__use_joystick_stop_label__ :
			# Stop using the joystick
			self.__use_joystick_event__.clear()
			self.joystick_process.join ()
			del self.__use_joystick_event__
			del self.joystick_process
			self.use_joystick_button.SetLabel (self.__use_joystick_start_label__)

		else : raise ValueError ("Unrecognized button label") 

	def do_scanning (self, event) :
		"""
		One of scanning or testing buttons was clicked
		"""
		# Extracting which button was clicked
		try :
			button = event.GetEventObject()
			# Mouse double clicking stops scanning
			if event.GetEventType() == wx.wxEVT_LEFT_DCLICK  : button.SetLabel (button.__stop_label__)
		except AttributeError : button = event

		if button.GetLabel() == button.__start_label__ :
			
			self.StopAllJobs ()
			
			if button.command in [CMD_SCAN, CMD_ADAPTIVE_SCAN, CMD_BOUNDARY_SCAN, CMD_BLEACHING_STUDY] :
				# Prepare for scanning the sample
			
				# Initiate Pico Harp
				self.input_pico_harp_queue.put( (CMD_INITIATE_PICO_HARP, self.get_scan_parametrs()) ) 
				result = self.output_pico_harp_queue.get ()
				# Exit, if initialization failed
				if result == RETURN_FAIL : return
				# otherwise, save the resolution
				self.PicoHarp_Resolution.SetValue (result)
			elif button.command in [CMD_TEST_SCANNING, CMD_TAKE_SCANNING_VIDEO] :
				# Turn on camera
				self.use_ccd_camera ()
			
			# Getting all scanning parameters 
			ScanParameters = self.get_scan_parametrs()
			
			# Create list where the disabled objects will be saved
			self.__to_be_enabled__ = []

			# Disable all controls so that no parameters can be changed interactively
			for key in dir(self) :
				obj = getattr (self, key)
				if isinstance(obj, (wx.SpinCtrl, wx.TextCtrl, wx.RadioButton)) and obj.IsEnabled() : 
					obj.Disable (); self.__to_be_enabled__.append (obj) 
			
			# Event indicating that scanning is continuing
			self.scannig_event = multiprocessing.Event()
			self.scannig_event.set ()
			
			# Event for pausing the scanning
			self.pause_scannig_event = multiprocessing.Event()
			self.pause_scannig_event.clear()

			self.scan_sample_process = multiprocessing.Process(target=scan_sample, \
				args=(self.scannig_event, self.pause_scannig_event, \
					self.input_pico_harp_queue,  self.output_pico_harp_queue, self.histogram_buffer, \
					self.write_serial_port_queue, self.read_serial_port_queue, 
					self.input_camera_queue, self.output_camera_queue, self.camera_img_buffer, 
					self.input_shutter_queue, self.output_shutter_queue, ScanParameters, button.command) )
			
			self.scan_sample_process.start ()
			
			# Start timer to monitor weather scanning is over
			TIMER_ID = wx.NewId()
			self.scanning_timer = wx.Timer (self, TIMER_ID)
			self.scanning_timer.Start (2000) # check every 2 seconds
			
			def check_weather_scanning_finished (event) : 
				if not self.scannig_event.is_set () : 
					button.SetLabel (button.__stop_label__); self.do_scanning (button)
			
			wx.EVT_TIMER (self, TIMER_ID, check_weather_scanning_finished)
			
			# Changing the button's label 
			button.SetLabel (button.__pause_label__)

		elif button.GetLabel() == button.__pause_label__ :
		# Pause scanning the sample
			self.pause_scannig_event.set()
			button.SetLabel (button.__resume_label__)

		elif button.GetLabel() == button.__resume_label__ :
		# Resume scanning 
			self.pause_scannig_event.clear()
			button.SetLabel (button.__pause_label__)

		elif button.GetLabel() == button.__stop_label__ : 
		# Stop scanning the sample
			self.scanning_timer.Stop ()
			self.scannig_event.clear ()
			self.pause_scannig_event.clear ()
			self.scan_sample_process.join ()
			del self.scannig_event
			del self.pause_scannig_event
			del self.scan_sample_process
			del self.scanning_timer
		
			# Enable all previously disabled controls except few controls are always disabled
			for obj in self.__to_be_enabled__ : obj.Enable ()
			del self.__to_be_enabled__

			button.SetLabel (button.__start_label__)

		else : raise ValueError ("Unrecognized button label")
		
	
	def draw_scanning_geometry (self, event) :
		"""
		Display scanning geometry
		"""
		self.StopAllJobs()
		visvis.clf()
	
		# Generating scanning geometry
		ScanParameters = self.get_scan_parametrs()
		voltages = get_scanning_polygon (ScanParameters)
		
		# Find axes offset
		ax1_offset = min(voltages, key=operator.itemgetter(0) )[0]
		ax2_offset = min(voltages, key=operator.itemgetter(1) )[1]
		ax3_offset = min(voltages, key=operator.itemgetter(2) )[2]

		# Display actual scanning geometry
		mm_per_V_ax1 = 1e-6*ScanParameters["unit_size"] 
		mm_per_V_ax2 = mm_per_V_ax1
		mm_per_V_ax3 = mm_per_V_ax1
		
		micron_per_V_ax1 = mm_per_V_ax1*1000 
		micron_per_V_ax2 = mm_per_V_ax2*1000
		micron_per_V_ax3 = mm_per_V_ax3*1000

		points = visvis.Pointset(3)
		for ax1, ax2, ax3 in voltages : 
			points.append( (ax1-ax1_offset)*micron_per_V_ax1, (ax2-ax2_offset)*micron_per_V_ax2, (ax3-ax3_offset)*micron_per_V_ax3)
		visvis.plot (points, ms='.', mc='r', mw='9', ls='', mew=0)
	
		# Dsiplay original data points used to generate this geometry
		del points
		points = visvis.Pointset(3)
		for x in self.points :
			points.append ( 1e3*(x[0]-mm_per_V_ax1*ax1_offset), 1e3*(x[1]-mm_per_V_ax2*ax2_offset), 0 )
		
		visvis.plot (points, ms='.', mc='b', mw='9', ls='', mew=0)

		axes = visvis.gca()
		axes.daspectAuto = False
		axes.cameraType = '3d'
		axes.axis.xLabel = 'Ax1 (micron)'
		axes.axis.yLabel = 'Ax2 (micron)'
		axes.axis.zLabel = 'Ax3 (micron)'


	def StopAllJobs (self) :
		"""
		Stop jobs running in other processes
		"""
		if self.use_joystick_button.GetLabel() == self.__use_joystick_stop_label__ : self.use_joystick ()
		if self.use_ccd_camera_button.GetLabel () == self.__stop_camera_label__ : self.use_ccd_camera ()
		if self.measure_single_histogram_button.GetLabel() == self.__stop_measuring_single_histogram_label__ : self.measure_single_histogram()

		for button in [self.bleaching_study_button, self.scanning_button, self.boundary_scanning_button, \
				self.adaptive_scanning_button, self.test_scanning_geometry_button, self.scanning_video_button] :
			if button.GetLabel() in [button.__stop_label__, button.__pause_label__] : self.do_scanning (button)

	def use_ccd_camera (self, event=None) :
		"""
		<self.use_ccd_camera_button> was clicked
		"""
		if self.use_ccd_camera_button.GetLabel () == self.__start_camera_label__ :
		# Start using the camera
			self.StopAllJobs ()

			# Start time to capture image
			TIMER_ID = wx.NewId()
			self.camera_timer = wx.Timer (self, TIMER_ID)
			self.camera_timer.Start (100) # check every 100 ms
		
			def draw_photo (event) :
				self.input_camera_queue.put (CMD_TAKE_PHOTO)
				if self.output_camera_queue.get () == RETURN_FAIL : return	
				# Displaying image
				self.camera_img_buffer.acquire ()
				try : 
					# Updating data, if <self.texture> was initiated
					img = camera_img_shape(self.camera_img_buffer)
					self.texture.SetData (img)
					self.texture.SetClim (img.min(), img.max())
					self.texture.Refresh ()
				except AttributeError :
					# Initiate <self.texture> and plot the first figure
					visvis.clf()
					self.texture = visvis.imshow (camera_img_shape(self.camera_img_buffer))
					visvis.title ('Camera view')
					ax = visvis.gca()
					ax.axis.xTicks = []
					ax.axis.yTicks = []
				self.camera_img_buffer.release ()
				
			wx.EVT_TIMER (self, TIMER_ID, draw_photo)

			# Changing button's label
			self.use_ccd_camera_button.SetLabel (self.__stop_camera_label__)

		elif self.use_ccd_camera_button.GetLabel () == self.__stop_camera_label__ :
		# Stop using the camera
		
			self.camera_timer.Stop ()
			del self.camera_timer

			# Cleaning the photo
			del self.texture
			visvis.clf(); self.fig.DrawNow()

			self.use_ccd_camera_button.SetLabel (self.__start_camera_label__)

		else : raise ValueError ("Unrecognized button label")


	def measure_single_histogram (self, event=None) :
		"""
		<self.measure_single_histogram_button> was clicked
		"""
		if self.measure_single_histogram_button.GetLabel() == self.__start_measuring_single_histogram_label__ : 
			self.StopAllJobs ()
			
			# Initiate Pico Harp
			self.input_pico_harp_queue.put( (CMD_INITIATE_PICO_HARP, self.get_scan_parametrs()) ) 
			result = self.output_pico_harp_queue.get ()
			if result == RETURN_FAIL : return
			self.PicoHarp_Resolution.SetValue (result)

			# Begin measuring histogram before we start the timer 
			self.input_pico_harp_queue.put (CMD_MEASURE_HISTOGRAM)
			
			# Clearing the figure
			visvis.clf()

			def draw_single_histogram (event) :
				" Timer function "
				if self.output_pico_harp_queue.get () == RETURN_FAIL : 
					# Asking to measure histogram one more time 
					self.input_pico_harp_queue.put (CMD_MEASURE_HISTOGRAM)
					return
			
				# Extract Pico Harp resolution
				self.input_pico_harp_queue.put (CMD_GET_RESOLUTION_PICO_HARP)			
				result = self.output_pico_harp_queue.get ()
				if result <> RETURN_FAIL : self.PicoHarp_Resolution.SetValue (result) 
	
				self.histogram_buffer.acquire()
				histogram = histogram_shape (self.histogram_buffer) 
				# time axis in nano sec
				time_ns = 1.e-3 * self.PicoHarp_Resolution.GetValue() * np.arange(histogram.size)
		
				ax = visvis.gca()
				ax.Clear()
			
				histogram_total_count = histogram.sum() 
				if histogram_total_count == 0 : 
					# there is no histogram 
					visvis.plot (time_ns, histogram, mw=0.1, ms='.')
					ax.axis.yLabel = 'Histogram, counts'
				else :
					data = np.log10(histogram)
					visvis.plot (time_ns, data, mw=0.1, ms='.')
					ax.axis.yLabel = 'Histogram, log(counts)'
				
				self.histogram_buffer.release ()
				
				ax.axis.xLabel = 'time (ns)'
							
				# Extract count rates PicoHarp
				self.input_pico_harp_queue.put (CMD_GET_COUNT_RATES_PICO_HARP)
				result = self.output_pico_harp_queue.get ()
				if result <> RETURN_FAIL :
					result += (histogram_total_count,)
					visvis.title ("Count rates: Ch0 %.2e / Ch1 %.2e / histogram %.2e" % result )
				else : visvis.title ("Histogram")
			
				# Measuring a histogram, which is to be displayed at the next call
				self.input_pico_harp_queue.put (CMD_MEASURE_HISTOGRAM)

		
			# Set up the timer
			TIMER_ID = wx.NewId()
			self.histogram_timer = wx.Timer (self, TIMER_ID)
			self.histogram_timer.Start (self.PicoHarp_Tacq.GetValue()) 
			wx.EVT_TIMER (self, TIMER_ID, draw_single_histogram)
			
			# Changing button's label
			self.measure_single_histogram_button.SetLabel (self.__stop_measuring_single_histogram_label__)

		elif self.measure_single_histogram_button.GetLabel() == self.__stop_measuring_single_histogram_label__ :
			##### Stop measuring single histograms #####
			self.histogram_timer.Stop ()
			del self.histogram_timer
			# Cleaning the pico harp queue 
			self.output_pico_harp_queue.get ()
			# Changing the button's label
			self.measure_single_histogram_button.SetLabel (self.__start_measuring_single_histogram_label__)

		else : raise ValueError ("Unrecognized button label")

	def get_current_location (self, event=None) :
		"""
		<self.get_current_location_button> was clicked
		"""
		# Save the current position New Port
		self.write_serial_port_queue.put( (True, "1WS0;1TP\r") )
		position_ax1 = float(self.read_serial_port_queue.get())
		
		self.write_serial_port_queue.put( (True, "2WS0;2TP\r") )
		position_ax2 = float(self.read_serial_port_queue.get())
		
		self.write_serial_port_queue.put( (True, "3WS0;3TP\r") )
		position_ax3 = float(self.read_serial_port_queue.get())

		# Converting mm to units
		unit_size = self.unit_size.GetValue()
		volt_ax1 = int(round( 1.e6*position_ax1 / unit_size ))
		volt_ax2 = int(round( 1.e6*position_ax2 / unit_size ))
		volt_ax3 = int(round( 1.e6*position_ax3 / unit_size ))

		# Setting the values
		self.go_to_position_ax1.SetValue(volt_ax1)
		self.go_to_position_ax2.SetValue(volt_ax2)
		self.go_to_position_ax3.SetValue(volt_ax3)

	def go_to_points_ABCDE (self, event) :
		"""
		<self.go_to_points_ABCDE> button was clicked. Move to one of points defying the scanning geometry
		"""
		# If there are no points, ignor this call
		if len(self.points) == 0 : return

		# Determine to which position to move
		position = int(event.GetEventObject().GetLabel().split('_')[-1])
		ax1, ax2 = self.points[position]

		# Update button's label (to move to the next point in the list)
		position += 1
		if position >= len(self.points) : event.GetEventObject().SetLabel("Go to point_0")
		else : event.GetEventObject().SetLabel("Go to point_%d" % position)
		
		# Moving to chosen point
		self.write_serial_port_queue.put( (False, "1PA%.9e\r" % ax1) )
		self.write_serial_port_queue.put( (False, "2PA%.9e\r" % ax2) )

	def go_to (self, event) :
		"""
		<self.go_to_button> was clicked
		"""
		new_position = ( self.go_to_position_ax1.GetValue(), self.go_to_position_ax2.GetValue(), self.go_to_position_ax3.GetValue() )
		moveto (new_position, self.write_serial_port_queue, self.read_serial_port_queue, self.unit_size.GetValue(), 
				self.input_shutter_queue, self.output_shutter_queue, False)

	def add_point (self, event) :
		"""
		<self.add_point_button> was clicked to characterize the scanning geometry.
		"""
		# Saving the current position of the New Port moving stage
		self.write_serial_port_queue.put( (True, "1TP\r") ); ax1 = float(self.read_serial_port_queue.get())
		self.write_serial_port_queue.put( (True, "2TP\r") ); ax2 = float(self.read_serial_port_queue.get())

		# Adding the point
		self.points.append( (ax1, ax2) )

		# Save a photo of that place
		self.input_camera_queue.put (CMD_TAKE_PHOTO_ACCURATELY)
		if self.output_camera_queue.get () <> RETURN_FAIL : 	
			self.camera_img_buffer.acquire ()
			# Saving the photo
			setattr(self, "photo_%d" % (len(self.points)-1), camera_img_shape(self.camera_img_buffer) ) 
			self.camera_img_buffer.release ()

		# Changing the colour
		self.add_point_button.SetBackgroundColour('red')
		
		# Update the values of number of scanning steps
		# based on the 2.5micron resolution 
		resolution = 2.5e-3
		
		min_ax1 = min(self.points, key=operator.itemgetter(0))[0]; max_ax1 = max(self.points, key=operator.itemgetter(0))[0] 
		min_ax2 = min(self.points, key=operator.itemgetter(1))[1]; max_ax2 = max(self.points, key=operator.itemgetter(1))[1]
		
		self.number_steps_ax1.SetValue( ( max_ax1 - min_ax1 )/resolution ) 
		self.number_steps_ax2.SetValue( ( max_ax2 - min_ax2 )/resolution ) 
		
	def get_scan_parametrs (self) :
		"""
		Return dictionay of parameters needed for scanning the sample.
		This function is closely related to <self.load_settings>
		"""
		ScanParameters = {
			# Scanning settings 
			"number_steps_ax1"	: self.number_steps_ax1.GetValue(),
			"number_steps_ax2" 	: self.number_steps_ax2.GetValue(),
			"number_steps_ax3" 	: self.number_steps_ax3.GetValue(),
			"filename"		: str(self.scan_filename.GetValue()),
			"info"			: str(self.scan_info.GetValue()),
			"points"		: self.points,
			"piezo_min_volt_ax3"	: self.piezo_min_volt_ax3.GetValue(),
			"piezo_max_volt_ax3"	: self.piezo_max_volt_ax3.GetValue(),
			"histogram_counts_threshold"		: self.histogram_counts_threshold.GetValue(),
			"requested_histogram_counts"		: self.requested_histogram_counts.GetValue(),
			"revisit_number" 			: self.revisit_number.GetValue(),
			"adaptive_scan_point_weight" 		: self.adaptive_scan_point_weight.GetValue(),
			"distance_between_adaptive_points"	: self.distance_between_adaptive_points.GetValue(),
			"unit_size" 				: self.unit_size.GetValue (),
			# Pico Harp settings
			"Offset"	: self.PicoHarp_Offset.GetValue(),
			"CFDZeroX0"	: self.PicoHarp_CFDZeroX0.GetValue(),
			"CFDLevel0"	: self.PicoHarp_CFDLevel0.GetValue(),
			"CFDZeroX1" 	: self.PicoHarp_CFDZeroX1.GetValue(),
			"CFDLevel1" 	: self.PicoHarp_CFDLevel1.GetValue(),
			"SyncDiv" 	: self.PicoHarp_SyncDiv.GetValue(),
			"Range"		: self.PicoHarp_Range.GetValue(),
			"Tacq"		: self.PicoHarp_Tacq.GetValue(),
			# Pico Harp resolution obtained after initializing 
			"Resolution"	: self.PicoHarp_Resolution.GetValue ()
		}

		# Adding photos taken while defying scanning geometry
		for n in xrange(len(self.points)) :
			key = "photo_%d" % n
			ScanParameters[key] = getattr(self, key)
		
		# Determine which axis is fast varying
		if self.ax1_slow.GetValue() 	: ScanParameters["slow_varying_axis"] = 0 
		elif self.ax2_slow.GetValue() 	: ScanParameters["slow_varying_axis"] = 1 
		elif self.ax3_slow.GetValue() 	: ScanParameters["slow_varying_axis"] = 2
		elif self.random_ax_slow.GetValue(): ScanParameters["slow_varying_axis"] = 3
		else : raise ValueError ("Not defined axis")

		return ScanParameters


	def load_settings (self, event) :
		"""
		<self.load_settings_button> was clicked. This function is closely related to <self.get_scan_parametrs>
		"""
		openFileDialog = wx.FileDialog(self, "Open HDF5 file to load settings", "", "",
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
		# Check whether user canceled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return	

		# Opening file
		file_settings = h5py.File (openFileDialog.GetPath(), 'r')

		for key, data in file_settings["parameters"].iteritems () :
			if key == "number_steps_ax1" 	: self.number_steps_ax1.SetValue(data[...])
			elif key == "number_steps_ax2" 	: self.number_steps_ax2.SetValue(data[...])
			elif key == "number_steps_ax3" 	: self.number_steps_ax3.SetValue(data[...])
			elif key == "filename"		: self.scan_filename.SetValue(str(data[...]))
			elif key == "info"			: self.scan_info.SetValue(str(data[...]))
			elif key == "histogram_counts_threshold"	: self.histogram_counts_threshold.SetValue(data[...])
			elif key == "requested_histogram_counts"		: self.requested_histogram_counts.SetValue(data[...])
			elif key == "revisit_number" 	: self.revisit_number.SetValue(data[...]) 
			elif key == "adaptive_scan_point_weight" 	: self.adaptive_scan_point_weight.SetValue(data[...])
			elif key == "distance_between_adaptive_points"	: self.distance_between_adaptive_points.SetValue(data[...])
			elif key == "points"		: self.points = list(data[...]); self.add_point_button.SetBackgroundColour('red')
			elif "photo_" in key		: setattr(self, key, data[...]) # Loading photos that correspond to <self.points>
			elif key == "piezo_min_volt_ax3" : self.piezo_min_volt_ax3.SetValue(data[...])
			elif key == "piezo_max_volt_ax3" : self.piezo_max_volt_ax3.SetValue(data[...])
			elif key == "unit_size" : 	self.unit_size.SetValue (data[...])
			# Pico Harp settings
			elif key == "Offset"		: self.PicoHarp_Offset.SetValue(data[...])
			elif key == "CFDZeroX0"		: self.PicoHarp_CFDZeroX0.SetValue(data[...])
			elif key == "CFDLevel0"		: self.PicoHarp_CFDLevel0.SetValue(data[...])
			elif key == "CFDZeroX1" 	: self.PicoHarp_CFDZeroX1.SetValue(data[...])
			elif key == "CFDLevel1" 	: self.PicoHarp_CFDLevel1.SetValue(data[...])
			elif key == "SyncDiv"	 	: self.PicoHarp_SyncDiv.SetValue(data[...])
			elif key == "Range"		: self.PicoHarp_Range.SetValue(data[...])
			elif key == "Tacq"		: self.PicoHarp_Tacq.SetValue(data[...])
			# Pico Harp resolution obtained after initializing 
			elif key == "Resolution"	: self.PicoHarp_Resolution.SetValue(data[...])
			# Vast varying axis
			elif key == "slow_varying_axis" :
				if data[...] == 0 : self.ax1_slow.SetValue (True)
				elif data[...] == 1 : self.ax2_slow.SetValue (True)
				elif data[...] == 2 : self.ax3_slow.SetValue (True)
				elif data[...] == 3 : self.random_ax_slow.Setvalue (True) 

		file_settings.close ()

	def save_settings (self, event) :
		"""
		<self.save_settings_button> was clicked
		"""
		default_filename = "settings_" + str(self.scan_filename.GetValue())
		openFileDialog = wx.FileDialog(self, "Open HDF5 file to load settings", "", default_filename,
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
		# Check whether user canceled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return	

		# Saving settings
		file_settings = h5py.File (openFileDialog.GetPath(), 'w')
		parameters_grp = file_settings.create_group("parameters")
		
		for key, value in self.get_scan_parametrs().items() :
			parameters_grp[key] = value
		
		file_settings.close ()
		

#########################################################################
if __name__ == '__main__' :
	app.Create()
	CMicroscopeConrol (None, title='FRET meter')
	app.Run()

