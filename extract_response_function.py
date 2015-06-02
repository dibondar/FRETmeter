"""
Convert file containing histograms into the response function
"""
import h5py
import wx
import numpy as np
import matplotlib.pyplot as plt

#############################################################################
# Select the file cantoning histograms, 
# which will be converted to response function
app = wx.App ()
openFileDialog = wx.FileDialog (None, "Chose file containing histograms to get repose function", "", "", \
		"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

# Check whether user canceled
if openFileDialog.ShowModal() == wx.ID_CANCEL : 
	raise ValueError("HDF5 file is not selected")
	
hist_filename = openFileDialog.GetPath()

del app
#############################################################################

# Adding up all histograms
with h5py.File(hist_filename, 'r') as F :
	for histogram in F["histograms"].values() :
		try :
			sum_histogram += histogram[...]
		except NameError :
			sum_histogram = histogram[...]

	# Loading resolution
	Resolution = F["parameters/Resolution"][...]
	
#############################################################################
# Remove zeros 
cut_off = np.nonzero(sum_histogram)[0].max()
sum_histogram = sum_histogram[:cut_off]

#############################################################################
# Normalize to max 1
sum_histogram = sum_histogram.astype(np.float) 
sum_histogram /= sum_histogram.max()

# Delete the background
sum_histogram[ np.nonzero(sum_histogram < 0.03) ] = 0

#############################################################################
# plot
plt.plot( 1e-3*Resolution*np.arange(sum_histogram.size), sum_histogram )
plt.title ("Total histogram with resolution %d ps" % Resolution )
plt.xlabel("time (ns)")
plt.ylabel("counts")
plt.show()

# Save
sum_histogram.tofile("response_function_%dps.dat" % Resolution)

