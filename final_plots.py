import numpy as np
import cPickle
import glob
import operator
from viewer import get_time_life_ns
from sklearn.cluster import k_means
import matplotlib.pyplot as plt
from ntpath import basename
import multiprocessing 

# Defying operators for fast extractions
get_fluoresence_rate 	= operator.itemgetter("fluoresence_rate")
get_fitted_data			= operator.itemgetter("fitted_data")
	
# Files containing information about working sensor 
vinculinTS_filenames = glob.glob("Data_Final_Plots\\vinculinTS\\vinculinTS_*_direct_sum_boundary.flim")
# Files for tailless control
TL_filenames = glob.glob("Data_Final_Plots\\TL\\TL_*_direct_sum_boundary.flim")
# Samples with no constructions
blank_eggs_filenames = glob.glob("Data_Final_Plots\\blank\\blank_*_direct_sum_boundary.flim")


def stat_from_flim_file (filename) :
	"""
	Extract the statistical characteristics of 
	fluorescence rates and time life form a FLIM file
	"""
	# loading data
	with open(filename, "rb") as in_file : fitted_data = cPickle.load(in_file)
	
	# Extract the measurements
	fluoresence_rate, extracted_fitted_data \
		= zip( *( (get_fluoresence_rate(value), get_fitted_data(value)) for value in fitted_data["points"].itervalues()) )
	
	fluoresence_rate	= np.array(fluoresence_rate)
	time_life			= get_time_life_ns(np.array(extracted_fitted_data), Resolution=fitted_data["Resolution"])
	
	return np.mean(fluoresence_rate), np.std(fluoresence_rate), np.mean(time_life), np.std(time_life)

def plot_cumulative_stat () :
	"""
	Display cumulative statistics (calculated by function <stat_from_flim_file>) for all the files.
	"""
	vinculinTS_mean_fluoresence_rate, vinculinTS_std_fluoresence_rate, vinculinTS_mean_time_life, vinculinTS_std_time_life \
		= zip( *map(stat_from_flim_file, vinculinTS_filenames) )


	TL_mean_fluoresence_rate, TL_std_fluoresence_rate, TL_mean_time_life, TL_std_time_life \
		= zip( *map(stat_from_flim_file, TL_filenames) )
		
	#plt.errorbar(vinculinTS_mean_fluoresence_rate, vinculinTS_mean_time_life, 
	#	xerr=vinculinTS_std_fluoresence_rate, yerr=vinculinTS_std_time_life, fmt='o')
		
	#plt.errorbar(TL_mean_fluoresence_rate, TL_mean_time_life, 
	#	xerr=TL_std_fluoresence_rate, yerr=TL_std_time_life, fmt='ro')
		
	plt.errorbar(vinculinTS_mean_fluoresence_rate, vinculinTS_mean_time_life, fmt='ro')
	plt.errorbar(TL_mean_fluoresence_rate, TL_mean_time_life, fmt='gd')

	plt.xlabel("Rate of fluorescence (counts per millisecond)")
	plt.ylabel("Time life (nanosecond)")

	plt.legend(["vinculinTS", "TL"])

	print "vinculinTS life-time %.2e +/- %.2e (ns) " % ( np.mean(vinculinTS_mean_time_life), np.std(vinculinTS_mean_time_life) )
	print "vinculinTS fluorescence rate %.2e +/- %.2e (counts per millisecond)" % ( np.mean(vinculinTS_mean_fluoresence_rate), np.std(vinculinTS_mean_fluoresence_rate) )
	print "TL life-time %.2e +/- %.2e (ns) " % ( np.mean(TL_mean_time_life), np.std(TL_mean_time_life) )
	print "TL fluorescence rate %.2e +/- %.2e (counts per millisecond)" % ( np.mean(TL_mean_fluoresence_rate), np.std(TL_mean_fluoresence_rate) )

	plt.show()
		
def sliced_sample (filename, nchunks = 15) :
	"""
	Slice a sample saved in file <filename> into <nchunks> chunks. 
	"""
	# loading data
	with open(filename, "rb") as in_file : fitted_data = cPickle.load(in_file)
	
	# Extract the measurements
	coordinates, fluoresence_rate, extracted_fitted_data = zip( *( 
		(key, get_fluoresence_rate(value), get_fitted_data(value)) for key,value in fitted_data["points"].iteritems()
	) )
	
	# Converting the data into numpy arrays 
	coordinates 			= np.array(coordinates, dtype=np.float)
	fluoresence_rate 		= np.array(fluoresence_rate)
	extracted_fitted_data	= np.array(extracted_fitted_data)
	
	# Get life-time from fitted data
	time_life = get_time_life_ns(extracted_fitted_data, Resolution=fitted_data["Resolution"])
	
	# Cluster coordinates
	_, chunk_labels, _ = k_means(coordinates, nchunks)
	
	# Display clustering for debugging purposes
	# plt.scatter ( coordinates[:,0], coordinates[:,1], c=chunk_labels)
	# plt.show()

	# Extract the incidences of each slice
	chunks_indx = ( np.nonzero(chunk_labels==chunk) for chunk in xrange(nchunks) )
	
	# Calculate the list of tuples of mean fluorescence rate and life-time for each slice
	return [ ( np.mean(fluoresence_rate[indx]), np.mean(time_life[indx]) ) for indx in chunks_indx ]
	
def plot_sliced_sample (filename, fig_folder="Data_Final_Plots\\Plots\\") :
	"""
	Generate plot on the basis of information obtained from function <sliced_sample>.
	<fig_folder> is a folder where the figure will be saved.
	"""
	# Extract fluorescence rate and time-life for different parts of the sample
	# saved in the file <filename>
	fluoresence_rate, time_life = zip(*sliced_sample(filename)) 
	
	# File name of the figure where the plot will be saved
	fig_filename = "%s %s.png" % tuple(basename(filename).split("_")[:2])
	
	# Generating the corresponding plot
	plt.clf()
	plt.plot (fluoresence_rate, time_life, 'ob')
	plt.title( "Chopped sample " + fig_filename)
	plt.xlabel("Rate of fluorescence (counts per millisecond)")
	plt.ylabel("Time life (nanosecond)")
	
	plt.savefig( fig_folder + fig_filename)

def plot_cumulative_sliced_samples () :
	"""
	
	"""
	pool = multiprocessing.Pool()
	
	# Chop each file and collect all the data
	TL_mean_fluoresence_rate, TL_mean_time_life \
			= zip(*reduce( operator.add, pool.map(sliced_sample, TL_filenames) ))
	vinculinTS_mean_fluoresence_rate, vinculinTS_mean_time_life \
			= zip(*reduce( operator.add, pool.map(sliced_sample, vinculinTS_filenames) ))
	blank_egg_mean_fluoresence_rate, blank_egg_mean_time_life \
			= zip(*reduce( operator.add, pool.map(sliced_sample, blank_eggs_filenames) ))
	
	print "vinculinTS life-time %.2e +/- %.2e (ns) " % ( np.mean(vinculinTS_mean_time_life), np.std(vinculinTS_mean_time_life) )
	print "vinculinTS fluorescence rate %.2e +/- %.2e (counts per millisecond)" % ( np.mean(vinculinTS_mean_fluoresence_rate), np.std(vinculinTS_mean_fluoresence_rate) )
	print "TL life-time %.2e +/- %.2e (ns) " % ( np.mean(TL_mean_time_life), np.std(TL_mean_time_life) )
	print "TL fluorescence rate %.2e +/- %.2e (counts per millisecond)" % ( np.mean(TL_mean_fluoresence_rate), np.std(TL_mean_fluoresence_rate) )
	
	plt.clf()
	
	plt.plot (vinculinTS_mean_fluoresence_rate, vinculinTS_mean_time_life, 'ro')
	plt.plot (TL_mean_fluoresence_rate, 		TL_mean_time_life, 'gd')
	plt.plot (blank_egg_mean_fluoresence_rate, 	blank_egg_mean_time_life, 'k*')
		
	plt.legend(['vinculinTS', 'TL', 'blank sample'])
	
	plt.title ('Choped samples and mearged')
	
	plt.xlabel("Rate of fluorescence (counts per millisecond)")
	plt.ylabel("Time life (nanosecond)")
	
	plt.show()
	
if __name__ == '__main__' :
	
	#plot_cumulative_stat()
	plot_cumulative_sliced_samples()
	
	# Chop individual sample and save each figure separately
	# pool = multiprocessing.Pool()
	#pool.map(plot_sliced_sample, vinculinTS_filenames + TL_filenames)
