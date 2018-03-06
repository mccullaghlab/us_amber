##!/Users/martinmccullagh/anaconda/bin/python
# NOTE: will have to point to version of python and have MDAnalysis library

# USAGE: dist_windows.py [config file]
# CONFIG FILE FORMAT:
# psffile = [psf, prmtop, gro, or pdb file]
# dcdfile = [traj file in format trr, dcd, etc]
# atom_sel_1 = [CHARMM style atom selection for first group of atoms]
# atom_sel_2 = [CHARMM style atom selection for second group of atoms]
# dist_max = [maximum distance]
# dist_min = [minimum distance]
# dist_delta = [distance bin size]


# load libraries
import sys
import os
import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



#################################################################################################################
##############################################     SUBROUTINEs     ##############################################
#################################################################################################################


# define subroutines
def ParseConfigFile(cfg_file):
	global dist_min,dist_max,dist_delta, output_file_name
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='dist_min':
				dist_min = float(value)
			elif option.lower()=='dist_max':
				dist_max = float(value)
			elif option.lower()=='dist_delta':
				dist_delta = float(value)
			elif option.lower()=='hist_out':
				output_file_name = value
			else :
				print "Option:", option, " is not recognized"
	
	f.close()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#################################################################################################################
##############################################     MAIN PROGRAM   ##############################################
#################################################################################################################

ParseConfigFile(sys.argv[1])

n_dist_bins = int((dist_max-dist_min)/dist_delta)

count = 0
window_value = []
for dirs, subdirs, files in os.walk('../'):
	if dirs[3:9] == "window" and is_number(dirs[9:]):
		window = float(dirs[9:])
		window_string = dirs[9:]
		if count == 0 or window < window_min:
			window_min = window
		if count == 0 or window > window_max:
			window_max = window
		check_file = dirs + "/ben2.run.r" + window_string + ".dat"
		if os.path.isfile(check_file):	
			# open file and create histogram
			data = np.loadtxt(check_file)
			if count == 0:
				window_data = np.histogram(data[:,1],bins=n_dist_bins,range=(dist_min,dist_max),density=True)[0]
			else:
				window_data = np.column_stack((window_data,np.histogram(data[:,1],bins=n_dist_bins,range=(dist_min,dist_max),density=True)[0]))
			count += 1
			window_value.append(window)


print window_data.shape
#print window_data
out = open(output_file_name,"w")
x_data = np.empty(window_data.shape[0],dtype=float)
for i in range(window_data.shape[0]):
	x_data[i] = i*dist_delta+dist_min
	out.write("%10.5f " % (i*dist_delta+dist_min))
	for j in range(window_data.shape[1]):
		out.write("%10.5f " % (window_data[i,j]))
	out.write("\n")
out.close()		

for i in range(window_data.shape[1]):
        plt.plot(x_data, window_data[:,i],lw=2)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel("distance", size=12)
plt.ylabel("probability density", size=12)
plt.savefig('hists.png')
plt.close()

