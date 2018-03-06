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
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
import numpy.linalg


#################################################################################################################
##############################################     SUBROUTINEs     ##############################################
#################################################################################################################


# define subroutines
def ParseConfigFile(cfg_file):
	global inp_psf_file,inp_dcd_file,user_sel_1,user_sel_2,user_sel_3,user_sel_4,dist_min,dist_max,dist_delta, job_script_file_name, force_constant
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
			if option.lower()=='psffile':
				inp_psf_file = value
			elif option.lower()=='dcdfile':
				inp_dcd_file = value
			elif option.lower()=='jobfile':
				job_script_file_name = value
			elif option.lower()=='atom_sel_1':
				user_sel_1 = value
			elif option.lower()=='atom_sel_2':
				user_sel_2 = value
			elif option.lower()=='dist_min':
				dist_min = float(value)
			elif option.lower()=='k':
				force_constant = float(value)
			elif option.lower()=='dist_max':
				dist_max = float(value)
			elif option.lower()=='dist_delta':
				dist_delta = float(value)
			else :
				print "Option:", option, " is not recognized"
	
	f.close()

# compute distance between two poiunts in R3
def computeDist(r1,r2):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

# compute distance between two poiunts in R3 using PBC
def computePbcDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

def computeDih(r1,r2,r3,r4):

	# define the distance vectors
	b1 = r1 - r2
	b2 = r2 - r3
	b3 = r3 - r4
	# compute the cross product vectors
	A = numpy.cross(b1,b2)
	A = A/math.sqrt(numpy.dot(A,A))
	B = numpy.cross(b2,b3)
	B = B/math.sqrt(numpy.dot(B,B))
	C = numpy.cross(b2,A)
	C = C/math.sqrt(numpy.dot(C,C))
	# now compute the dihedral
	dih = -numpy.arctan2(numpy.dot(C,B),numpy.dot(A,B))

	# convert to degrees
	dih *= 180.0/3.1415926535
	return dih

def wrapPositions(old_positions,n_atoms,box,target):
	new_positions = numpy.empty( (n_atoms,3), dtype=float)
	for i in range(0,n_atoms):
		for j in range(0,3):
			if (old_positions[i,j]-target[j]) < -box[j]/2.0:
				#print "Changing position of atom", i+1
				new_positions[i,j] = old_positions[i,j] + box[j]
			elif (old_positions[i,j]-r1[j]) > box[j]/2.0:
				#print "Changing position of atom", i+1
				new_positions[i,j] = old_positions[i,j] - box[j]
			else:
				new_positions[i,j] = old_positions[i,j]
	# update the MDAnalysis positions
	return new_positions

def writeInpcrd(positions,filename,box):

	n_atoms = positions.shape[0]
	out = open(filename,'w')
	# first write a blank line
	out.write("\n")
	out.write("%6d\n" % (n_atoms))
	count = 0
	# write atom positions
	for i in range(n_atoms):
		for j in range(3):
			out.write("%12.7f" % (positions[i,j]))
			count += 1
			if (count%6==0):
				out.write("\n")
	if (count%6 != 0):
		out.write("\n")
	# write box info
	for j in range(3):
			out.write("%12.7f" % (box[j]))
	for j in range(3):
			out.write("%12.7f" % (90.0))
	out.close()



def writeJobScriptHeader(job_script,job_script_file_name):

	file_root = os.path.splitext(job_script_file_name)[0]
	job_script.write("#!/bin/bash\n") 
	job_script.write("#SBATCH --job-name=%s\n" % (file_root))
	job_script.write("#SBATCH --output=%s.out\n" % (file_root))
	job_script.write("#SBATCH --time=48:00:00\n")
	job_script.write("#SBATCH --nodes=1\n")
	job_script.write("#SBATCH --partition=mccullagh-gpu\n")
	job_script.write("#SBATCH --reservation=workshop_gpu\n")
	job_script.write("#SBATCH --gres=gpu:titan:1\n")

	job_script.write("export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/software/usr/gcc-4.9.2/lib64\"\n")
	job_script.write("export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/software/usr/hpcx-v1.2.0-292-gcc-MLNX_OFED_LINUX-2.4-1.0.0-redhat6.6/ompi-mellanox-v1.8/lib\"\n")
	job_script.write("export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/usr/local/cuda-7.5/lib64\"\n")
	job_script.write("export AMBERHOME=\"/mnt/lustre_fs/users/mjmcc/apps/amber14\"\n")
	job_script.write("export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:$AMBERHOME/lib\"\n")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#################################################################################################################
##############################################     MAIN PROGRAM   ##############################################
#################################################################################################################

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "PSF file:", inp_psf_file
print "Coord DCD file:", inp_dcd_file

# start MDAnalysis with a universal
#coordinate universe
coord = MDAnalysis.Universe(inp_psf_file,inp_dcd_file)

# print some general log info
print "Numer of time steps in coordinate trajectory:", len(coord.trajectory)

# define the first selection
sel1 = coord.select_atoms(user_sel_1)
print "First atom selection:", user_sel_1
print sel1
# define the second atom selection
sel2 = coord.select_atoms(user_sel_2)
print "Second atom selection:", user_sel_2
print sel2

# allocate COM position vectors
r1 = numpy.empty(3,dtype=float)
r2 = numpy.empty(3,dtype=float)

# define ATOMGROUP that contains all atoms for printing purposes
all_atoms = coord.select_atoms("all")

# determine number of distance bins
n_dist_bins = int((dist_max-dist_min)/dist_delta)
dist_pop = numpy.zeros(n_dist_bins,dtype=int)

for ts in coord.trajectory:
	# get box dimensions
	box = coord.trajectory.ts.dimensions[:3]
	# compute new centers of mass 
	r1 = sel1.center_of_mass()
	r2 = sel2.center_of_mass()
	# compute distance
	dist = computeDist(r1,r2)
	# determine distance bin
	dist_bin = int((dist-dist_min)/dist_delta)
	# check to see if it falls in range and if we have already populated that window
	if dist_bin >=0 and dist_bin < n_dist_bins and dist_pop[dist_bin]==0:
		dist_pop[dist_bin]=1
		window = dist_bin*dist_delta + dist_min
		window_dir = "window"+str(window)
		command = "mkdir " + window_dir
		os.system(command)
		window_pdb_file = window_dir + "/window"+str(window)+".inpcrd"
		print "Writing pdb file:", window_pdb_file, "with a distance of:", dist
		writeInpcrd(all_atoms.positions,window_pdb_file,box)


# create job submit script and copy necessary files
job_script = open(job_script_file_name, "w")
writeJobScriptHeader(job_script, job_script_file_name)
# loop through bins
for dist_bin in range(n_dist_bins):
	if dist_pop[dist_bin] > 0:
		window = dist_bin*dist_delta + dist_min
		window_dir = "window"+str(window)
		window_inpcrd_file = "window"+str(window)+".inpcrd"
		window_equil_rst_file = "window"+str(window)+".equil.rst"
		window_equil_log_file = "window"+str(window)+".equil.log"
		window_equil_traj_file = "window"+str(window)+".equil.ncdf"
		window_run_rst_file = "window"+str(window)+".run.rst"
		window_run_log_file = "window"+str(window)+".run.log"
		window_run_traj_file = "window"+str(window)+".run.ncdf"
		# write to job script
		job_script.write("# run window %s\n" % (str(window)))
		job_script.write("cd %s\n" % (window_dir))
		job_script.write("$AMBERHOME/bin/pmemd.cuda -O -i ben2.umb.equil.in -o %s -p ben2.prmtop -c %s -r %s -x %s\n" % (window_equil_log_file, window_inpcrd_file, window_equil_rst_file, window_equil_traj_file))
		job_script.write("$AMBERHOME/bin/pmemd.cuda -O -i ben2.umb.run.in -o %s -p ben2.prmtop -c %s -r %s -x %s\n" % (window_run_log_file, window_equil_rst_file, window_run_rst_file, window_run_traj_file))
		job_script.write("cd ..\n")
		# copy files
		command = "cp window_files/ben2.prmtop " + window_dir
		os.system(command)
		command = "sed -e s/XX/" + str(window) + "/g < window_files/ben2.umb.equil.in > " + window_dir + "/ben2.umb.equil.in"
		os.system(command)
		command = "sed -e s/XX/" + str(window) + "/g < window_files/ben2.umb.run.in > " + window_dir + "/ben2.umb.run.in"
		os.system(command)
		command = "sed -e s/XX/" + str(window) + "/g -e s/kk/" + str(force_constant/2.0) + "/g < window_files/dist.rst > " + window_dir + "/dist_r" + str(window) + ".rst"
		os.system(command)
#
job_script.write("/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python hists.py hist.cfg\n")
# create wham directory
command = "mkdir wham"
os.system(command)
wham_file = "wham/wham.config"
wham = open(wham_file, 'w')
count = 0
for dirs, subdirs, files in os.walk('./'):
	if dirs[2:8] == "window" and is_number(dirs[8:]):
		window = float(dirs[8:])
		window_string = dirs[8:]
		if count == 0 or window < window_min:
			window_min = window
		if count == 0 or window > window_max:
			window_max = window
		check_file = dirs + "/ben2.run.r" + window_string + ".dat"
		print check_file
#		if os.path.isfile(check_file):	
		wham.write("../window%s/ben2.run.r%s.dat %8.3f %8.3f\n" %(window_string, window_string, window, force_constant))
		count += 1
wham.close()
job_script.write("cd wham\n")
job_script.write("wham %f %f 100 0.0001 298.0 0 wham.config wham.out 500 12345\n" % (window_min, window_max))
job_script.write("/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python plot_wham.py\n")
job_script.close()
	
command = "sbatch " + job_script_file_name
os.system(command)

