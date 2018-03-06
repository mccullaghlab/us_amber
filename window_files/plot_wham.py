
# load libraries
import sys
import os
import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

kT = 0.593

data = np.loadtxt("wham.out")

data[:,1] += 2*kT*np.log(data[:,0])
data[:,1] -= np.mean(data[-10:-1,1])
#data[:,2] += 2*kT*np.log(data[:,0])
plt.errorbar(data[:,0],data[:,1],yerr=data[:,2],lw=2)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel("distance", size=12)
plt.ylabel("free energy", size=12)
plt.savefig('freeEnergy.png')
plt.close()

