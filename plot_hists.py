
import numpy as np
import sys
import os
import matplotlib.pyplot as plt


datalist = np.loadtxt("us_hists.dat")

for i in range(1,datalist.shape[1]):
	plt.plot(datalist[:,0], datalist[:,i],lw=2)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel("distance", size=12)
plt.ylabel("probability density", size=12)
plt.savefig('hists.png')
plt.close()

