# -*- coding: utf-8 -*-
"""
Vojtěch Kulvait
"""

import numpy as np
import matplotlib.pyplot as plt
dat = np.genfromtxt('convergence',dtype=None, skip_header=1)
normb = dat[0,0]
dat = dat*100/normb
#dat = np.delete(dat, 0,axis=0)

plt.rcParams.update({'font.size': 12})
#plt.figure(figsize=(6, 3.9), dpi=72)
plt.figure(figsize=(6,3.9))
plt.tight_layout()
plt.plot(dat[:,0], label="CGLS")
plt.plot(dat[:,1], label="PSIRT")
plt.xlabel("Number of iterations")
plt.ylabel("Relative discrepancy [%]")
plt.legend()
plt.yscale("log")
#plt.ylim([-1,25])
plt.show()
plt.savefig('convergenceLOG.png', dpi=300)

plt.figure(figsize=(6,3.9))
plt.tight_layout()
plt.plot(dat[:,0], label="CGLS")
plt.plot(dat[:,1], label="PSIRT")
plt.xlabel("Number of iterations")
plt.ylabel("Relative discrepancy [%]")
plt.legend(loc=1)
plt.xlim([-1,40])
plt.ylim([-1,25])
#plt.figure(figsize=(3.5, 2.2), dpi=150)
plt.show()
plt.savefig('convergenceA.png', dpi=300)
