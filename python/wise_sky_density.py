# Compute the WISE density of sources in a fixed lat-lon square

from __future__ import division
from matplotlib import pyplot as plt
import numpy as np

nr = 2453

d = np.arange(10)+1
matches_all = np.array([8,27,60,97,143,189,264,331,413,497])
matches_w1w2w3 = np.array([3,10,20,30,49,70,109,132,165,195])       # Come from TOPCAT matches in increasing cone radii of 1 to 10 arcsec

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(d,matches_all/2453,color='b',label='All WISE')
ax.scatter(d,matches_w1w2w3/2453,color='r',label='W1,W2,W3 S/N > 2')
ax.set_xlabel('Matching distance [arcsec]',fontsize=16)
ax.set_ylabel('Likelihood of random match',fontsize=16)
ax.legend(loc='upper left')


plt.show()
