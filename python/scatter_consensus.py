import numpy as np
from matplotlib import pyplot as plt

# Data
# mostly from Google Drive spreadsheet: https://docs.google.com/spreadsheets/d/1yydM_FmGU21uMzC5BNNAd4IDS-GepSLoskmloxPF4gI/edit#gid=0
# and output of consensus.py

cv = np.array([0.31, 0.90, 0.36, 1.00, 0.71, 0.92, 1.00, 0.38, 1.00, 0.89, 0.89, 0.67, 0.80, 0.33, 0.82, 0.50, 0.66, 0.89, 0.33, 0.69, 1.00, 0.21, 0.74, 1.00, 0.64, 0.59, 0.33, 0.31, 0.39, 0.82, 0.83, 0.67, 1.00, 0.93, 0.56, 0.59, 0.64, 0.66, 0.61, 0.86, 0.71, 0.43, 1.00, 0.31, 0.63, 0.89, 0.91, 0.20, 0.71, 0.53, 0.61, 0.58, 0.57, 0.50, 0.64, 0.85, 0.88, 0.33, 1.00, 0.47, 0.92, 0.41, 0.56, 0.90, 0.44, 0.81, 1.00, 0.54, 0.93, 0.67, 1.00, 0.36, 0.70, 0.47, 0.58, 0.20, 1.00, 0.38, 0.38, 0.55, 0.78, 0.50, 0.31, 0.86, 0.83, 0.56, 0.57, 0.77, 0.87, 0.64, 0.96, 0.33, 0.36, 0.50, 0.88, 0.89, 1.00, 0.44, 0.90, 0.83])

ce = np.array([0.5, 1.0, 0.625, 1.0, 0.888888888889, 1.0, 1.0, 0.875, 1.0, 1.0, 1.0, 0.75, 1.0, 0.75, 1.0, 0.555555555556, 1.0, 1.0, 0.375, 1.0, 1.0, 0.375, 1.0, 1.0, 0.625, 0.75, 1.0, 0.75, 0.875, 1.0, 1.0, 0.777777777778, 1.0, 1.0, 0.875, 1.0, 0.875, 0.875, 1.0, 0.875, 0.888888888889, 0.5, 1.0, 0.75, 0.888888888889, 1.0, 1.0, 0.375, 0.75, 0.75, 1.0, 0.875, 0.888888888889, 0.5, 0.888888888889, 1.0, 1.0, 0.75, 1.0, 1.0, 1.0, 0.75, 1.0, 1.0, 0.75, 1.0, 0.75, 0.625, 1.0, 0.888888888889, 1.0, 0.5, 0.777777777778, 0.75, 0.666666666667, 0.5, 1.0, 0.5, 0.625, 0.666666666667, 0.75, 0.888888888889, 0.875, 1.0, 1.0, 0.875, 0.888888888889, 0.875, 1.0, 0.625, 1.0, 0.75, 0.875, 0.555555555556, 0.875, 1.0, 1.0, 0.875, 0.888888888889, 1.0])

ns_expert = np.array([2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1])

expert_class = np.array(['c', 'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'c', 'a', 'a', 'b', 'a', 'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'b', 'a', 'a', 'a', 'a', 'a', 'c', 'b', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'c', 'a', 'a', 'a', 'b', 'a', 'a', 'b', 'a', 'a', 'c', 'a', 'a', 'a', 'c', 'b', 'b', 'c', 'c', 'a', 'c', 'a', 'c', 'b', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'c', 'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a'])

vol_class = np.array(['c', 'a', 'c', 'a', 'a', 'a', 'a', 'c', 'a', 'a', 'a', 'b', 'b', 'c', 'b', 'c', 'b', 'a', 'c', 'b', 'a', 'c', 'b', 'a', 'c', 'b', 'a', 'c', 'b', 'a', 'b', 'b', 'a', 'a', 'c', 'b', 'a', 'b', 'b', 'a', 'b', 'b', 'a', 'b', 'b', 'a', 'a', 'c', 'b', 'b', 'b', 'b', 'c', 'b', 'b', 'a', 'a', 'c', 'a', 'a', 'a', 'b', 'b', 'a', 'b', 'a', 'a', 'c', 'a', 'b', 'a', 'c', 'a', 'a', 'c', 'c', 'a', 'c', 'b', 'c', 'b', 'c', 'a', 'a', 'a', 'b', 'b', 'b', 'a', 'a', 'a', 'c', 'c', 'c', 'a', 'a', 'a', 'c', 'a', 'a'])

ecn = np.array([3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 2, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 1, 3, 2, 2, 3, 3, 1, 3, 1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 1, 1, 1, 1])

vcn = np.array([3, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 2, 2, 3, 2, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 3, 2, 1, 2, 2, 1, 1, 3, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 3, 2, 2, 2, 2, 3, 2, 2, 1, 1, 3, 1, 1, 1, 2, 2, 1, 2, 1, 1, 3, 1, 2, 1, 3, 1, 1, 3, 3, 1, 3, 2, 3, 2, 3, 1, 1, 1, 2, 2, 2, 1, 1, 1, 3, 3, 3, 1, 1, 1, 3, 1, 1])

ve = np.array([False, True, True, True, True, True, True, False, True, True, False, False, False, True, True, True, True, True, True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, True, False, True, True, True, True, True, True, False, True, True, True, True, True, False, True, True, True, True, False, False, True, True, True, False, True, True, True, True, True, True, False, True, True, False, True, True, True, False, False, False, False, False, True, True, True, True, True, True, False, True, True, True, True, True, False, False, True, True, True, False, True, True, True, True, True, True])

# Plot the agreement in consensus for the 100 galaxy control sample, experts vs. volunteers

fig = plt.figure(1,(8,8))
ax1 = fig.add_subplot(111)

a_agree = ve & (expert_class == 'a')
b_agree = ve & (expert_class == 'b')
c_agree = ve & (expert_class == 'c')
a_disagree = np.logical_not(ve) & (expert_class == 'a')
b_disagree = np.logical_not(ve) & (expert_class == 'b')
c_disagree = np.logical_not(ve) & (expert_class == 'c')

sa1 = ax1.scatter(ce[a_disagree],cv[a_disagree],s=80,facecolors='none',edgecolors='r')
sa2 = ax1.scatter(ce[b_disagree],cv[b_disagree],s=80,facecolors='none',edgecolors='b')
sa3 = ax1.scatter(ce[c_disagree],cv[c_disagree],s=80,facecolors='none',edgecolors='y')
sd1 = ax1.scatter(ce[a_agree],cv[a_agree],s=80,facecolors='r')
sd2 = ax1.scatter(ce[b_agree],cv[b_agree],s=80,facecolors='b')
sd3 = ax1.scatter(ce[c_agree],cv[c_agree],s=80,facecolors='y')

ax1.set_xlim(0,1.05)
ax1.set_ylim(0,1.05)
ax1.set_xlabel('Expert RGZ consensus',fontsize=16)
ax1.set_ylabel('Volunteer RGZ consensus',fontsize=16)

ax1.legend((sa1,sd1),('Agree','Disagree'),loc='upper left',fontsize=14,scatterpoints=1)

fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/figures/scatter_consensus.eps')

