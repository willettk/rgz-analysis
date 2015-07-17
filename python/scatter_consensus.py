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

ir_position_agrees = np.ones(len(cv),dtype=bool)
for i in (10,11,12,23,25,52,53,72,73,82,88,89):
    ir_position_agrees[i] = False

names = ['ARG00000b0', 'ARG00003y9', 'ARG00004rl', 'ARG00004um', 'ARG000051b', 'ARG000051q', 'ARG00008v5', 'ARG000090l', 'ARG0000cps', 'ARG0000czo', 'ARG0000d56', 'ARG0000dju', 'ARG0000fle', 'ARG0000jqj', 'ARG0000lv1', 'ARG0000mpd', 'ARG0000nbh', 'ARG0000opo', 'ARG0000p2c', 'ARG0000pfu', 'ARG0000slz', 'ARG0000t1w', 'ARG0000v25', 'ARG0000vcr', 'ARG0000y8j', 'ARG0000yhw', 'ARG0000yw0', 'ARG0000zaz', 'ARG0000zeh', 'ARG00011aj', 'ARG000180p', 'ARG00018qp', 'ARG0001bnj', 'ARG0001c6k', 'ARG0001d12', 'ARG0001e8e', 'ARG0001fsk', 'ARG0001fsn', 'ARG0001gxs', 'ARG0001gxz', 'ARG0001kqx', 'ARG0001lau', 'ARG0001mzk', 'ARG0001n7u', 'ARG0001nbb', 'ARG0001o4l', 'ARG0001p9x', 'ARG0001qy6', 'ARG0001r2h', 'ARG0001rcw', 'ARG0001u0z', 'ARG0001ugd', 'ARG0001wnd', 'ARG0001ynm', 'ARG00020o8', 'ARG00021gd', 'ARG00022wv', 'ARG00025a1', 'ARG000269p', 'ARG00029d6', 'ARG0002ato', 'ARG0002c2p', 'ARG0002ehu', 'ARG0002gz5', 'ARG0002jz8', 'ARG0002k8k', 'ARG0002kcp', 'ARG0002kve', 'ARG0002mfx', 'ARG0002nl1', 'ARG0002pf3', 'ARG0002qfa', 'ARG0002qj7', 'ARG0002r8g', 'ARG0002rwh', 'ARG0002rx3', 'ARG0002w9n', 'ARG0002wel', 'ARG0002wze', 'ARG0002xat', 'ARG0002y2z', 'ARG0002yeu', 'ARG0002yms', 'ARG0002yzd', 'ARG000316y', 'ARG00031s5', 'ARG00031s6', 'ARG00037z5', 'ARG00039w1', 'ARG0003aob', 'ARG0003bq8', 'ARG0003bul', 'ARG0003ccm', 'ARG0003dzd', 'ARG0003glu', 'ARG0003ioj', 'ARG0003j0a', 'ARG0003j4t', 'ARG0003jlu', 'ARG0003knd']


# Axes definitions
nullfmt = plt.NullFormatter()
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# Generate initial figure, scatter plot, and histogram quadrants
fig = plt.figure(221, figsize=(10,10))

axScatter = fig.add_subplot(223, position=rect_scatter)

red = '#e41a1c'
blue = '#377eb8'
green = '#4daf4a'

axHistX = fig.add_subplot(221, position=rect_histx)
axHistX.set_xlim(0,1.05)
axHistX.set_ylim(0, 50.)

axHistY = fig.add_subplot(224, position=rect_histy)
axHistY.set_xlim(0, 25.)
axHistY.set_ylim(0,1.05)

# Remove labels from histogram edges touching scatter plot
axHistX.xaxis.set_major_formatter(nullfmt)
axHistY.yaxis.set_major_formatter(nullfmt)

# Draw scatter plot
#axScatter.scatter(x, y, marker='o', color = 'darkblue', edgecolor='none', s=5, alpha=1)

#fig = plt.figure(1,(8,8))
#ax1 = fig.add_subplot(111)

'''
This includes manual classifications, after experts discussed mistakes or techniques and the A,B,C classifications were adjusted

It also includes things like tiny components that are clearly not part of the overall structure, and that would be removed with a size/flux cut

If these are to be used, then the vote fractions for experts should also be changed.


a_agree = ve & (expert_class == 'a')
b_agree = ve & (expert_class == 'b')
c_agree = ve & (expert_class == 'c')
a_disagree = np.logical_not(ve) & (expert_class == 'a')
b_disagree = np.logical_not(ve) & (expert_class == 'b')
c_disagree = np.logical_not(ve) & (expert_class == 'c')
'''

# Cutoffs assume 9(10) classifiers and A = 8-9 agree, B = 6-7 agree, C = 0-5 agree

a_agree = ve & (ce >= 0.88)
b_agree = ve & (ce < 0.88) & (ce >= 0.66)
c_agree = ve & (ce < 0.66)
a_disagree = np.logical_not(ve) & (ce >= 0.88)
b_disagree = np.logical_not(ve) & (ce < 0.88) & (ce >= 0.66)
c_disagree = np.logical_not(ve) & (ce < 0.66)

sa_blank = axScatter.scatter(ce[a_agree][0],cv[a_agree][0],s=80,facecolors='k')
sd_blank = axScatter.scatter(ce[a_disagree][0],cv[a_disagree][0],s=80,facecolors='none',edgecolors='k')
sa1 = axScatter.scatter(ce[a_disagree],cv[a_disagree],s=80,facecolors='none',edgecolors=red)
sa2 = axScatter.scatter(ce[b_disagree],cv[b_disagree],s=80,facecolors='none',edgecolors=blue)
sa3 = axScatter.scatter(ce[c_disagree],cv[c_disagree],s=80,facecolors='none',edgecolors=green)
sd1 = axScatter.scatter(ce[a_agree],cv[a_agree],s=80,facecolors=red)
sd2 = axScatter.scatter(ce[b_agree],cv[b_agree],s=80,facecolors=blue)
sd3 = axScatter.scatter(ce[c_agree],cv[c_agree],s=80,facecolors=green)

axScatter.set_xlim(0,1.05)
axScatter.set_ylim(0,1.05)
axScatter.set_xlabel('Expert RGZ consensus',fontsize=20)
axScatter.set_ylabel('Volunteer RGZ consensus',fontsize=20)

axScatter.legend((sa_blank,sd_blank),('Agree','Disagree'),loc='upper left',fontsize=17,scatterpoints=1)

# Calculate number of bins based on binsize for both x and y
x1 = ce[ce >= 0.88]
x2 = ce[(ce < 0.88) & (ce >= 0.66)]
x3 = ce[ce < 0.66]

y1 = cv[ce >= 0.88]
y2 = cv[(ce < 0.88) & (ce >= 0.66)]
y3 = cv[ce < 0.66]

num_bins = 10

'''
axHistX.hist(x1, num_bins, range=(0,1), ec='r', fc='r', histtype='stepfilled',lw=2,alpha=0.5)
axHistX.hist(x2, num_bins, range=(0,1), ec='b', fc='none', histtype='stepfilled',lw=4,alpha=1.0)
axHistX.hist(x3, num_bins, range=(0,1), ec='g', fc='g', histtype='stepfilled',lw=2,alpha=1.0)
axHistY.hist(y1, num_bins, range=(0,1), ec='r', fc='r', histtype='stepfilled',lw=2,alpha=0.5,orientation='horizontal')
axHistY.hist(y2, num_bins, range=(0,1), ec='b', fc='none', histtype='stepfilled',lw=4,alpha=1.0,orientation='horizontal')
axHistY.hist(y3, num_bins, range=(0,1), ec='g', fc='g', histtype='stepfilled',lw=2,alpha=1.0,orientation='horizontal')
'''

axHistX.hist((x1,x2,x3), num_bins, range=(0,1), color=(red,blue,green), stacked=True, histtype='stepfilled',lw=2,alpha=1.0)
axHistY.hist((y1,y2,y3), num_bins, range=(0,1), color=(red,blue,green), stacked=True, histtype='stepfilled',lw=2,alpha=1.0,orientation='horizontal')

axScatter.tick_params(axis='x', labelsize=16)
axScatter.tick_params(axis='y', labelsize=16)

'''
axHistY.set_xlabel("Count",fontsize=10)
axHistX.set_ylabel("Count",fontsize=10)
'''

#plt.show()

fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/figures/scatter_hist_consensus.eps')
fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/paper/fig8.eps')

# Print results

ae_av = (ce >= 0.88) & (cv >= 0.88)
ae_bv = (ce >= 0.88) & (cv < 0.88) & (cv >= 0.66)
ae_cv = (ce >= 0.88) & (cv < 0.66)
be_av = (ce < 0.88) & (ce >= 0.66) & (cv >= 0.88)
be_bv = (ce < 0.88) & (ce >= 0.66) & (cv < 0.88) & (cv >= 0.66)
be_cv = (ce < 0.88) & (ce >= 0.66) & (cv < 0.66)
ce_av = (ce < 0.66) & (cv >= 0.88)
ce_bv = (ce < 0.66) & (cv < 0.88) & (cv >= 0.66)
ce_cv = (ce < 0.66) & (cv < 0.66)

ae_av_agree = ve & (ce >= 0.88) & (cv >= 0.88)
ae_bv_agree = ve & (ce >= 0.88) & (cv < 0.88) & (cv >= 0.66)
ae_cv_agree = ve & (ce >= 0.88) & (cv < 0.66)
be_av_agree = ve & (ce < 0.88) & (ce >= 0.66) & (cv >= 0.88)
be_bv_agree = ve & (ce < 0.88) & (ce >= 0.66) & (cv < 0.88) & (cv >= 0.66)
be_cv_agree = ve & (ce < 0.88) & (ce >= 0.66) & (cv < 0.66)
ce_av_agree = ve & (ce < 0.66) & (cv >= 0.88)
ce_bv_agree = ve & (ce < 0.66) & (cv < 0.88) & (cv >= 0.66)
ce_cv_agree = ve & (ce < 0.66) & (cv < 0.66)

ae_av_disagree = np.logical_not(ve) & (ce >= 0.88) & (cv >= 0.88)
ae_bv_disagree = np.logical_not(ve) & (ce >= 0.88) & (cv < 0.88) & (cv >= 0.66)
ae_cv_disagree = np.logical_not(ve) & (ce >= 0.88) & (cv < 0.66)
be_av_disagree = np.logical_not(ve) & (ce < 0.88) & (ce >= 0.66) & (cv >= 0.88)
be_bv_disagree = np.logical_not(ve) & (ce < 0.88) & (ce >= 0.66) & (cv < 0.88) & (cv >= 0.66)
be_cv_disagree = np.logical_not(ve) & (ce < 0.88) & (ce >= 0.66) & (cv < 0.66)
ce_av_disagree = np.logical_not(ve) & (ce < 0.66) & (cv >= 0.88)
ce_bv_disagree = np.logical_not(ve) & (ce < 0.66) & (cv < 0.88) & (cv >= 0.66)
ce_cv_disagree = np.logical_not(ve) & (ce < 0.66) & (cv < 0.66)

print 'A & &      %2i      &   %2i      &      %2i \\\\'    % (ae_av.sum(),ae_bv.sum(),ae_cv.sum())
print 'B & &      %2i      &   %2i      &      %2i \\\\'    % (be_av.sum(),be_bv.sum(),be_cv.sum())
print 'C & &      %2i      &   %2i      &      %2i \\\\ \n' % (ce_av.sum(),ce_bv.sum(),ce_cv.sum())

print 'A & &      %2i      &   %2i      &      %2i \\\\'    % (ae_av_agree.sum(),ae_bv_agree.sum(),ae_cv_agree.sum())
print 'B & &      %2i      &   %2i      &      %2i \\\\'    % (be_av_agree.sum(),be_bv_agree.sum(),be_cv_agree.sum())
print 'C & &      %2i      &   %2i      &      %2i \\\\'    % (ce_av_agree.sum(),ce_bv_agree.sum(),ce_cv_agree.sum())
print '\\hline'
print 'A & &      %2i      &   %2i      &      %2i \\\\'    % (ae_av_disagree.sum(),ae_bv_disagree.sum(),ae_cv_disagree.sum())
print 'B & &      %2i      &   %2i      &      %2i \\\\'    % (be_av_disagree.sum(),be_bv_disagree.sum(),be_cv_disagree.sum())
print 'C & &      %2i      &   %2i      &      %2i \\\\ \n' % (ce_av_disagree.sum(),ce_bv_disagree.sum(),ce_cv_disagree.sum())

print 'Total galaxies: %i' % np.sum((ae_av.sum(),ae_bv.sum(),ae_cv.sum(), be_av.sum(),be_bv.sum(),be_cv.sum(), ce_av.sum(),ce_bv.sum(),ce_cv.sum()))
print 'Number that agreed: %i' % np.sum(ve)
print 'Number that agreed in A or B: %i' % np.sum(ve[ce >= 0.66])

print '\nTotal A galaxies: %i' % np.sum([ae_av.sum(),ae_bv.sum(),ae_cv.sum()])
print 'Total B galaxies: %i' % np.sum([be_av.sum(),be_bv.sum(),be_cv.sum()])
print 'Total C galaxies: %i' % np.sum([ce_av.sum(),ce_bv.sum(),ce_cv.sum()])

ab_no_ir = np.logical_not(ir_position_agrees[ae_av_disagree | ae_bv_disagree | ae_cv_disagree | be_av_disagree | be_bv_disagree | be_cv_disagree]).sum()
total = np.sum([(ae_av_disagree.sum(),ae_bv_disagree.sum(),ae_cv_disagree.sum(),be_av_disagree.sum(),be_bv_disagree.sum(),be_cv_disagree.sum())])
print '\n%i/%i of the expert A/B sample with which volunteers disagreed were due to IR differences, not radio' % (ab_no_ir,total)
