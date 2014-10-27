# What is the distribution of the number of radio components per image? 

import rgz

subjects,classifications,users = rgz.load_rgz_data()

cc_list = list(subjects.find())

cc_all = [x['metadata']['contour_count'] for x in cc_list]
cc_complete = [x['metadata']['contour_count'] if x['state'] == 'complete' else -999 for x in cc_list]

from matplotlib import pyplot as plt

ncomp = 40
fig = plt.figure(1,figsize=(8,6))
ax = fig.add_subplot(111)
ax.hist(cc_all,bins=ncomp,range=(0,ncomp),label='All RGZ subjects')
ax.hist(cc_complete,bins=ncomp,range=(0,ncomp),label='Completed')
ax.set_yscale('log')
ax.set_xlabel('Number of radio components in image',fontsize=16)
ax.set_ylabel('Count',fontsize=16)
legend = ax.legend(loc='upper right', shadow=True)
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('medium')
for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
    plt.show()

plt.show()

fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/plots/hist_radiocomponents.pdf')

# There are 182 subjects with zero radio contours in the image. 
# 180 of them are labeled as "disabled" (will not be classified), 1 is complete and 1 is inactive. 
# Make sure to eliminate these if any show up in final classifications.

# Interestingly, though, these are not images without radio emission; the shading indicates that there is FIRST data there. 
# For some reason, the contours weren't appropriately generated in the offline processing. 
