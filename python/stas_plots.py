import rgz
from matplotlib import pyplot as plt

subjects,classifications = rgz.load_rgz_data()
catalog = rgz.load_catalog()

cc_list = list(subjects.find({'metadata.survey':'first'}))

cc_all = [x['metadata']['contour_count'] for x in cc_list]
cc_complete = [x['metadata']['contour_count'] if x['state'] == 'complete' else -999 for x in cc_list]
all_list = list(catalog.find({'consensus.level':{'$gte':0.50}}))
cc_consensus = [x['radio']['numberComponents'] for x in all_list]

ncomp = 20

def set_legend_props(legend):

    for label in legend.get_texts():
        label.set_fontsize('x-small')
    for label in legend.get_lines():
        label.set_linewidth(1.5)

    return None

fig = plt.figure(1,figsize=(8,8))

# Plot: histogram of number of components per source overall

ax1 = fig.add_subplot(221)
ax1.hist(cc_all,bins=ncomp,range=(0,ncomp),label='All subjects')
ax1.hist(cc_complete,bins=ncomp,range=(0,ncomp),label='Completed subjects')
ax1.hist(cc_consensus,bins=ncomp,range=(0,ncomp),label=r'$\geq$'+'50% consensus')
ax1.set_yscale('log')
ax1.set_xlabel('Number of radio components in image',fontsize=10)
ax1.set_ylabel('Count',fontsize=10)
ax1.set_title('All FIRST subjects',fontsize=11)

legend1 = ax1.legend(loc='upper right', shadow=True)
set_legend_props(legend1)

# Plot: histogram of number of components per source with an IR counterpart

ir_list = list(catalog.find({'AllWISE':{'$exists':True},'consensus.level':{'$gte':0.50}}))
ir_all = [x['radio']['numberComponents'] for x in ir_list]

optical_list = list(catalog.find({'SDSS':{'$exists':True},'consensus.level':{'$gte':0.50}}))
optical_all = [x['radio']['numberComponents'] for x in optical_list]

ax2 = fig.add_subplot(222)
ax2.hist(ir_all,bins=ncomp,range=(0,ncomp),label='AllWISE IR host')
ax2.hist(optical_all,bins=ncomp,range=(0,ncomp),label='SDSS optical host')
ax2.set_yscale('log')
ax2.set_xlabel('Number of radio components in image',fontsize=10)
ax2.set_ylabel('Count',fontsize=10)
ax2.set_title('FIRST subjects with an IDed host',fontsize=11)

legend2 = ax2.legend(loc='upper right', shadow=True)
set_legend_props(legend2)

# Plot: number of RGZ sources per image

temp = [catalog.find({'consensus.label':letter}).count() for letter in 'abcdefghijklmnopqrstuvwxyz']
sourcecount = []
for idx,t in enumerate(temp):
    sourcecount.extend([idx+1]*t)

ax4 = fig.add_subplot(224)
ax4.hist(sourcecount,bins=11,range=(0,11))
ax4.set_yscale('log')
ax4.set_xlabel('Number of discrete sources in RGZ image',fontsize=10)
ax4.set_ylabel('Count',fontsize=10)

fig.tight_layout()

#plt.show()

fig.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/plots/stas_rgz_plots.pdf')

