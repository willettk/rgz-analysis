# How many sources per field?

import json
from matplotlib import pyplot as plt

with open("/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/json/consensus_rgz_first_75.json",'r') as f:
    first = json.load(f)

# Eliminate the empty fields (should be done in the catalog step)
fsources = []
for x in first:
    if len(x['answer']) > 0:
        fsources.append(len(x['answer']))


with open("/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/json/consensus_rgz_atlas_75.json",'r') as f:
    atlas = json.load(f)

asources = [len(x['answer']) for x in atlas]

# Changes: plot CDF, not histogram. Also scale by the size of the field (and possibly sensitivity?)

fig,axarr = plt.subplots(nrows=1,ncols=2,figsize=(10,4))

ax1,ax2 = axarr.ravel()

ax1.hist(fsources,bins=range(10),histtype='stepfilled',label='FIRST')
ax1.hist(asources,bins=range(10),histtype='step',label='ATLAS',lw=2)
ax1.set_xlabel("Sources per image w/75% consensus")
ax1.set_yscale('log')
ax1.legend(loc='best')

ax2.hist(fsources,bins=range(10),normed=True,histtype='stepfilled',label='FIRST')
ax2.hist(asources,bins=range(10),normed=True,histtype='step',label='ATLAS',lw=2)
ax2.set_xlabel("Sources per image w/75% consensus")
ax2.legend(loc='best')

plt.show()
