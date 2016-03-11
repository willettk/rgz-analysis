import rgz

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/plots'

subjects,classifications = rgz.load_rgz_data()

atlas_counts = []

for s in subjects.find({'metadata.survey':'atlas'}):
    atlas_counts.append(s['classification_count'])

from matplotlib import pyplot as plt

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

n,bins,patches = ax.hist(atlas_counts,bins=20,range=(0,20))

ax.set_title('ATLAS subjects - 24 Nov 2015')
ax.set_xlabel('Classification counts')

fig.savefig('%s/atlas_classification_counts.png' % rgz_dir)

# Over the last seven days, what is our average classification rate per day?

import datetime
import numpy as np

lastdate = datetime.datetime(2015,11,16)
startdate = lastdate - datetime.timedelta(7)

fcarr = []
acarr = []

for day in range(7):
    fc,ac = 0,0
    enddate = startdate - datetime.timedelta(1)
    c = classifications.find({"created_at":{"$lte":startdate,"$gt":enddate}})
    for cl in c:
        sid = cl['subject_ids'][0]
        s = subjects.find_one({'_id':sid})
        if s['metadata']['survey'] == 'atlas':
            ac += 1
        else:
            fc += 1
    #print enddate,startdate,fc,ac
    fcarr.append(fc)
    acarr.append(ac)
    startdate = startdate + datetime.timedelta(1)

fcarr = np.array(fcarr)
acarr = np.array(acarr)
print '\nAverage no. of total classifications/day: {0:.0f} +- {1:.0f}'.format(np.mean(fcarr + acarr),np.std(fcarr + acarr))
print '\tFIRST: {0:.0f}'.format(np.mean(fcarr))
print '\tATLAS: {0:.0f}'.format(np.mean(acarr))

# How many until we're done?

remaining_atlas = [x*(20-b) for x,b in zip(n[6:20],bins[6:20])]
print '\n{0:.0f} ATLAS classifications still required.'.format(np.sum(remaining_atlas))
print "At this rate, RGZ will finish the ATLAS sample in {0:.1f} days.\n".format(np.sum(remaining_atlas) * 1./np.mean(acarr))
