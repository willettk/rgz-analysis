from pymongo import MongoClient
from collections import Counter
import datetime
import numpy as np
import functools,operator
import pprint

client = MongoClient('localhost', 27017)
db = client['radio'] 

main_release_date = datetime.datetime(2013, 12, 17, 0, 0, 0, 0)
subjects = db['radio_subjects'] 		# subjects = images
classifications = db['radio_classifications']	# classifications = classifications of each subject per user

sub = subjects.find_one({'state':'complete','metadata.contour_count':3,'classification_count':20})
imgid = sub['_id']
clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": main_release_date}}))

cdict = {}
clengths = []
for c in clist:
    clengths.append(len(c['annotations']))

gals_unique = Counter(clengths)

print ''
print gals_unique,'\n'

'''
OK - so this isn't good enough at the start, since some users have these metadata annotations and others don't
First step is to find out how many of these are REAL marks:
  - not date, user_agent, or lang
  - not empty in both IR and radio contours

dict: {3:[checksum1,checksum2],4:[checksum1,checksum2,checksum3]}
iterate over keys run Counter on each of these
find which has the largest value
'''

checksum_list = []

for c in clist:
    # Want most popular combo for each NUMBER of galaxies
    
    sumlist = []    # List of the checksums over all possible combinations
    annlist = c['annotations']
    goodann = 0
    for ann in annlist:

        if ann.has_key('started_at') or ann.has_key('finished_at') or ann.has_key('user_agent') or ann.has_key('lang'):
            continue

        xmaxlist = []
        radio_comps = ann['radio']

        # loop over all the radio components within an galaxy
        if radio_comps != 'No Contours':
            for rc in radio_comps:
                xmaxlist.append(float(radio_comps[rc]['xmax']))
        # or make the value -99 if there are no contours
        else:
            xmaxlist.append(-99)

        product = reduce(operator.mul, xmaxlist, 1)
        #print 'xmaxlist,product',xmaxlist,product

        goodann += 1

        sumlist.append(round(product,3))
        #print 'sumlist',sumlist

    checksum = sum(sumlist)
    checksum_list.append(checksum)
    c['checksum'] = checksum

    if cdict.has_key(goodann):
        cdict[goodann].append(checksum)
    else:
        cdict[goodann] = [checksum]

print cdict,'\n'

maxval=0
mc_checksum = 0.
ngals = 0
for k,v in cdict.iteritems():
    mc = Counter(v).most_common()
    print '%i galaxies: %s' % (k,mc)
    if mc[0][0] == -99.0:
        if len(mc) > 1:
            mc_best = mc[1]
        else:
            continue
    else:
        mc_best = mc[0]
    if mc_best[1] > maxval:
        maxval = mc_best[1]
        mc_checksum = mc_best[0]
        ngals = k

suffix = 'y' if ngals == 1 else 'ies'
print 'Most common selection (%i/%i users) is %i galax%s\n' % (maxval,len(clist),ngals,suffix)

# Find a galaxy that matches the checksum (easier to keep track as a list)

cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
pprint.pprint(cmatch['annotations'])
print ''
print 'http://radiotalk.galaxyzoo.org/#/subjects/%s\n' % sub['zooniverse_id']

print imgid
