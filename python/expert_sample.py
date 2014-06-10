# Create the expert sample for Radio Galaxy Zoo. #
# Start out as 100 or so galaxies, to be classified by Larry Rudnick and Kyle Willett. 
# We'll see how this compares to the user classifications done so far.

import rgz
import webbrowser
import itertools
import numpy as np
from numpy.random import random
from random import shuffle

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
curated_zid = open('%s/expert/expert_curated_zooniverse_ids.txt' % rgz_dir).read().splitlines()
    
# Select equal number of random galaxies

# Note - Chris Snyder recommends chunking it into 25 or so at once. Not sure if the browser will like preloading and caching 100 subjects all at once.
# Formal max limit from the API is 100 subjects

n_random = 100 - len(curated_zid)

subjects,classifications,users = rgz.load_rgz_data()

batch = subjects.find({'state':'complete','classification_count':20}).limit(n_random)

random_zid = []
random_numbers = random(n_random)

for r in random_numbers:
    sf = subjects.find({'random':{'$gte':random()}}).sort([("random", 1)]).limit(1)
    s = list(sf)[0]
    random_zid.append(s['zooniverse_id'])
    url = "http://radiotalk.galaxyzoo.org/#/subjects/%s" % s['zooniverse_id']
    #print "%8.6f, %s" % (r,url)
    #webbrowser.open(url)

# Combine the list, randomize their order

'''
rgz_chain = itertools.chain(*zip(curated_zid,random_zid))
result = list(rgz_chain)
'''

result = curated_zid + random_zid
shuffle(result)

with open('%s/expert_urls_test.txt' % rgz_dir,'wb') as f:
    for i in np.arange(4):
        url_str = 'http://radio.galaxyzoo.org/?subjects='
        for r in result[i::4]:
            url_str += '%s,' % r
            
        f.write(url_str[:-1]+'\n')


