# Create the expert sample for Radio Galaxy Zoo. #
# Start out as 100 or so galaxies, to be classified by Larry Rudnick and Kyle Willett. 
# We'll see how this compares to the user classifications done so far.

import rgz
import webbrowser
import itertools
import numpy as np
from numpy.random import random

curated_zid = ( 'ARG0003j4t', 'ARG0002wze', 'ARG0000slz', 'ARG0001nbb', 'ARG0000lv1', 'ARG0002rwh', 'ARG0001gxs', 'ARG0002rx3', 'ARG0001e8e', 'ARG0000nbh', 'ARG0001d12', 'ARG0002kve', 'ARG0000v25', 'ARG0000jqj', 'ARG0001r2h', 'ARG0002r8g', 'ARG00011aj', 'ARG0001fsn', 'ARG0000y8j', 'ARG0003bq8', 'ARG0001kqx', 'ARG0001qy6', 'ARG0002ato', 'ARG000180p', 'ARG0002xat', 'ARG0001rcw', 'ARG0000dju', 'ARG00029d6', 'ARG0002pf3', 'ARG0001u0z', 'ARG0003bul', 'ARG000090l', 'ARG0000czo', 'ARG000051b', 'ARG0000p2c', 'ARG0003dzd', 'ARG0002yeu', 'ARG00022wv', 'ARG0002wel', 'ARG00025a1', 'ARG0001n7u', 'ARG0003j0a', 'ARG00004rl', 'ARG0000t1w', 'ARG0000pfu', 'ARG0002ehu', 'ARG0002mfx', 'ARG0000zeh', 'ARG0002c2p', 'ARG0002qfa', 'ARG0002jz8', 'ARG0001p9x', 'ARG0003aob', 'ARG0001fsk', 'ARG0000yhw', 'ARG0002yms', 'ARG0000vcr')

# Select equal number of random galaxies

# Note - Chris Snyder recommends chunking it into 25 or so at once. Not sure if the browser will like preloading and caching 100 subjects all at once.
# Formal max limit from the API is 100 subjects

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
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

rgz_chain = itertools.chain(*zip(curated_zid,random_zid))
result = list(rgz_chain)

with open('%s/expert_urls.txt' % rgz_dir,'wb') as f:
    for i in np.arange(4):
        url_str = 'http://radio.galaxyzoo.org/?subjects='
        for r in result[i::4]:
            url_str += '%s,' % r
            
        f.write(url_str[:-1]+'\n')


