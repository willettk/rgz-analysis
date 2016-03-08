# Create the expert sample for Radio Galaxy Zoo.
# We'll see how this compares to the user classifications done so far.

# Kyle Willett, 10 Jun 2014

import rgz
from numpy.random import random
from random import shuffle

# Paths
rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

# Curated subjects from http://radiogalaxyzoo.pbworks.com/w/page/81017963/Expert%20classifications
curated_zid = open('%s/expert/expert_curated_zooniverse_ids.txt' % rgz_dir).read().splitlines()
    
# Note - Chris Snyder recommends chunking it into 25 or so at once. 
# Not sure if the browser will like preloading and caching 100 subjects all at once.
# Formal max limit from the API is 100 subjects

# Select random galaxies to fill out sample; total will be 100 subjects

n_random = 100 - len(curated_zid)

# Load in the RGZ database through pymongo
subjects,classifications = rgz.load_rgz_data()

# Completed RGZ subjects
batch = subjects.find({'state':'complete','classification_count':20}).limit(n_random)

# Select random subjects from the completed database

random_zid = []
random_numbers = random(n_random)

for r in random_numbers:
    sf = subjects.find({'random':{'$gte':random()}}).sort([("random", 1)]).limit(1)
    s = list(sf)[0]
    random_zid.append(s['zooniverse_id'])

# Combine the curated and randomly selected subjects
result = curated_zid + random_zid

# Randomize the order that they appear in
shuffle(result)

# Write out the URLs to a text file, batched in units of 25
# Note: these are later put through a URL shortener for ease in copying/pasting

with open('%s/expert/expert_urls.txt' % rgz_dir,'wb') as f:
    for i in range(4):
        url_str = 'http://radio.galaxyzoo.org/?subjects='
        for r in result[i::4]:
            url_str += '%s,' % r
            
        f.write(url_str[:-1]+'\n')

# Write the Zooniverse IDs to a text file as a sorted list

with open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir,'wb') as f:
    sorted_result = result.sort()
    for r in result:
        f.write(r+'\n')

