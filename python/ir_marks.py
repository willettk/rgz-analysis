from pymongo import MongoClient
from bson import ObjectId
from collections import Counter


client = MongoClient('localhost', 27017)
db = client['radio'] 
subjects = db['radio_subjects']
classifications = db['radio_classifications']

'''
atlas,first,other = 0,0,0
wf = open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/ir0_coo.csv','w')
print >> wf,"ra,dec,n"
with open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/ir0_subjects.csv','r') as f:
    for line in f.readlines()[1:]:
        objectid = ObjectId(line.strip()[9:-1])
        sub = subjects.find_one({"_id":objectid})
        print >> wf,"{0:.6f},{1:.6f},{2:d}".format(sub['coords'][0],sub['coords'][1],sub['metadata']['contour_count'])
        if sub['metadata']['survey'] == 'atlas':
            atlas += 1
        elif sub['metadata']['survey'] == 'first':
            first += 1
        else:
            other += 1
'''

#print "ATLAS: {0}".format(atlas)
#print "FIRST: {0}".format(first)
#print "other: {0}".format(other)

'''
user_agents = []
with open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/ir0_exterior.csv','r') as f:
    for line in f.readlines()[1:]:
        objectid = ObjectId(line.strip().split(',')[-1][9:-1])
        c = classifications.find_one({"_id":objectid})
        try:
            user_agents.append(c['annotations'][-1]['user_agent'])
        except:
            print "No user agent found for {0}".format(objectid)

c = Counter(user_agents)
for k in c:
    print c[k],'\t',k

'''

from subprocess import call
ir_file = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/ir0_all.csv'
with open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/goldstandard/gs_zids.txt','r') as gs:
    for line in gs.readlines():
        sub = subjects.find_one({'zooniverse_id':line.strip()})
        call(["grep"," -v {0} {1} > temp_rm".format(sub['_id'],ir_file)])
        call(["mv","temp_rm {0}".format(ir_file)])
