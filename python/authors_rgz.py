# DEPRECATED
# As of 7 Feb 2016, RGZ data dumps do not include the users collection.

from pymongo import MongoClient
import json

client = MongoClient('localhost', 27017)
db = client['radio'] 

classifications = db['radio_classifications']
users = db['radio_users']	# volunteers doing each classification (can be anonymous)

mrc = classifications.find().sort([("updated_at", -1)]).limit(1)
most_recent_date = [x for x in mrc][0]['updated_at']
rgz_projectid = '52afdb804d69636532000001'

rgzusers = users.find({'projects.%s.classification_count' % rgz_projectid : {'$gt':0}})
nusers = rgzusers.count()

user_json = {'rgz_authors':[]}

for u in rgzusers:
    user_json['rgz_authors'].append(u['name'].strip())

user_json['rgz_authors'].sort()

print '\n%i Radio Galaxy Zoo users as of %s, from %s to %s\n' % (nusers,most_recent_date,user_json['rgz_authors'][0],user_json['rgz_authors'][-1])

with open('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/rgzauthors.json','w') as f:
    json.dump(user_json,f)
