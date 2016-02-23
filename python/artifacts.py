from pymongo import MongoClient 
import json

path = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'

project_name = 'radio'

client = MongoClient('localhost', 27017)
db = client['ouroboros'] 

discussions = db['discussions']
projects = db['projects']

project = projects.find_one({'name':project_name})
artifact_discussions = discussions.find({'project_id':project['_id'], "focus.type":"RadioSubject", "comments.tags":{"$in":["artifact","artifacts","artefact","artefacts"]}})

with open('%s/csv/artifacts_tagged.csv' % path,'wb') as f:
    for ad in artifact_discussions:
        f.write('%s\n' % ad['focus']['_id'])

