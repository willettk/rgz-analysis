# Quick stats on the project

from pymongo import MongoClient
import numpy as np

client = MongoClient('localhost',27017)
radio = client['radio']
subjects = radio['radio_subjects']
classifications = radio['radio_classifications']

def gini(list_of_values):
    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(list_of_values) / 2.
    return (fair_area - area) / fair_area

print "State"
aggstate = subjects.aggregate([{"$group":{"_id":"$state","count":{"$sum":1}}}])
if aggstate['ok']:
    for k in aggstate['result']:
        print "\t id: {0:20} count: {1:6}".format(k['_id'],k['count'])

print "Survey"
aggstate = subjects.aggregate([{"$group":{"_id":"$metadata.survey","count":{"$sum":1}}}])
if aggstate['ok']:
    for k in aggstate['result']:
        print "\t id: {0:20} count: {1:6}".format(k['_id'],k['count'])

print "Volunteers"
aggstate = classifications.aggregate([{"$group":{"_id":"$user_name","count":{"$sum":1}}}])
if aggstate['ok']:
    aggclass = aggstate['result']
    rc = len(aggstate['result'])

print "Total number of subjects: {0}".format(subjects.count())
print "Total number of classifications: {0}".format(classifications.count())
print "Total number of registered users (<1 classification): {0}".format(rc)

counts = []
for x in aggclass:
    if x['_id'] != None:
        counts.append(x['count'])

print "Gini coefficient: {0:.3f}".format(gini(counts))

