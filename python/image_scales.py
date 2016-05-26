# See what the actual distribution of points in the images are for ATLAS, FIRST images

import pymongo 
from pymongo import MongoClient

client = MongoClient('localhost', 27017)
db = client['radio'] 

subjects = db['radio_subjects'] # subjects = images
classifications = db['radio_classifications']# classifications = classifications of each subject per user

import pandas as pd

data = pd.read_csv('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/
