# What percentage of the classifications are contributed by the superusers?

import rgz
import cPickle as pkl
import pandas as pd

subjects,classifications,users = rgz.load_rgz_data()

projid = classifications.find_one()['project_id']

names = []
ccu = []
ccc = []

users = classifications.distinct('user_name')

for idx,u in enumerate(users):

    names.append(u['name'])

    # Classification count in the user document doesn't match the actual count in the classifications. Record both.

    ccc.append(classifications.find({'user_name':u['name']}).count())

    try:
        ccu.append(u['projects'][str(projid)]['classification_count'])
    except KeyError:
        ccu.append(0)

    if not idx % 100:
        print '%i/%i' % (idx,len(users))

d = {'names':names,'cc_usercount':ccu,'cc_classcount':ccc}
df = pd.DataFrame(d)

# This is time-intensive, so do the search once and save the result to a file

with open('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/superusers.pkl','wb') as pkl_file:
    pkl.dump(df,pkl_file)

# Analysis

total_classifications = df[df['cc_classcount'] > 0].sort('cc_classcount',ascending=0)['cc_classcount'].sum()
total_classifiers = len(df[df['cc_classcount'] > 0].sort('cc_classcount',ascending=0)['cc_classcount'])

top_5_percent = total_classifiers * 0.05

classifications_by_top_5 = df[df['cc_classcount'] > 0].sort('cc_classcount',ascending=0)['cc_classcount'][:int(top_5_percent)].sum() 
superuser_contributions = classifications_by_top_5/total_classifications
