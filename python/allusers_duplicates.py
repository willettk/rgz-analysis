import rgz
import numpy as np
import pandas as pd
import datetime

subjects,classifications,users = rgz.load_rgz_data()

javascript_fix_date = datetime.datetime(2014, 3, 24, 0, 0, 0, 0)

# Search for duplicates among all users

recent_date = datetime.datetime(2014, 6, 6, 0, 0, 0, 0)
user_list = []
for c in classifications.find({'user_name':{'$exists':'true'},'created_at':{'$gt':recent_date}},{'user_name':1}):
    user_list.append(c['user_name'])

unique_users = set(user_list)

for idx,tu in enumerate(unique_users):
    batch = classifications.find({'user_name':tu,'subjects.zooniverse_id':{'$exists':'true'},'created_at':{'$gt':recent_date}},{'created_at':1,'subjects.zooniverse_id':1,'annotations':1})
    
    zlist,dlist,alist = [],[],[]
    for cc in list(batch):
        zlist.append(cc['subjects'][0]['zooniverse_id'])
        dlist.append(cc['created_at'])
        alist.append(cc['annotations'])
    
    d = {'created_at':dlist,'zooniverse_id':zlist,'annotations':alist}
    dfa = pd.DataFrame(d) 
    
    dups = dfa.duplicated(cols='zooniverse_id')
    dup_counts = dfa[dups]['zooniverse_id'].value_counts()
    dup2_zid = dup_counts[dup_counts == 1].index
    
    count_nodups = 0
    nodup_users = []
    if len(dup2_zid) > 0:
        dfa_2dups = dfa[dfa['zooniverse_id'].isin(dup2_zid)] 
        
        # Sort on zooniverse_id
        
        c1 = dfa_2dups.sort(columns='zooniverse_id')['created_at'][::2]
        c2 = dfa_2dups.sort(columns='zooniverse_id')['created_at'][1::2]
        td = c1.values - c2.values
        
        td_log = np.log10(np.abs(td*1e-9).astype(float) + 1)
        
        print '%30s had %5i duplicates out of %5i classifications; avg. gap was %5.3f seconds' % (tu,len(dup2_zid),len(dfa),np.mean(td)*1e-9)

    else:
        count_nodups += 1
        nodup_users.append(tu)

    # Was the annotation on the duplicate classification the same?

    if len(dup2_zid) == 1:
        dfa_2dups = dfa[dfa['zooniverse_id'].isin(dup2_zid)] 
        
        # Sort on zooniverse_id
        
        a1 = dfa_2dups.sort(columns='zooniverse_id')['annotations'][::2]
        a2 = dfa_2dups.sort(columns='zooniverse_id')['annotations'][1::2]
        if a1 != a2:
            print a1,a2
        
        print '%30s annotation is different: %s, %s' % (tu,a1,a2)



print '\n%i users had no duplicates.' % count_nodups
print nodup_users


