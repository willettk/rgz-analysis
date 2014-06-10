import rgz
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
from matplotlib import rc

subjects,classifications,users = rgz.load_rgz_data()

# Try for range of users

top_users = ('antikodon','planetari7','pamelaann','WizardHowl','JeanTate','Dolorous Edd','DocR','xDocR','KWillett')

fig = plt.figure(1,(8,7))
fig.clf()

javascript_fix_date = datetime.datetime(2014, 6, 1, 0, 0, 0, 0)

for idx,tu in enumerate(top_users):
    batch = classifications.find({'user_name':tu,'subjects.zooniverse_id':{'$exists':'true'},'created_at':{'$gt':javascript_fix_date}},{'created_at':1,'subjects.zooniverse_id':1})
    #batch = classifications.find({'user_name':tu,'subjects.zooniverse_id':{'$exists':'true'}},{'created_at':1,'subjects.zooniverse_id':1})
    
    zlist,dlist = [],[]
    for cc in list(batch):
        zlist.append(cc['subjects'][0]['zooniverse_id'])
        dlist.append(cc['created_at'])
    
    d = {'created_at':dlist,'zooniverse_id':zlist}
    dfa = pd.DataFrame(d) # 41901 classifications
    
    dups = dfa.duplicated(cols='zooniverse_id')
    dup_counts = dfa[dups]['zooniverse_id'].value_counts()
    dup2_zid = dup_counts[dup_counts == 1].index
    
    dfa_2dups = dfa[dfa['zooniverse_id'].isin(dup2_zid)] # 5338 classifications
    
    # Sort on zooniverse_id
    
    c1 = dfa_2dups.sort(columns='zooniverse_id')['created_at'][::2]
    c2 = dfa_2dups.sort(columns='zooniverse_id')['created_at'][1::2]
    td = c1.values - c2.values
    
    td_log = np.log10(np.abs(td*1e-9).astype(float) + 1)
    
    ax = fig.add_subplot(3,3,idx+1)
    timerange = (0,3)
    ax.hist(td_log,bins=25,range=timerange)
    ax.set_xlim(timerange[0],timerange[1])
    ax.set_xlabel(r'log $\Delta t$ [sec]')
    ax.set_ylabel('Counts')
    ax.set_title(r'%s - $N_{dup}=$%i' % (tu,len(dfa_2dups)),fontsize=10)
    rc(('xtick','ytick'), labelsize=10)

    print '%s\n' % tu

fig.set_tight_layout(True)
plt.savefig('/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis/plots/dups_since_jun07.png')
plt.clf()


