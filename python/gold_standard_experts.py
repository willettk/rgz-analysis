# How many galaxies from the gold standard sample are necessary to establish whether a user agrees with the experts?

# Hm - not many have done more than 1. Might have to expand the experiment to the full expert data set of 100 galaxies ...

# Plot (N/2 - M/2) vs. (M+N) and see where agreement lies!!

# Read in the list of gold standard images

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
import consensus
import rgz
import operator
import numpy as np
import cPickle as pickle
from matplotlib import pyplot as plt
import scipy.stats.distributions as dist
from collections import Counter

subjects,classifications,users = rgz.load_rgz_data()

def get_galaxies():
    
    #with open('%s/goldstandard/gs_zids.txt' % rgz_dir) as f:
    with open('%s/expert/expert_all_zooniverse_ids.txt' % rgz_dir) as f:
        gs = f.readlines()
    
    gals = [g.strip() for g in gs]

    return gals

def get_top_volunteers(gs):

    experts=("42jkb", "ivywong", "stasmanian", "klmasters", "Kevin", "akapinska", "enno.middelberg", "xDocR", "DocR", "vrooje", "KWillett")
    
    # How many unique users have done each galaxy in the expert sample?
    
    allusers = []
    
    for g in gs:
    
        s = subjects.find_one({'zooniverse_id':g})
        c = classifications.find({'subject_ids.0':s['_id']})
        
        ulist = []
        for cc in c:
            if cc.has_key('user_name'):
                ulist.append(cc['user_name'])
            else:
                ulist.append('Anonymous')
    
        setu = set(ulist)
    
        for u in setu:
            if u not in experts:
                allusers.append(u)
    
    cntr_users = Counter(allusers)
    
    # Limit to users who have done 10 or more of the expert sample
    
    ud = {}
    for k,v in cntr_users.iteritems():
        if v >= 10 and k != 'Anonymous':
            ud[k] = {}
            ud[k]['ngals'] = v
    
    return ud

def get_experts(gals):

    experts=("42jkb", "ivywong", "stasmanian", "klmasters", "akapinska", "enno.middelberg", "xDocR", "vrooje", "KWillett")

    ud = {}
    for ex in experts:
        ud[ex] = {}
        ud[ex]['ngals'] = 0
        for g in gals:
            s = subjects.find_one({'zooniverse_id':g})
            ccount = classifications.find({'subject_ids.0':s['_id'],'user_name':ex}).count()
            if ccount > 0:
                ud[ex]['ngals'] += 1
                
    return ud
    
def compare(zid,user_name,verbose=True,exclude_user=False):

    # Compare the classification of a single galaxy by a single user to that of the expert science team

    their_answer = consensus.one_answer(zid,user_name)
    if exclude_user:
        expert_answer = consensus.checksum(zid,experts_only=True,excluded=user_name)
    else:
        expert_answer = consensus.checksum(zid,experts_only=True)

    # Total number of galaxies identified in each answer
    nt = len(their_answer['answer'])
    nx = len(expert_answer['answer'])

    suff_x = '' if nx == 1 else 's'
    suff_t = '' if nt == 1 else 's'

    if nt == nx:
        match_temp = True
        for tk in their_answer['answer'].iterkeys():
            match_temp = True and expert_answer['answer'].has_key(tk)

        if verbose:
            if match_temp:
                print "Yes - %s's answer for %s (%i source%s) matches that of the experts." % (user_name,zid,nt,suff_t)
            else:
                print "No  - %s's answer for %s has the same number of source%s (%i) as the experts, but not the same arrangement." % (user_name,zid,suff_t,nt)

        match_expert = match_temp
    else:
        if verbose:
            print "No  - %s had %i source%s in %s, compared to %i source%s by the experts" % (user_name,nt,suff_t,zid,nx,suff_x)
        match_expert = False

    return match_expert

def all_compare(gals,user_name):

    # Run the comparison between a single volunteer and the expert RGZ science team for all galaxies in a list

    matches = []
    for zid in gals:

        sub = subjects.find_one({'zooniverse_id':zid})
        imgid = sub['_id']

        c = classifications.find_one({"subject_ids": imgid, "updated_at": {"$gt": consensus.main_release_date},'user_name':user_name})
        
        if c is not None:
            match_expert = compare(zid,user_name,verbose=False)
            matches.append(match_expert)

    return matches

def ac_data():

    # Compute and save the levels of agreement for volunteers and experts for all overlapping galaxies

    gals = get_galaxies()

    # Experts only

    udx = get_experts(gals)

    # Run comparison for all matching galaxies and sum answers together to get total success ratio
    for user_name in udx:
        matches = all_compare(gals,user_name)
        udx[user_name]['match_ratio'] = np.sum(matches,dtype=float)/len(matches)

    # Write dict to file.
    with open('%s/goldstandard/%s.pkl' % (rgz_dir,'gs_expert_results'), 'wb') as output:
        pickle.dump(udx, output)
    
    # Same again, but for the top volunteers

    udv = get_top_volunteers(gals)

    for user_name in udv:
        matches = all_compare(gals,user_name)
        udv[user_name]['match_ratio'] = np.sum(matches,dtype=float)/len(matches)

    with open('%s/goldstandard/%s.pkl' % (rgz_dir,'gs_topvol_results'), 'wb') as output:
        pickle.dump(udv, output)
    
    return None

def barchart(expert=True,savefig=False):

    # Load saved data from pickled files

    if expert:
        filename = 'gs_expert_results'
        title_stub = 'Expert agreement'
        ylabel_stub = 'rest of '
        barcolor = '#e41a1c'
        figname = 'experts'
    else:
        filename = 'gs_topvol_results'
        title_stub = 'Volunteer agreement with experts'
        ylabel_stub = ''
        barcolor = '#377eb8'
        figname = 'topvols'


    with open('%s/goldstandard/%s.pkl' % (rgz_dir,filename), 'rb') as pkl_file:
        data = pickle.load(pkl_file)

    ind = np.arange(len(data))  # the x locations for the groups
    width = 0.75                # the width of the bars

    # Unpack dictionary into arrays
    match_ratios = []
    usernames = []
    ngals = []
    for k,v in data.iteritems():
        usernames.append(k)
        match_ratios.append(v['match_ratio'])
        ngals.append(v['ngals'])


    # Parameters for binomial uncertainty bars (see Cameron et al. 2011)
    c = 0.683
    n = ngals
    k = [int(nn*mr) for nn,mr in zip(n,match_ratios)]
    p_lower = [dist.beta.ppf((1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    p_upper = [dist.beta.ppf(1-(1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    err_lower = [d-pl for d,pl in zip(match_ratios,p_lower)]
    err_upper = [pu-d for d,pu in zip(match_ratios,p_upper)]

    # Set up bar chart
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, match_ratios, width,yerr=[err_lower,err_upper],color=barcolor)
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Percentage matched with %sexperts' % ylabel_stub)
    ax.set_title('%s for sample of 100 galaxies' % title_stub)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(usernames,rotation=45,ha='right' )
    ax.set_ylim(0,1)
    
    # Plot number of total galaxies above bars
    
    for rect,ngal in zip(rects1,n):
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%i'%int(ngal),
                ha='center', va='bottom',fontsize=8)

    # Show the average level of agreement, weighted by total number of galaxies classified

    wavg = np.average(match_ratios,weights=n)
    ax.hlines(wavg,ax.get_xlim()[0],ax.get_xlim()[1],lw=2,linestyle='--',color='k')
    
    fig.subplots_adjust(bottom=0.2)
    if savefig:
        fig.savefig('%s/goldstandard/barchart_%s.png' % (rgz_dir,figname))
    else:
        plt.show() 

    return None

if __name__ == '__main__':

    #ac_data()
    barchart(expert=True,savefig=True)
    barchart(expert=False,savefig=True)
