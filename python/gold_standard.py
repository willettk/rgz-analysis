# How many galaxies from the gold standard sample are necessary to establish whether a user agrees with the experts?

# Hm - not many have done more than 1. Might have to expand the experiment to the full expert data set of 100 galaxies ...

# Plot (N/2 - M/2) vs. (M+N) and see where agreement lies!!

# Read in the list of gold standard images

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
import consensus_4panel as cons
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
    
def top_volunteers(gs):

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
            ud[k] = v
    
    return ud

def one_answer(zid,user_name):

    sub = subjects.find_one({'zooniverse_id':zid})
    imgid = sub['_id']
    zid = sub['zooniverse_id']

    # Classifications for this subject after launch date
    clist = list(classifications.find({"subject_ids": imgid, "updated_at": {"$gt": cons.main_release_date},'user_name':user_name}))
    
    # Empty dicts and lists 
    cdict = {}
    checksum_list = []
    
    for c in clist:
        # Want most popular combination for each NUMBER of galaxies identified in image
        
        sumlist = []    # List of the checksums over all possible combinations
        # Only find data that was an actual marking, not metadata
        goodann = [x for x in c['annotations'] if x.keys()[0] not in cons.bad_keys]
        n_galaxies = len(goodann)
    
        for idx,ann in enumerate(goodann):
    
            xmaxlist = []
            radio_comps = ann['radio']
    
            # loop over all the radio components within an galaxy
            if radio_comps != 'No Contours':
                for rc in radio_comps:
                    xmaxlist.append(float(radio_comps[rc]['xmax']))
            # or make the value -99 if there are no contours
            else:
                xmaxlist.append(-99)
    
            # To create a unique ID for the combination of radio components,
            # take the product of all the xmax coordinates and sum them together.
            product = reduce(operator.mul, xmaxlist, 1)
            sumlist.append(round(product,3))
    
        checksum = sum(sumlist)
        checksum_list.append(checksum)
        c['checksum'] = checksum
    
        # Insert checksum into dictionary with number of galaxies as the index
        if cdict.has_key(n_galaxies):
            cdict[n_galaxies].append(checksum)
        else:
            cdict[n_galaxies] = [checksum]
    
    #print cdict,'\n'
    
    maxval=0
    mc_checksum = 0.
    ngals = 0

    # Find the number of galaxies that has the highest number of consensus classifications

    for k,v in cdict.iteritems():
        mc = Counter(v).most_common()
        # Check if the most common selection coordinate was for no radio contours
        if mc[0][0] == -99.0:
            if len(mc) > 1:
                # If so, take the selection with the next-highest number of counts
                mc_best = mc[1]
            else:
                continue
        # Selection with the highest number of counts
        else:
            mc_best = mc[0]
        # If the new selection has more counts than the previous one, choose it as the best match;
        # if tied or less than this, remain with the current consensus number of galaxies
        if mc_best[1] > maxval:
            maxval = mc_best[1]
            mc_checksum = mc_best[0]
            ngals = k
    
    # Find a galaxy that matches the checksum (easier to keep track as a list)
    
    cmatch = next(i for i in clist if i['checksum'] == mc_checksum)
   
    # Find IR peak for the checksummed galaxies
    
    goodann = [x for x in cmatch['annotations'] if x.keys()[0] not in cons.bad_keys]

    # Find the sum of the xmax coordinates for each galaxy. This gives the index to search on.
    
    consensus = {}
    ir_x,ir_y = {},{}
    for k,gal in enumerate(goodann):
        xmax_temp = []
        try:
            for v in gal['radio'].itervalues():
                xmax_temp.append(float(v['xmax']))
        except AttributeError:
            xmax_temp.append(-99)

        checksum2 = round(sum(xmax_temp),3)
        consensus[checksum2] = {}
        consensus[checksum2]['ind'] = k
        consensus[checksum2]['xmax'] = xmax_temp
    
        # Make empty copy of next dict in same loop
        ir_x[k] = []
        ir_y[k] = []
    
    # Now loop over the galaxies themselves
    for c in clist:
        if c['checksum'] == mc_checksum:
    
            annlist = [ann for ann in c['annotations'] if ann.keys()[0] not in cons.bad_keys]
            for ann in annlist:
                if 'ir' in ann.keys():
                    # Find the index k that this corresponds to
                    try:
                        xmax_checksum = round(sum([float(ann['radio'][a]['xmax']) for a in ann['radio']]),3)
                    except TypeError:
                        xmax_checksum = -99
                    k = consensus[xmax_checksum]['ind']
    
                    if ann['ir'] == 'No Sources':
                        ir_x[k].append(-99)
                        ir_y[k].append(-99)
                    else:
                        # Only takes the first IR source right now; NEEDS TO BE MODIFIED.
    
                        ir_x[k].append(float(ann['ir']['0']['x']))
                        ir_y[k].append(float(ann['ir']['0']['y']))


                    for k,v in consensus.iteritems():
                        if v['ind'] == k:
                            consensus[k]['ir'] = (xpeak,ypeak)
        
    return consensus

def compare(zid,user_name,verbose=True):

    their_answer = one_answer(zid,user_name)
    expert_answer = cons.checksum(zid,save_fig=True)

    nt = len(their_answer)
    nx = len(expert_answer)

    if nt == nx:
        match_temp = True
        for tk in their_answer.iterkeys():
            if expert_answer.has_key(tk):
                match_temp = True and expert_answer.has_key(tk)

        if verbose:
            if match_temp:
                print "%s's answer for %s (%i sources) matches that of the experts." % (user_name,zid,nt)
            else:
                print "%s's answer for %s (%i sources) does not match that of the experts." % (user_name,zid,nt)

        match_expert = match_temp
    else:
        if verbose:
            print "%s's number of sources for %s does not match (%i vs %i)" % (user_name,zid,nt,nx)
        match_expert = False

    return match_expert

def all_compare(gs,user_name):

    matches = []
    for zid in gs:

        sub = subjects.find_one({'zooniverse_id':zid})
        imgid = sub['_id']

        c = classifications.find_one({"subject_ids": imgid, "updated_at": {"$gt": cons.main_release_date},'user_name':user_name})
        
        if c is not None:
            match_expert = compare(zid,user_name,verbose=False)
            matches.append(match_expert)

    return matches

def ac_data():

    gs = get_galaxies()
    ud = top_volunteers(gs)

    match_ratio = []
    for user_name in ud:
        matches = all_compare(gs,user_name)
        match_ratio.append(np.sum(matches,dtype=float)/len(matches))

    with open('%s/goldstandard/gs_compare.pkl' % rgz_dir, 'wb') as output:
        pickle.dump(match_ratio, output)
    
    with open('%s/goldstandard/gs_topusers.pkl' % rgz_dir, 'wb') as output:
        pickle.dump(ud, output)
    
    return None

def barchart():

    with open('%s/goldstandard/gs_compare.pkl' % rgz_dir, 'rb') as pkl_file:
        data = pickle.load(pkl_file)

    with open('%s/goldstandard/gs_topusers.pkl' % rgz_dir, 'rb') as output:
        ud = pickle.load(output)

    ind = np.arange(len(data))  # the x locations for the groups
    width = 0.75       # the width of the bars
    c = 0.683
    n = ud.values()
    k = [int(nn*mr) for nn,mr in zip(n,data)]
    p_lower = [dist.beta.ppf((1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    p_upper = [dist.beta.ppf(1-(1-c)/2., kk+1,nn-kk+1) for kk,nn in zip(k,n)]
    err_lower = [d-pl for d,pl in zip(data,p_lower)]
    err_upper = [pu-d for d,pu in zip(data,p_upper)]

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, data, width,yerr=[err_lower,err_upper])
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Percentage matched with experts')
    ax.set_title('Volunteer agreement on expert sample of 100 galaxies')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ud.keys(),rotation=45,ha='right' )
    ax.set_ylim(0,1)
    
    #ax.legend( (rects1[0],), ('RGZ user',) )
    
    for rect,ngal in zip(rects1,ud.itervalues()):
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%i'%int(ngal),
                ha='center', va='bottom',fontsize=8)
    
    fig.subplots_adjust(bottom=0.2)
    plt.show() 

    return None
