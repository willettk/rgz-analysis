from astropy import coordinates as coord, units as u
import itertools, os
from consensus import db, rgz_path

catalog_0 = db['catalog']
catalog_5 = db['catalog_5']
catalog_10 = db['catalog_10']
catalog_20 = db['catalog_20']
catalogs = {'catalog_0':catalog_0, 'catalog_5':catalog_5, 'catalog_10':catalog_10, 'catalog_20':catalog_20}

def cids_for_removal():
    
    # Find duplicate sources (so ones not in the static catalog)
    
    cids = {}
    
    for weight in [0, 5, 10, 20]:
        cids['cids_{}'.format(weight)] = []

        # Old naming convention
        if catalogs['catalog_{}'.format(weight)].find_one({'duplicateSources':{'$exists':True}}):
            for c in catalogs['catalog_{}'.format(weight)].find({'duplicateSources.exactDuplicate':{'$exists':True}}):
                if c['catalog_id'] != min(c['duplicateSources']['exactDuplicate']):
                    cids['cids_{}'.format(weight)].append(c['catalog_id'])

        # New naming convention
        elif catalogs['catalog_{}'.format(weight)].find_one({'duplicate_sources':{'$exists':True}}):
            for c in catalogs['catalog_{}'.format(weight)].find({'duplicate_sources.exact_duplicate':{'$exists':True}}):
                if c['catalog_id'] != min(c['duplicate_sources']['exact_duplicate']):
                    cids['cids_{}'.format(weight)].append(c['catalog_id'])
    
    return cids

######### Begin IR block #########

def get_ir_mismatch(cids_unwt, cids_wt, weight):
    
    radio_pos_unwt, radio_pos_wt = set(), set()
    
    # Get set of all radio positions in each catalog
    print '    Finding radio positions'
    for c in catalog_0.find({'catalog_id':{'$nin':cids_unwt}}):
        radio_pos_unwt.add(c['rgz_name'])
    for c in catalogs['catalog_{}'.format(weight)].find({'catalog_id':{'$nin':cids_wt}}):
        radio_pos_wt.add(c['rgz_name'])
    
    # Find radio positions that are only in weighted catalogs
    unique_to_wt = radio_pos_wt.difference(radio_pos_unwt)
    
    # Iterate over every entry in the unweighted catalog
    no_counterpart_unwt_wt, ir_mismatch = compare_ir(names=radio_pos_unwt, cat_1=catalog_0, cat_2=catalogs['catalog_{}'.format(weight)], \
                                                     cids_1=cids_unwt, cids_2=cids_wt)
    
    # Iterate over every entry that's only in the weighted catalogs
    no_counterpart_wt_unwt = no_counterparts(names=unique_to_wt, cat_1=catalogs['catalog_{}'.format(weight)], cids_1=cids_wt)
    
    return ir_mismatch, no_counterpart_unwt_wt, no_counterpart_wt_unwt

def get_ir(c):
    if 'IR_ra' in c['consensus']:
        return coord.SkyCoord(c['consensus']['IR_ra'], c['consensus']['IR_dec'], unit=(u.deg,u.deg), frame='icrs')
    else:
        return coord.SkyCoord(0, 0, unit=(u.deg,u.deg), frame='icrs')

def compare_ir(names, cat_1, cat_2, cids_1, cids_2):
    
    # Iterate over every entry in catalog 1 for matches in catalog 2
    # This returns the lists no_counterpart and ir_mismatch
    
    print '    Comparing {} and {} IR positions'.format(cat_1.name, cat_2.name)
    count = 0
    
    no_counterpart, ir_mismatch = [], []
    
    for rgz_name in names:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(names))
        
        # Get all entries with this position from each catalog
        c1s, c2s = [], []
        for c1 in cat_1.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_1}}):
            c1s.append(c1)
        for c2 in cat_2.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_2}}):
            c2s.append(c2)
        
        # If the weighted catalog has no match, log in 'no_counterpart'
        if not len(c2s):
            for c1 in c1s:
                no_counterpart.append(c1['zooniverse_id'])
            
        # Otherwise check each combination for an IR match
        else:
            for c1, c2 in itertools.product(c1s, c2s):
                ir_1 = get_ir(c1)
                ir_2 = get_ir(c2)
                
                # When a pair is found with same IR position, remove them from consideration
                if ir_1.separation(ir_2).arcsecond <= 3.:
                    try:
                        c1s.remove(c1)
                    except ValueError:
                        pass
                    try:
                        c2s.remove(c2)
                    except ValueError:
                        pass
            
            # For each combination that didn't match, log in 'ir_mismatch'
            for c1, c2 in itertools.product(c1s, c2s):
                ir_mismatch.append((c1['catalog_id'], c2['catalog_id']))
    
    return no_counterpart, ir_mismatch

def no_counterparts(names, cat_1, cids_1):
    
    # Iterate over every entry that's only in the weighted catalogs
    # These can just be logged in 'no_counterpart' without checking (by definition)
    
    print '    No counterpart between catalog and {}'.format(cat_1.name)
    count = 0
    
    no_counterpart = []
    
    for rgz_name in names:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(names))
        
        for c1 in cat_1.find({'rgz_name':rgz_name, 'catalog_id':{'$nin':cids_1}}):
            no_counterpart.append(c1['zooniverse_id'])
    
    return no_counterpart

def print_ir(filename, ir_mismatch):
    
    # Print the list of IR mismatches to a CSV
    
    with open(filename, 'w') as f:
        print >> f, 'unweighted_catalog_id,weighted_catalog_id'
        for c1, c2 in ir_mismatch:
            print >> f, '{},{}'.format(c1, c2)
    
    return None
                
######### End IR block; begin radio block #########

def get_radio_mismatch(no_counterpart_1, no_counterpart_2, cat_1, cat_2, cids_1, cids_2):
    
    print '    Find radio mismatch between {} and {}'.format(cat_1.name, cat_2.name)
    count = 0
    
    zids = set().union(no_counterpart_1, no_counterpart_2)
    
    # Find each source in each catalog from these subjects
    czids_1, czids_2 = {}, {}
    for zid in zids:
        
        count += 1
        if not count%1000:
            print '        {}/{}'.format(count, len(zids))
        
        for c1 in cat_1.find({'zooniverse_id':zid, 'catalog_id':{'$nin':cids_1}}):
            if zid in czids_1:
                czids_1[zid].append(c1['catalog_id'])
            else:
                czids_1[zid] = [c1['catalog_id']]
        
        for c2 in cat_2.find({'zooniverse_id':zid, 'catalog_id':{'$nin':cids_2}}):
            if zid in czids_2:
                czids_2[zid].append(c2['catalog_id'])
            else:
                czids_2[zid] = [c2['catalog_id']]
    
    return czids_1, czids_2

def print_radio(filename, czid_unweighted, czid_weighted, weight):
    
    # Print the subject fields that contain radio sources without matches
    
    zids = set().union(czid_unweighted.keys(), czid_weighted.keys())
    
    with open(filename, 'w') as f:
        print >> f, 'zooniverse_id,catalog_id,catalog_weight'
        for zid in zids:
            if zid in czid_unweighted:
                for cid in czid_unweighted[zid]:
                    print >> f, '{},{},0'.format(zid, cid)
            if zid in czid_weighted:
                for cid in czid_weighted[zid]:
                    print >> f, '{},{},{}'.format(zid, cid, weight)
    
    return None

######### End radio block #########

if __name__ == "__main__":

    output_dir = '{}/csv/weighting_checks'.format(rgz_path)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print 'Finding duplicate entries'
    cids = cids_for_removal()
    print '    Removed {}, {}, {}, and {} duplicate entries'.format(len(cids['cids_0']), len(cids['cids_5']), len(cids['cids_10']), len(cids['cids_20']))

    for weight in [5, 10, 20]:

        print 'Testing catalog with weight {}'.format(weight)

        # IR
        print 'Finding IR mismatches'
        ir_mismatch, no_counterpart_unwt_wt, no_counterpart_wt_unwt = get_ir_mismatch(cids_unwt=cids['cids_0'], cids_wt=cids['cids_{}'.format(weight)], weight=weight)
        ir_filename = '{}/ir_mismatch_{}.csv'.format(output_dir, weight)
        print 'IR checks complete; printing {} IR mismatches from catalog_{} to {}'.format(len(ir_mismatch), weight, ir_filename)
        print_ir(ir_filename, ir_mismatch)

        # Radio
        print 'Finding radio mismatches'
        czid_weighted, czid_unweighted = get_radio_mismatch(no_counterpart_1=no_counterpart_wt_unwt, no_counterpart_2=no_counterpart_unwt_wt, \
                                                            cat_1=catalogs['catalog_{}'.format(weight)], cat_2=catalogs['catalog_0'], \
                                                            cids_1=cids['cids_{}'.format(weight)], cids_2=cids['cids_0'])
        radio_filename = '{}/no_counterpart_{}.csv'.format(output_dir, weight)
        print 'Radio checks complete; printing {} radio mismatches from catalog_{} to {}'.format(len(czid_weighted)+len(czid_unweighted), weight, radio_filename)
        print_radio(radio_filename, czid_unweighted, czid_weighted, weight)
