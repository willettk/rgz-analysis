# Processing the RGZ catalog to deal with the fact that some images have the
# same radio galaxy, but centered on different positions

from astropy import coordinates as coord, units as u
from astropy.io import fits
from consensus import rgz_path, db, version
import itertools, logging

#contains groups of subjects within 3' of each other, determined in TopCat
internal = fits.getdata("{0}/fits/internal_matches.fits".format(rgz_path),1)

catalog = db['catalog{}'.format(version)]

# Test case

def coordinate_match(ra0,dec0,ra1,dec1,tol_ra=5e-4,tol_dec=5e-4):
    
    return (abs(ra0 - ra1) < tol_ra) & (abs(dec0 - dec1) < tol_dec)

def component_match(c1_components,c2_components):
    
    # Bounding box surrounding a radio component in Image 2
    c2_ra0,c2_ra1 = c2_components['ra_range']
    c2_dec0,c2_dec1 = c2_components['dec_range']
    # Bounding box surrounding a radio component in Image 1
    c1_ra0,c1_ra1 = c1_components['ra_range']
    c1_dec0,c1_dec1 = c1_components['dec_range']
    
    # Does the lower left corner of the component match for Image 1 and Image 2?
    ll_match = coordinate_match(c2_ra0,c2_dec0,c1_ra0,c1_dec0)
    # Does the upper right corner of the component match for Image 1 and Image 2?
    ur_match = coordinate_match(c2_ra1,c2_dec1,c1_ra1,c1_dec1)
    
    return ll_match,ur_match

def append_overlaps(c1,c2,field):
    
    newField = 'duplicate_sources.{0}'.format(field)
    c1_cid = c1['catalog_id']
    c2_cid = c2['catalog_id']
    
    if c1.has_key('duplicate_sources'):
        #print "{0} has key 'duplicateSources'".format(c1_cid)
        if c1['duplicate_sources'].has_key(field):
            #print "Has {0}".format(newField)
            if c2_cid not in c1['duplicate_sources'][field]:
                #print "Appending {0} to list".format(c2_cid)
                catalog.update({'catalog_id': c1_cid},{'$addToSet':{newField:c2_cid}})
        else:
            #print "Didn't have '{2}'; creating and appending {0} to {1}".format(c2_cid,c1_cid,newField)
            catalog.update({'catalog_id': c1_cid},{'$set':{newField:[c2_cid,]}})
        if field == 'match_components':
            catalog.update({'catalog_id': c1_cid},{'$pull':{'duplicate_sources.share_components':c2_cid}})
    
    else:
        #print "{1} didn't have key 'duplicateSources'; creating and appending {0} under '{2}'".format(c2_cid,c1_cid,field)
        catalog.update({'catalog_id': c1_cid},{'$set':{newField:[c2_cid,]}})
    
    if c2.has_key('duplicate_sources'):
        #print "{0} has key 'duplicateSources'".format(c2_cid)
        if c2['duplicate_sources'].has_key(field):
            #print "Has {0}".format(newField)
            if c1_cid not in c2['duplicate_sources'][field]:
                #print "Appending {0} to list".format(c1_cid)
                catalog.update({'catalog_id': c2_cid},{'$addToSet':{newField:c1_cid}})
        else:
            #print "Didn't have '{2}'; creating and appending {0} to {1}".format(c1_cid,c2_cid,newField)
            catalog.update({'catalog_id': c2_cid},{'$set':{newField:[c1_cid,]}})
        if field=='match_components':
            catalog.update({'catalog_id': c2_cid},{'$pull':{'duplicate_sources.share_components':c1_cid}})
    
    else:
        #print "{1} didn't have key 'duplicateSources'; creating and appending {0} under '{2}'".format(c1_cid,c2_cid,field)
        catalog.update({'catalog_id': c2_cid},{'$set':{newField:[c1_cid,]}})
    
    if field == 'exact_duplicate':
        catalog.update({'catalog_id': c1_cid},{'$addToSet':{newField:c1_cid}})
        catalog.update({'catalog_id': c2_cid},{'$addToSet':{newField:c2_cid}})
    
    return None

def check_symmetry(zooniverse_ids):
    returnVal = True
    for z in zooniverse_ids:
        for subject1 in catalog.find({ 'zooniverse_id':z, 'duplicate_sources':{'$exists':True} }):
            cid1 = subject1['catalog_id']
            for duplicateCategory in subject1['duplicate_sources']:
                for cid2 in subject1['duplicate_sources'][duplicateCategory]:
                    subject2 = catalog.find_one({'catalog_id':cid2})
                    if ('duplicate_sources' not in subject2) or \
                       (duplicateCategory not in subject2['duplicate_sources']) or \
                       (cid1 not in subject2['duplicate_sources'][duplicateCategory]):
                        print "{1}'s 'duplicate_sources.{2}' should contain {0} but doesn't".format(cid1, cid2, duplicateCategory)
                        returnVal = False
    #print "All symmetric"
    return returnVal

def find_duplicates(zid):
    
    groupID = internal[internal['zooniverse_id']==zid]['GroupID'][0]
    if groupID > 0:
        
        i = (internal['GroupID'] == groupID)
        zooniverse_ids = internal[i]['zooniverse_id']
        
        for z1,z2 in itertools.combinations(zooniverse_ids,2):
            
            count1 = catalog.find({'zooniverse_id':z1}).count()
            count2 = catalog.find({'zooniverse_id':z2}).count()
            
            if (count1 > 0) & (count2 > 0):
                for c1 in catalog.find({'zooniverse_id':z1}):
                    for c2 in catalog.find({'zooniverse_id':z2}):
                        
                        c1_id = c1['catalog_id']
                        c2_id = c2['catalog_id']
                        
                        # Do any of the components match?
                        
                        anymatch = False
                        for c1r,c2r in itertools.product(c1['radio']['components'],c2['radio']['components']):
                            ll_match,ur_match = component_match(c1r,c2r)
                            anymatch |= (ll_match and ur_match)
                        
                        if anymatch:
                            print "Shared components found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id)
                            logging.info("Shared components found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id))
                            append_overlaps(c1,c2,'share_components')
                            c1 = catalog.find_one({'catalog_id':c1_id})
                            c2 = catalog.find_one({'catalog_id':c2_id})
                        
                        if c1['radio']['number_components'] == c2['radio']['number_components']:
                            
                            # For all possible orders of the components in both lists, 
                            # is there a set that matches within the tolerance?
                            
                            c2permutations = itertools.permutations(c2['radio']['components'])
                            
                            # Loop over permutations of radio components
                            for c2perm in c2permutations:
                                results = True
                                
                                # Loop over individual components in main and comparison images
                                # and see if their positions match
                                
                                c1fix = c1['radio']['components']
                                for c1_components,c2_components in zip(c1fix,c2perm):
                                    
                                    ll_match,ur_match = component_match(c1_components,c2_components)
                                    
                                    # If both corners match, then this is the same component in Image 1 and Image 2. 
                                    results &= (ll_match and ur_match)
                                    
                                if results:
                                    print "Matching components found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id)
                                    logging.info("Matching components found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id))
                                    
                                    # Record the answer by appending it to the catalog entry
                                    # Only do it once, otherwise it'll repeat when the second one rolls around.
                                    
                                    append_overlaps(c1,c2,'match_components')
                                    c1 = catalog.find_one({'catalog_id':c1_id})
                                    c2 = catalog.find_one({'catalog_id':c2_id})
                                    
                                    #check if IR matches as well
                                    if 'ir_ra' not in c1['consensus'] and 'ir_ra' not in c2['consensus']:
                                        print "Exact match (no IR) found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id)
                                        logging.info("Exact match (no IR) found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id))
                                        append_overlaps(c1,c2,'exact_duplicate')
                                        c1 = catalog.find_one({'catalog_id':c1_id})
                                        c2 = catalog.find_one({'catalog_id':c2_id})
                                        
                                    elif 'ir_ra' in c1['consensus'] and 'ir_ra' in c2['consensus']:
                                        if ('AllWISE' in c1 and 'AllWISE' in c2 and c1['AllWISE']==c2['AllWISE']) or \
                                           ('AllWISE' not in c1 and 'AllWISE' not in c2):
                                            print "Exact match (same WISE catalog ID) found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id)
                                            logging.info("Exact match (same WISE catalog ID) found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id))
                                            append_overlaps(c1,c2,'exact_duplicate')
                                            c1 = catalog.find_one({'catalog_id':c1_id})
                                            c2 = catalog.find_one({'catalog_id':c2_id})
                                            
                                        else:
                                            c1_ir = coord.SkyCoord(c1['consensus']['ir_ra'], c1['consensus']['ir_dec'], unit=(u.deg,u.deg), frame='icrs')
                                            c2_ir = coord.SkyCoord(c2['consensus']['ir_ra'], c2['consensus']['ir_dec'], unit=(u.deg,u.deg), frame='icrs')
                                            if c1_ir.separation(c2_ir).arcsecond < 3.0:
                                                print "WISE catalog mismatch found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id)
                                                logging.info("WISE catalog mismatch found for GroupID {0}: sources {1} and {2}".format(groupID,c1_id,c2_id))
                                                append_overlaps(c1,c2,'WISE_cat_mismatch')
                                                c1 = catalog.find_one({'catalog_id':c1_id})
                                                c2 = catalog.find_one({'catalog_id':c2_id})
                                    
                                    # Don't need to bother doing the rest of the permutations if we find an answer
                                    break
        
        if not check_symmetry(zooniverse_ids):
            raise RuntimeError("zooniverse_ids {0} aren't symmetric".format(str(zooniverse_ids)))
