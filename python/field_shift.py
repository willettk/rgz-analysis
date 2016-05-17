# Processing the RGZ catalog to deal with the fact that some images have the
# same radio galaxy, but centered on different positions

from astropy.io import fits
from consensus import rgz_path,db
import itertools
import pprint

internal = fits.getdata("{0}/fits/internal_matches.fits".format(rgz_path),1)

catalog = db['catalog']

# Test case

def coordinate_match(ra0,dec0,ra1,dec1,tol_ra=5e-4,tol_dec=5e-5):

    if (abs(ra0 - ra1) < tol_ra) & (abs(dec0 - dec1) < tol_dec):
        return True
    else:
        return False

def component_match(c1_components,c2_components):

    # Bounding box surrounding a radio component in Image 2
    c2_ra0,c2_ra1 = c2_components['raRange']
    c2_dec0,c2_dec1 = c2_components['decRange']
    # Bounding box surrounding a radio component in Image 1
    c1_ra0,c1_ra1 = c1_components['raRange']
    c1_dec0,c1_dec1 = c1_components['decRange']
    
    # Does the lower left corner of the component match for Image 1 and Image 2?
    ll_match = coordinate_match(c2_ra0,c2_dec0,c1_ra0,c1_dec0)
    # Does the upper right corner of the component match for Image 1 and Image 2?
    ur_match = coordinate_match(c2_ra1,c2_dec1,c1_ra1,c1_dec1)

    return ll_match,ur_match
                
def append_overlaps(c1,c2,field):

    newField = 'overlapImages.{0}'.format(field)
    c1_id = c1['catalog_id']
    c2_id = c2['catalog_id']

    if c1.has_key('overlapImages'):
        print "Has key overlapImages"
        if c1['overlapImages'].has_key(field):
            print "Has {0}".format(newField)
            if c2_id not in c1['overlapImages'][field]:
                print "ID already exists in list, appending {0}".format(c2_id)
                catalog.update({"zooniverse_id": z1},{'$addToSet':{newField:c2_id}},multi=True)
        else:
            print "Didn't have {2}; creating and appending {0} to {1}".format(c2_id,c1_id,newField)
            catalog.update({"zooniverse_id": z1},{'$set':{newField:[c2_id,]}},multi=True)
            print "{0} didn't have overlapImages key".format(c1_id)

    else:
        print "Didn't have {2}; creating and appending {0} to {1}".format(c2_id,c1_id,newField)
        catalog.update({"zooniverse_id": z1},{'$set':{newField:[c2_id,]}},multi=True)

    return None

threes = internal[internal['GroupSize'] == 3]
imax = threes['GroupID'].max()
#imax = 5

answers = []
for m in range(imax):
    i = (threes['GroupID'] == m+1)
    answer = False
    
    zooniverse_ids = threes[i]['zooniverse_id']

    for z1,z2 in itertools.combinations(zooniverse_ids,2):
    
        count1 = catalog.find({'zooniverse_id':z1}).count()
        count2 = catalog.find({'zooniverse_id':z2}).count()
    
        if (count1 > 0) & (count2 > 0):
    
            c1 = catalog.find_one({'zooniverse_id':z1})
            c2 = catalog.find_one({'zooniverse_id':z2})

            c1_id = c1['catalog_id']
            c2_id = c2['catalog_id']

            # Do any of the components match?

            anymatch = False
            for c1r,c2r in itertools.product(c1['radio']['components'],c2['radio']['components']):
                ll_match,ur_match = component_match(c1r,c2r)
                anymatch |= (ll_match and ur_match)

            if anymatch:
                append_overlaps(c1,c2,'shareComponents')
    
            # Test on two-component images only. TO BE REMOVED.
            if (c1['radio']['numberComponents'] == c2['radio']['numberComponents']) and (c1['radio']['numberComponents'] == 2):
            #if c1['radio']['numberComponents'] == c2['radio']['numberComponents']:
            
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
                        if ll_match and ur_match:
                            results &= True
                        else:
                            results &= False
            
                    if results:
                        print "Match found for GroupID {0}: images {1} and {2}".format(m+1,c1_id,c2_id)
                        answer |= True

                        # Record the answer by appending it to the catalog entry
                        # Only do it once, otherwise it'll repeat when the second one rolls around.

                        append_overlaps(c1,c2,'matchComponents')

                        # Don't need to bother doing the rest of the permutations if we find an answer
                        break

    if answer:
        # Breaks when it finds set of two matching components
        break
        #answers.append(answer)

# print sum(answers),len(answers),len(internal)
# 10778 40270 177461
