# Processing the RGZ catalog to deal with the fact that some images have the same radio galaxy, but centered
# on different positions

from astropy.io import fits
from consensus import rgz_path,data_path,plot_path,db
import numpy as np
import itertools

internal = fits.getdata("{0}/internal_matches.fits".format(rgz_path),1)

catalog = db['catalog']

# Test case

def coomatch(ra0,dec0,ra1,dec1,tol_ra=5e-4,tol_dec=5e-5):

    if (abs(ra0 - ra1) < tol_ra) & (abs(dec0 - dec1) < tol_dec):
        return True
    else:
        return False

imax = internal['GroupID'].max()

answers = []
for m in range(imax):
    i = (internal['GroupID'] == m+1)
    answer = False
    
    zooniverse_ids = internal[i]['zooniverse_id']

    for z1,z2 in itertools.combinations(zooniverse_ids,2):
    
        c1 = catalog.find({'zooniverse_id':z1}).count()
        c2 = catalog.find({'zooniverse_id':z2}).count()

        if (c1 > 0) & (c2 > 0):

            c1 = catalog.find_one({'zooniverse_id':z1})
            c2 = catalog.find_one({'zooniverse_id':z2})

            if c1['radio']['numberComponents'] != c2['radio']['numberComponents']:
                #print "Number of radio components in {0},{1} is different.".format(z1,z2)
                pass
            else:
            
                # For all possible orders of the components in both lists, 
                # is there a set that matches within the tolerance?
            
                # (c1['radio']['components'])
            
                c2perm = itertools.permutations(c2['radio']['components'])
                for x in c2perm:
                    results = True
                    for yy,xx in zip(c1['radio']['components'],x):
            
                        c2_ra0,c2_ra1 = xx['raRange']
                        c2_dec0,c2_dec1 = xx['decRange']
                        c1_ra0,c1_ra1 = yy['raRange']
                        c1_dec0,c1_dec1 = yy['decRange']
                
                        top = coomatch(c2_ra0,c2_dec0,c1_ra0,c1_dec0)
                        bot = coomatch(c2_ra1,c2_dec1,c1_ra1,c1_dec1)
                
                        if top and bot:
                            results &= True
                        else:
                            results &= False
            
                    if results:
                        #print "Match found for GroupID {0}".format(m+1)
                        answer |= True
    answers.append(answer)

print sum(answers),len(answers),len(internal)

# 10778 40270 177461
