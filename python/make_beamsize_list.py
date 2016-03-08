# Make a CSV file of the beam sizes for all the Radio Galaxy Zoo sources

# Get list of the subjects from MongoDB

import rgz

subjects,classifications = rgz.load_rgz_data()

first_ids = []
for s in subjects.find():
    first_ids.append(s['metadata']['source'])

atlas_ids = [x for x in first_ids if x[0] == 'C']

# Kluge. Make sure tutorial subject is also cut out. 
first_ids = first_ids[:175000]

'''
FIRST data
'''

# Match to the appropriate header in the FITS files

fitslist_file = '/Volumes/3TB/rgz/first_fits.txt'

from astropy.io import fits

with open(fitslist_file,'r') as f:
    lines = f.readlines()

biglist_dirno = [int(x.split()[0]) for x in lines]
biglist_first = [x.split()[1] for x in lines]

with open('/Users/willettk/Astronomy/Research/GalaxyZoo/radiogalaxyzoo/rgz_beamsizes.txt','w') as writefile:

    for fi in first_ids:
        path = '/Volumes/3TB/rgz/raw_images/RGZ-full.%i/FIRST-IMGS' % biglist_dirno[biglist_first.index(fi)]
        with fits.open('%s/%s.fits' % (path,fi)) as f:
            hdr = f[0].header
    
        # Write as CSV

        bmaj,bmin = hdr['BMAJ']*3600,hdr['BMIN']*3600
        bpa = float(hdr['HISTORY'][10].split()[-1])
        s = subjects.find_one({'metadata.source':fi})
        print >> writefile,s['_id'],fi,s['zooniverse_id'],bmaj,bmin,bpa

    '''
    ATLAS data
    '''

    for ai in atlas_ids:
        path = '/Volumes/3TB/rgz/raw_images/ATLAS/2x2'
        with fits.open('%s/%s_radio.fits' % (path,ai)) as f:
            hdr = f[0].header
    
        # Write as CSV

        bmaj,bmin,bpa = hdr['BMAJ']*3600,hdr['BMIN']*3600,hdr['BPA']
        s = subjects.find_one({'metadata.source':ai})
        print >> writefile,s['_id'],'%8s' % ai,s['zooniverse_id'],bmaj,bmin,bpa
    
