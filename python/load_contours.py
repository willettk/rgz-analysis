'''
Find paths for RGZ data and use it to load raw data such as contour files

Looks in local directories first; if not found, sends request over network to AWS

Kyle Willett, UMN
'''

import os
import json
import requests

def make_pathdict(local=False):

    # Get dictionary for finding the path to FITS files and WCS headers

    paths = ['/Volumes/REISEPASS', \
             '/Volumes/3TB', \
             '/data/extragal/willett']

    data_found = False
    for dp in paths:
        if os.path.exists(dp):
            rgz_filepath = "{0:s}/rgz".format(dp)
            img_filepath = "{0:s}/raw_images".format(rgz_filepath)
            print "Using RGZ data from {0:s}".format(dp)
            data_found = True
            break 
    if not data_found:
        print "No external drives found for RGZ data"
        return None

    pathdict = {}

    # FIRST

    with open('%s/first_fits.txt' % rgz_filepath) as f:
        lines = f.readlines()
    
    for l in lines:
        spl = l.split(' ')
        dirno = int(spl[0])
        source = spl[1].strip()
        pathdict[source] = {}
        pathdict[source]['contours'] = '%s/RGZ-full.%i/CONTOURS/%s.json' % (img_filepath,dirno,source)
        pathdict[source]['radio'] = '%s/RGZ-full.%i/FIRST-IMGS/%s.fits' % (img_filepath,dirno,source)

    # ATLAS

    with open('%s/atlas_subjects.txt' % rgz_filepath) as f:
        lines = f.readlines()
    
    for l in lines:
        spl = l.rstrip()
        pathdict[spl] = {}
        pathdict[spl]['contours'] = '%s/ATLAS/CONTOURS/%s.json' % (img_filepath,spl)
        pathdict[spl]['radio'] = '%s/ATLAS/2x2/%s_radio.fits' % (img_filepath,spl)
        pathdict[spl]['ir'] = '%s/ATLAS/2x2/%s_ir.fits' % (img_filepath,spl)

    return pathdict

def get_contours(subject,pathdict):

    # Return the radio contours for an RGZ subject, input as a dictionary from MongoDB

    source = subject['metadata']['source']
    
    try:
        with open(pathdict[source]['contours']) as f:
            contours = json.load(f)
    # If pathdict is None, try to download over network
    except TypeError:
        r = requests.get(subject['location']['contours'])
        contours = r.json()
    
    return contours

