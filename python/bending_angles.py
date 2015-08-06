# Measure bending angles in the RGZ data for double components
#
# Kyle Willett, UMN

from __future__ import division
import rgz
import consensus

import pandas as pd

from astropy.io import ascii,fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

import requests
from PIL import Image
import cStringIO
import urllib

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.path import Path
import matplotlib.patches as patches

from scipy.interpolate import griddata
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

import time
import random

from astroML.plotting import hist as histML
    
from collections import Counter
from pymongo.errors import CursorNotFound

rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
local_fits_path = '%s/RGZ-full.44/FIRST-IMGs' % rgz_dir

IMG_HEIGHT_OLD = 424.0			# number of pixels in the original JPG image along the y axis
IMG_WIDTH_OLD = 424.0			# number of pixels in the original JPG image along the x axis
IMG_HEIGHT_NEW = 500.0			# number of pixels in the downloaded JPG image along the y axis
IMG_WIDTH_NEW = 500.0			# number of pixels in the downloaded JPG image along the x axis
FITS_HEIGHT = 301.0			# number of pixels in the FITS image (?) along the y axis
FITS_WIDTH = 301.0			# number of pixels in the FITS image (?) along the x axis
FIRST_FITS_HEIGHT = 132.0       # number of pixels in the FITS image along the y axis
FIRST_FITS_WIDTH = 132.0       # number of pixels in the FITS image along the y axis

first_ir_scale_x = FIRST_FITS_WIDTH / IMG_WIDTH_NEW
first_ir_scale_y = FIRST_FITS_HEIGHT / IMG_HEIGHT_NEW

PIXEL_SIZE = 0.00016667#/3600.0		# the number of arcseconds per pixel in the FITS image
xmin = 1.
xmax = IMG_HEIGHT_NEW
ymin = 1.
ymax = IMG_WIDTH_NEW

subjects,classifications,users = rgz.load_rgz_data()

rgz_75_file = '%s/csv/consensus_rgz_75.csv' % rgz_dir

def get_doubles():

    # Load in all examples of double-lobed radio galaxies from RGZ

    rgz75 = ascii.read(rgz_75_file,format='csv')

    dblidx = (rgz75['n_radio'] == 2) & (rgz75['consensus_level'] >= 0.75)
    doubles = rgz75[dblidx]

    return doubles

def get_triples():

    # Find all sources with three radio components

    rgz75 = ascii.read(rgz_75_file,format='csv')

    trpidx = (rgz75['n_radio'] == 3) & (rgz75['consensus_level'] >= 0.75)
    triples = rgz75[trpidx]

    return triples

def get_singles_large():

    rgz75 = ascii.read(rgz_75_file,format='csv')

    sngidx = (rgz75['n_radio'] == 1) & (rgz75['consensus_level'] >= 0.75)
    singles_large = rgz75[sngidx]

    return singles_large

def make_pathdict():

    with open(rgz_75_file,'r') as f:
        zid_read = f.readlines()

    zooniverse_ids = [z.split()[0].strip() for z in zid_read]
    
    # Get dictionary for finding the path to FITS files and WCS headers
    
    with open('%s/first_fits.txt' % rgz_dir) as f:
        lines = f.readlines()
    
    pathdict = {}
    for l in lines:
        spl = l.split(' ')
        pathdict[spl[1].strip()] = '/Volumes/3TB/rgz/raw_images/RGZ-full.%i/FIRST-IMGS/%s.fits' % (int(spl[0]),spl[1].strip())

    return pathdict

def all_doubles_pixradio(doubles,pathdict):

    with open('%s/bending_angles/angles_double_pixradio.csv' % rgz_dir,'w') as f:
        print >> f,'zooniverse_id,bending_angle,position_angle'
        for double in doubles:
            irx,iry,radio_components = pix_convert(double,pathdict)
            if irx is not None:
                xc,yc,radio_centroids = pix_radio(irx,iry,radio_components)
                alpha = bending_angle(xc,yc,radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                alpha_deg = alpha * 180./np.pi
                phi = position_angle(xc,yc,radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                phi_deg = phi * 180./np.pi
                if alpha is not None:
                    print >> f,'%s,%.3f' % (double['zooniverse_id'],alpha_deg,phi_deg)

    return None

def all_doubles_plot(doubles,pathdict):
    
    '''
    Deprecated in favor of all_doubles_pixradio
    Distribution of angles is identical, but I'm still worried about small differences. Could be due to rounding errors?
    '''
    
    with open('%s/bending_angles/angles_double_plot.csv' % rgz_dir,'w') as f:
        for dbl in doubles:
            irx,iry,radio_components = pix_convert(dbl,pathdict)
            cons = consensus.checksum(dbl['zooniverse_id'])
            if irx is not None:
                xc,yc,radio_centroids = pix_radio(irx,iry,radio_components)
                answer = cons['answer']
                if len(answer) > 0: # At least one galaxy was identified
                    for idx,ans in enumerate(answer.itervalues()):
                        try:
                            xc = ans['ir_peak'][0]
                            yc = ans['ir_peak'][1]
                        except KeyError:
                            xc,yc = ans['ir']
                alpha = bending_angle(xc,yc,radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                if alpha is not None:
                    print >> f,alpha*180/np.pi

    return None

def dblid(doubles,zooniverse_id):

    dbl = doubles[doubles['zooniverse_id'] == zooniverse_id][0]

    return dbl

def pix_convert(galaxy,pathdict,local=False):

    # Convert IR coordinates from RA/dec into pixel

    sub = subjects.find_one({'zooniverse_id':galaxy['zooniverse_id']})

    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    radio_components = contours['contours']

    try:
        assert len(radio_components) > 1, \
            'Radio data only has %i component for %s' % (len(radio_components),galaxy['zooniverse_id'])
    except AssertionError:
        return None,None,None

    # Keep everything in pixel coordinates. Reverse what's done in consensus.py; 
    # how can I transform an RA/dec pair into the radio component x/y pair?

    # Convert the pixel coordinates into RA,dec using the WCS object from the header

    if local:
        hdulist = fits.open('%s/%s' % (local_fits_path,pathdict[galaxy['first_id']].split('/')[-1]))
    else:
        hdulist = fits.open(pathdict[galaxy['first_id']])

    w = wcs.WCS(hdulist[0].header)

    worldcrd = np.array([[galaxy['ra'],galaxy['dec']]],np.float_)
    pix = w.wcs_world2pix(worldcrd,0)

    irx,iry = pix[0]
    irx_first,iry_first = np.round(pix[0][0] / first_ir_scale_x), np.round(IMG_HEIGHT_NEW - pix[0][1] / first_ir_scale_y)

    return irx,iry,radio_components

def radio_to_wcs(galaxy,pathdict):

    sub = subjects.find_one({'zooniverse_id':galaxy['zooniverse_id']})

    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    radio_components = contours['contours']

    try:
        assert len(radio_components) > 1, \
            'Radio data only has %i component for %s' % (len(radio_components),galaxy['zooniverse_id'])
    except AssertionError:
        return None,None,None

    # Keep everything in pixel coordinates. Reverse what's done in consensus.py; 
    # how can I transform an RA/dec pair into the radio component x/y pair?

    # Convert the pixel coordinates into RA,dec using the WCS object from the header

    hdulist = fits.open(pathdict[galaxy['first_id']])
    w = wcs.WCS(hdulist[0].header)

    worldcrd = np.array([[galaxy['ra'],galaxy['dec']]],np.float_)
    pix = w.wcs_world2pix(worldcrd,0)

    irx,iry = pix[0]
    irx_first,iry_first = np.round(pix[0][0] / first_ir_scale_x), np.round(IMG_HEIGHT_NEW - pix[0][1] / first_ir_scale_y)

    return irx,iry,radio_components

def centroid(bbox):

    cent_x = np.median((bbox[0],bbox[2]))
    cent_y = np.median((bbox[1],bbox[3]))

    return cent_x,cent_y

def bending_angle(xc,yc,cx1,cy1,cx2,cy2):

    # Compute the bending angle between three points

    '''
    Points are:

        - xc,yc:    x,y center of IR counterpart
        - cx1,cy1:  x,y center of 1st radio lobe
        - cx2,cy2:  x,y center of 2nd radio lobe
    '''


    # Alternate method of computing bending angle

    r1 = np.array([cx1,cy1])
    r2 = np.array([cx2,cy2])
    center = np.array([xc,yc])

    r1diff = r1-center
    r2diff = r2-center

    r1len = np.hypot(r1diff[0],r1diff[1])
    r2len = np.hypot(r2diff[0],r2diff[1])

    alpha = np.arccos(np.dot(r1diff,r2diff) / (r1len*r2len))

    return alpha

def position_angle(xc,yc,cx1,cy1,cx2,cy2):

    # Compute the bending angle (in radians) between three points

    '''
    Points are:

        - xc,yc:    x,y center of bending angle
        - cx1,cy1:  x,y center of 1st component
        - cx2,cy2:  x,y center of 2nd component
    '''

    r1 = np.array([cx1,cy1])
    r2 = np.array([cx2,cy2])
    center = np.array([xc,yc])

    r12sum = (r1-center) + (r2-center)
    r12len = np.hypot(r12sum[0],r12sum[1])

    north = np.array([0,1])
    northlen = np.hypot(north[0],north[1])

    alpha = np.arccos(np.dot(r12sum,north) / (r12len*northlen))

    # Measure CCW from north

    if r12sum[0] > 0.:
        alpha = 2*np.pi - alpha

    return alpha

def pix_radio(irx,iry,radio_components):

    # Add centers of bounding boxes

    radio_centroids = []
    for comp in radio_components:
        bbox = comp[0]['bbox']

        '''
        cx = np.median((bbox[0],bbox[2])) / first_ir_scale_x
        cy = np.median(((bbox[1]/first_ir_scale_y),(bbox[3]/first_ir_scale_y)))
        '''

        '''
        cxu,cyu = centroid(bbox)
        '''

        cxu = np.median((bbox[0],bbox[2]))
        cyu = np.median((bbox[1],bbox[3]))

        cx,cy = cxu/first_ir_scale_x,cyu/first_ir_scale_y
        radio_centroids.append((cx,cy))

    xc = irx / first_ir_scale_x
    yc = IMG_HEIGHT_NEW - iry / first_ir_scale_y

    return xc,yc,radio_centroids

def load_angles(filename):

    with open('%s/bending_angles/%s.csv' % (rgz_dir,filename),'r') as f:
        angstr = f.readlines()

    ba = [float(x.split(',')[1]) for x in angstr[1:]]
    #pa = [float(x.split(',')[2]) for x in angstr[1:]]

    return ba#,pa

def plothist(savefig=False):

    angles_double_pixradio = load_angles('angles_double_pixradio')
    angles_triple_pixradio = load_angles('angles_triple_pixradio')
    angles_double_mps = load_angles('angles_multipeaked_singles')
    angles_triple_mps = load_angles('angles_multipeaked_singles_no_optical')

    fig = plt.figure(2,(8,8))
    ax1 = fig.add_subplot(111)
    
    '''
    c1 = '#e41a1c'
    c2 = '#a6cee3'
    c3 = '#377eb8'
    c4 = '#386cb0'
    '''
    
    c1 = '#377eb8'
    c2 = '#e41a1c'
    c3 = '#4daf4a'
    c4 = '#984ea3'

    histML(angles_double_pixradio, bins=20, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c1, range=(0,180),label='double multi-contour')
    histML(angles_triple_pixradio, bins=20, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c2, range=(0,180),label='triple multi-contour')
    histML(angles_double_mps,      bins=20, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c3, range=(0,180),label='double single-contour')
    histML(angles_triple_mps,      bins=20, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c4, range=(0,180),label='triple single-contour')
    ax1.set_xlim(0,180)
    ax1.vlines(x=np.median(angles_double_pixradio),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c1,linestyle='--')
    ax1.vlines(x=np.median(angles_triple_pixradio),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c2,linestyle='--')
    ax1.vlines(x=np.median(angles_double_mps),     ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c3,linestyle='--')
    ax1.vlines(x=np.median(angles_triple_mps),     ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c4,linestyle='--')
    ax1.set_xlabel(r'bending angle [deg]',fontsize=24)
    ax1.set_ylabel('count',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    
    ax1.legend(loc='upper left')
    fig.tight_layout()

    if savefig:
        fig.savefig('%s/bending_angles/plots/bending_angles_hist.pdf' % rgz_dir)
    else:
        plt.show()

    return None

def plot_one_double(dbl,pathdict,figno=1,save_fig=False):

    cons = consensus.checksum(dbl['zooniverse_id'])

    '''
    if pathdict[dbl['first_id']].split('/')[5].split('.')[-1] == '35':
        irx,iry,radio_components = pix_convert(dbl,pathdict,local=False)
    else:
        irx,iry,radio_components = pix_convert(dbl,pathdict)
    '''

    irx = float(dbl['ir_peak'][1:-1].split(',')[0]) / first_ir_scale_x
    iry = float(dbl['ir_peak'][1:-1].split(',')[1]) / first_ir_scale_y

    sub = subjects.find_one({'zooniverse_id':dbl['zooniverse_id']})
    r = requests.get(sub['location']['contours'])
    contours = r.json()
    radio_components = contours['contours']

    # Plot image
    
    zooniverse_id = cons['zid']
    answer = cons['answer']
    sub = subjects.find_one({'zooniverse_id':zooniverse_id})

    # Download contour data
    
    r = requests.get(sub['location']['contours'])
    contours = r.json()
    
    sf_x = 500./contours['width']
    sf_y = 500./contours['height']
    
    verts_all = []
    codes_all = []
    components = contours['contours']

    for comp in components:
    
        # Order of bounding box components is (xmax,ymax,xmin,ymin)
        comp_xmax,comp_ymax,comp_xmin,comp_ymin = comp[0]['bbox']
        
        # Only plot radio components identified by the users as the consensus;
        # check on the xmax value to make sure
        for v in answer.itervalues():
            if comp_xmax in v['xmax']:
    
                for idx,level in enumerate(comp):
                    verts = [((p['x'])*sf_x,(p['y']-1)*sf_y) for p in level['arr']]
                    
                    codes = np.ones(len(verts),int) * Path.LINETO
                    codes[0] = Path.MOVETO
                
                    verts_all.extend(verts)
                    codes_all.extend(codes)
    
    try:
        path = Path(verts_all, codes_all)
        patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)
    except AssertionError:
        print 'Users found no components for consensus match of %s' % zooniverse_id
    
    # Plot the infrared results
    
    fig = plt.figure(figno,(15,4))
    fig.clf()
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)
    
    colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
    colorarr = ['r','b','m','g','c','y','k'][::-1]
    
    if len(answer) > 0: # At least one galaxy was identified
        for idx,ans in enumerate(answer.itervalues()):

            if ans.has_key('peak_data'):

                # Plot the KDE map
                colormap = colormaparr.pop()
                ax3.imshow(np.rot90(ans['peak_data']['Z']), cmap=colormap,extent=[xmin, xmax, ymin, ymax])
        
                # Plot individual sources
                color = colorarr.pop()
                '''
                x_plot = [xt * 500./424 for xt in ans['ir_x'] if xt != -99.0]
                y_plot = [yt * 500./424 for yt in ans['ir_y'] if yt != -99.0]
                '''
                x_plot,y_plot = ans['ir_x'],ans['ir_y']
                ax3.scatter(x_plot, y_plot, c=color, marker='o', s=10, alpha=1./len(x_plot))

                ax4.plot([ans['ir_peak'][0]],[ans['ir_peak'][1]],color=color,marker='*',markersize=12)
    
            elif ans.has_key('ir'):
                color = colorarr.pop()
                x_plot,y_plot = ans['ir']
                ax3.plot([x_plot],[y_plot],color=color,marker='o',markersize=2)
                ax4.plot([x_plot],[y_plot],color=color,marker='*',markersize=12)
            else:
                ax4.text(550,idx*25,'#%i - no IR host' % idx,fontsize=11)

    
    ax3.set_xlim([0, 500])
    ax3.set_ylim([500, 0])
    ax3.set_title(zooniverse_id)
    ax3.set_aspect('equal')
    
    ax4.set_xlim([0, 500])
    ax4.set_ylim([500, 0])
    ax4.set_title('Consensus (%i/%i users)' % (cons['n_users'],cons['n_total']))
    
    ax4.set_aspect('equal')
    
    # Display IR and radio images
    
    url_standard = sub['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.set_title('WISE')

    url_radio = sub['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title(sub['metadata']['source'])
    ax2.get_yaxis().set_ticklabels([])

    ax3.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(answer) > 0:
        ax4.add_patch(patch_black)
        # Add centers of bounding boxes
        bbox1 = components[0][0]['bbox']
        bbox2 = components[1][0]['bbox']
        cx1 = np.median((bbox1[0],bbox1[2])) / first_ir_scale_x
        cx2 = np.median((bbox2[0],bbox2[2])) / first_ir_scale_x
        cy1 = np.median(((bbox1[1]/first_ir_scale_y),(bbox1[3]/first_ir_scale_y)))
        cy2 = np.median(((bbox2[1]/first_ir_scale_y),(bbox2[3]/first_ir_scale_y)))

        ax4.scatter(cx1,cy1, c='g', marker='s', s=15, alpha=1)
        ax4.scatter(cx2,cy2, c='g', marker='s', s=15, alpha=1)

        # Draw bounding box

        xmin1,xmin2 = bbox1[2] / first_ir_scale_x,bbox2[2] / first_ir_scale_x
        xmax1,xmax2 = bbox1[0] / first_ir_scale_x,bbox2[0] / first_ir_scale_x
        ymin1,ymin2 = bbox1[3]/first_ir_scale_y,bbox2[3]/first_ir_scale_y
        ymax1,ymax2 = bbox1[1]/first_ir_scale_y,bbox2[1]/first_ir_scale_y
        ax4.plot([xmin1,xmin1,xmax1,xmax1,xmin1],[ymin1,ymax1,ymax1,ymin1,ymin1],color='g')
        ax4.plot([xmin2,xmin2,xmax2,xmax2,xmin2],[ymin2,ymax2,ymax2,ymin2,ymin2],color='g')

        xc = ans['ir_peak'][0]
        yc = ans['ir_peak'][1]
        m1 = (cy1 - yc) / (cx1 - xc)
        b1 = yc - m1*xc
        m2 = (cy2 - yc) / (cx2 - xc)
        b2 = yc - m2*xc

        if cy1 < yc:
            xt,yt = -b1/m1,0
            xb,yb = (500.-b2)/m2,500
        else:
            xt,yt = -b2/m2,0
            xb,yb = (500.-b1)/m1,500

        ax4.plot([xt,xc],[yt,yc],color='orange',linestyle='--')
        ax4.plot([xb,xc],[yb,yc],color='orange',linestyle='--')

    alpha_deg = bending_angle(xc,yc,cx1,cy1,cx2,cy2) * 180/np.pi
    ax4.text(550,0,r'$\alpha$ = %.1f deg' % alpha_deg,fontsize=11)
    phi_deg = position_angle(xc,500-yc,cx1,500-cy1,cx2,500-cy2) * 180/np.pi
    ax4.arrow(xc,yc,0,-250,head_width=20, head_length=40, fc='grey', ec='grey',ls='dotted')
    ax4.text(550,50,r'$\phi$ = %.1f deg' % phi_deg,fontsize=11)

    # Draw the bisector vector

    print xt,yt
    print xb,yb
    print xc,yc
    center = np.array([xc,yc])
    r1 = np.array([xb,yb]) - center
    r2 = np.array([xt,yt]) - center
    r1norm = r1/np.hypot(r1[0],r1[1])
    r2norm = r2/np.hypot(r2[0],r2[1])
    r12 = r1+r2
    r12norm = r1norm+r2norm
    print r1norm
    print r2norm
    print r12norm
    print r1
    print r2
    print r1+r2

    ax4.arrow(xc,yc,r12[0],r12[1],head_width=20, head_length=40, fc='blue', ec='blue')



    ax4.yaxis.tick_right()
    
    ax1.get_xaxis().set_ticks([0,100,200,300,400])
    ax2.get_xaxis().set_ticks([0,100,200,300,400])
    ax3.get_xaxis().set_ticks([0,100,200,300,400])
    ax4.get_xaxis().set_ticks([0,100,200,300,400,500])

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if save_fig:
        fig.savefig('%s/bending_angles/plots/individual/ba_%s.pdf' % (rgz_dir,zooniverse_id))
        plt.close()
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    
    return None

def plot_some(n=100):

    doubles = get_doubles()
    pathdict = make_pathdict()
    somedoubles = random.sample(doubles,n)

    for dbl in doubles:
        plot_one_double(dbl,pathdict,save_fig=True)

    return None

def pix_dist(x1,y1,x2,y2):

    dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    return dist

def all_triples_pixradio(triples,pathdict):

    # Limit the triple set to those with an optical ID within 1 beam size of the center

    radiobeamsize = 5. # arcsec
    imagesize = 3. # arcmin
    imagescale = IMG_HEIGHT_NEW/imagesize / 60. # pixel / arcsec
    radio_tol = radiobeamsize * imagescale

    with open('%s/bending_angles/angles_triple_pixradio.csv' % rgz_dir,'w') as f:
        print >> f,'zooniverse_id,bending_angle,position_angle'
        for triple in triples:
            irx,iry,radio_components = pix_convert(triple,pathdict)
            xc,yc,radio_centroids = pix_radio(irx,iry,radio_components)
    
            maxdist = 0
            for centroid in radio_centroids:
                pd = pix_dist(xc,yc,centroid[0],centroid[1])
                maxdist = pd if pd > maxdist else maxdist
                if pd <= radio_tol:
                    middle_radio = centroid
                    radio_centroids.remove(middle_radio)
    
                    alpha = bending_angle(middle_radio[0],middle_radio[1],radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                    alpha_deg = alpha * 180./np.pi
                    phi = position_angle(middle_radio[0],middle_radio[1],radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                    phi_deg = phi * 180./np.pi
    
                    print >> f,'%s,%.3f' % (triple['zooniverse_id'],alpha_deg,phi_deg)
                    break
            '''
            else:
                print '%s; No radio match to optical ID; closest is %i pixels' % (triple['first_id'],maxdist)
            '''
    
    return None

def find_multipeaked_singles(subject,plot=False,verbose=True):

    # Find multi-peaked single component sources via binary kernels

    '''
    zooniverse_id = 'ARG000079s'
    zooniverse_id = 'ARG0002qyh'
    zooniverse_id = 'ARG0003l84'
    sub = subjects.find_one({'zooniverse_id':zooniverse_id})
    '''
    # Download contour data
    
    r = requests.get(subject['location']['contours'])
    contours = r.json()
    
    lobe = contours['contours'][0]
    # Order of bounding box components is (xmax,ymax,xmin,ymin)
    xmax,ymax,xmin,ymin = lobe[0]['bbox']
    xsize,ysize = 0.1,0.1
    X,Y = np.mgrid[xmin:xmax:xsize,ymin:ymax:ysize]

    parr = []
    valarr = []
    for cl in lobe:
        parr.extend([(p['x'],p['y']) for p in cl['arr']])
        valarr.extend(np.ones(len(cl['arr']))+cl['level'])

    points = np.array([(px,py) for px,py in parr])
    values = np.array(valarr)

    grid_z0 = griddata(points,values,(X,Y),method='nearest')
    grid_z1 = griddata(points,values,(X,Y),method='linear')
    grid_z2 = griddata(points,values,(X,Y),method='cubic')

    # Find the number of peaks
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array

    #neighborhood = generate_binary_structure(2,2)
    kernelsize = 50
    neighborhood = np.ones((kernelsize,kernelsize))

    '''
    Z = np.copy(grid_z2)
    Z[np.isnan(grid_z2)] = 1.
    '''
    Z = grid_z2

    local_max = maximum_filter(Z, footprint=neighborhood)==Z
    background = np.isnan(Z)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
    all_peaks = local_max - eroded_background

    # Check if peak is in the background

    detected_peaks = np.isfinite(Z) & all_peaks
    npeaks = detected_peaks.sum()

    xdp = X[detected_peaks]
    ydp = Y[detected_peaks]

    if verbose:
        print '%i peaks detected' % npeaks
        print xdp,ydp

    if plot:
        plt.subplot(231,aspect='equal')
        plt.plot(points[:,0], points[:,1], 'k.', ms=1)
        plt.title(subject['zooniverse_id'])

        plt.subplot(232)
        plt.imshow(grid_z0.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower',interpolation='none')
        plt.title('Nearest')

        plt.subplot(233)
        plt.imshow(grid_z1.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower',interpolation='none')
        plt.title('Linear')

        plt.subplot(234)
        plt.imshow(grid_z2.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower',interpolation='none')
        plt.title('Cubic')

        plt.subplot(235)
        plt.imshow(Z.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower')#,vmin=0.999,vmax=1.012)
        plt.title('Z')

        '''
        plt.subplot(235)
        plt.imshow(background.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower',interpolation='none')
        plt.title('Background')

        plt.subplot(235)
        plt.imshow(eroded_background.T, extent=(xmin,xmax,ymin,ymax), cmap = cm.cubehelix, origin='lower',interpolation='none')
        plt.title('Eroded background')
        '''

        plt.subplot(236,aspect='equal')
        plt.plot(points[:,0], points[:,1], 'k.', ms=1)
        plt.plot(xdp,ydp,'ro')
        plt.title('Detected peaks')
        plt.gcf().set_size_inches(18, 12)
        plt.show()

    return xdp,ydp

def centroid(arr):

    x = [l['x'] for l in arr]
    y = [l['y'] for l in arr]

    xmean = np.mean(x)
    ymean = np.mean(y)

    return xmean,ymean

def point_in_poly(x,y,poly):

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def make_polygon(arr):

    x = [l['x'] for l in arr]
    y = [l['y'] for l in arr]

    polygon = [(xx,yy) for xx,yy in zip(x,y)]

    return polygon
    
def mps_cc(subject,plot=True,verbose=True):

    # Find multi-peaked single component sources via contour counting

    '''
    zooniverse_id = 'ARG0002qyh'
    zooniverse_id = 'ARG000079s'
    zooniverse_id = 'ARG0003l84'
    subject = subjects.find_one({'zooniverse_id':zooniverse_id})
    '''

    # Download contour data
    
    r = requests.get(subject['location']['contours'])
    contours = r.json()
    
    lobe = contours['contours'][0]
    xmax,ymax,xmin,ymin = lobe[0]['bbox']

    parr = []
    valarr = []
    for cl in lobe:
        parr.extend([(p['x'],p['y']) for p in cl['arr']])
        valarr.extend(np.ones(len(cl['arr']))+cl['level'])

    points = np.array([(px,py) for px,py in parr])
    values = np.array(valarr)

    # Find levels with multiple contours
    # For each of those levels, check if next level up has geometric center within that contour
    # If no, then that level's geometric center is a local maximum
    # If yes, then move up one level and repeat

    k = [l['k'] for l in lobe]
    ck = Counter(k)

    mlarr = []
    for x,y in ck.iteritems():
        if y > 1:
            mlarr.append(x)

    if max(k) not in mlarr:
        mlarr.append(max(k))

    mlarr.sort()

    local_maxima = []
    for m in mlarr:
        levels = [l for l in lobe if l['k'] == m]

        # Is there a higher level?
        if m < max(k):

            upper_levels = [l for l in lobe if l['k'] == m+1]

            for level in levels:
                within = False
                for ul in upper_levels:

                    gc = centroid(ul['arr'])
                    polygon = make_polygon(level['arr'])
                    result = point_in_poly(gc[0],gc[1],polygon)
                    within += result

                if not within:
                    gc = centroid(level['arr'])
                    local_maxima.append((m,gc))
                    if verbose:
                        print 'Point in poly, m=%i, center=(%.1f,%.1f)' % (m,gc[0],gc[1])

        # If no higher level, centroids = local max
        else:
            for level in levels:
                gc = centroid(level['arr'])
                local_maxima.append((m,gc))
                if verbose:
                    print 'No higher levels, m=%i, center=(%.1f,%.1f)' % (m,gc[0],gc[1])

    # Plot locations of peaks

    npeaks = len(local_maxima)

    if plot:
        
        xc = [x[1][0] for x in local_maxima]
        yc = [x[1][1] for x in local_maxima]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        verts_all = []
        codes_all = []
        components = contours['contours']

        for comp in components:
        
            # Order of bounding box components is (xmax,ymax,xmin,ymin)
            comp_xmax,comp_ymax,comp_xmin,comp_ymin = comp[0]['bbox']
            
            # Only plot radio components identified by the users as the consensus;
            # check on the xmax value to make sure
        
            for idx,level in enumerate(comp):
                verts = [(p['x'],p['y']) for p in level['arr']]
                
                codes = np.ones(len(verts),int) * Path.LINETO
                codes[0] = Path.MOVETO
            
                verts_all.extend(verts)
                codes_all.extend(codes)
        
        try:
            path = Path(verts_all, codes_all)
            patch_black = patches.PathPatch(path, facecolor = 'none', edgecolor='black', lw=1)
        except AssertionError:
            print 'Users found no components for consensus match of %s' % zooniverse_id
        
        # Plot contours identified as the consensus
        ax.add_patch(patch_black)
        
        ax.plot(xc,yc,'r*',ms=10)
        ax.set_xlim(0,FIRST_FITS_WIDTH)
        ax.set_ylim(FIRST_FITS_HEIGHT,0)
        ax.set_aspect('equal')
        #ax.title(subject['zooniverse_id'])

        plt.show()

    return local_maxima

def batch_mps_cc():

    # Run all multi-peaked singles using contour-counting technique

    '''
    5013.01 seconds for 38750 images
    7.73 images per second
    '''

    tstart = time.time()
    mps = subjects.find({'state':'complete','metadata.contour_count':1},timeout=False)
    n = mps.count()

    with open('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,'w') as f:
        for subject in mps:
            local_maxima = mps_cc(subject,plot=False,verbose=False)
            if len(local_maxima) > 0:
                for idx,lm in enumerate(local_maxima):
                    print >> f,subject['zooniverse_id'],idx+1,len(local_maxima),lm[1][0],lm[1][1]
        
    mps.close()
    tend = time.time()

    print '%.2f seconds for %i images' % (tend - tstart,n)
    print '%.2f images per second' % (n/(tend - tstart))

    return None

def hist_mps_cc():

    # How many of the single-contour sources actually have more than one peak?

    data = ascii.read('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,delimiter=' ',data_start=1,header_start=0)

    c = Counter(data['zooniverse_id'])

    fig = plt.figure(1,(12,6))
    ax1 = fig.add_subplot(121)

    histML(c.values(), bins=range(10), ax=ax1, histtype='step', lw=2, alpha=1.0, color='#377eb8',log=True)
    ax1.set_xlim(1,10)
    ax1.set_xlabel(r'$N_{peaks}$',fontsize=18)
    ax1.set_ylabel('Count')
    ax1.set_title('RGZ 1-contour sources')

    ax2 = fig.add_subplot(122)
    histML(c.values(), bins=range(10), ax=ax2, histtype='step', lw=2, alpha=1.0, color='#e41a1c',cumulative=True,normed=True)
    ax2.set_xlabel(r'$N_{peaks}$',fontsize=18)
    ax1.set_title('RGZ 1-contour sources')
    ax2.set_ylabel('Cumulative fraction')

    fig.savefig('%s/bending_angles/plots/mps_cc.pdf' % rgz_dir)

    plt.show() 

    return None

def mps_cc_bendingangle():

    # What's the distribution of bending angles for the multipeaked double sources?

    '''
    Two sets:
        radio triples, no IR souces
        radio doubles, IR counterpart
    '''

    # Radio triples

    data = ascii.read('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,delimiter=' ',data_start=1,header_start=0)

    mps_triples = data[data['ntotal'] == 3]
    mps_triples_ids = set(mps_triples['zooniverse_id'])

    # Check if they had IR counterpart

    consensus_data = ascii.read('%s/csv/consensus.csv' % rgz_dir,delimiter=',',data_start=1,header_start=0)
    alpha_deg = []
    for id3 in mps_triples_ids:
        cons = consensus_data[consensus_data['zooniverse_id'] == id3]
        # Only keep ones without IR counterpart
        if cons['ir_peak'] == '(-99, -99)':
            t3 = mps_triples[mps_triples['zooniverse_id'] == id3]
            alpha = bending_angle(t3[0]['xc'],t3[0]['yc'],t3[1]['xc'],t3[1]['yc'],t3[2]['xc'],t3[2]['yc'])
            alpha_deg.append(alpha * 180./np.pi)
        else:
            print cons['ir_peak']

    # OK - currently claims there aren't any

    if len(alpha_deg) > 0:
        fig = plt.figure(2,(8,8))
        ax1 = fig.add_subplot(111)
        
        c1 = '#e41a1c'
        c2 = '#a6cee3'
        c3 = '#386cb0'
        
        histML(alpha_deg, bins=20, ax=ax1, histtype='stepfilled', alpha=1.0, color=c1, range=(0,180),label='Single-contour, triple-peaked radio sources w/o IR')
        ax1.set_xlim(0,180)
        ax1.vlines(x=np.median(alpha_deg),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color='k',linestyle='--')
        ax1.set_xlabel(r'bending angle [deg]',fontsize=24)
        ax1.set_ylabel('count',fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=20)
        
        ax1.legend(loc='upper left')
        
        fig.tight_layout()
        fig.savefig('%s/bending_angles/plots/mps_cc_bendingangle.pdf' % rgz_dir)


    return None

def batch_mps_kernel():

    # Run all multi-peaked singles using kernel technique

    tstart = time.time()
    mps = subjects.find({'state':'complete','metadata.contour_count':1})
    n = mps.count()

    with open('%s/bending_angles/multipeaked_singles.csv' % rgz_dir,'w') as f:
        for subject in mps:
            xdp,ydp = find_multipeaked_singles(subject,plot=False,verbose=False)
            if len(xdp) > 0 and len(ydp) > 0:
                for idx,(xsubpeak,ysubpeak) in enumerate(zip(xdp,ydp)):
                    print >> f,subject['zooniverse_id'],idx+1,len(xdp),xsubpeak,ysubpeak
        
    tend = time.time()

    print '%.2f seconds for %i images' % (tend - tstart,n)
    print '%.2f images per second' % (n/(tend - tstart))

    return None

def mps_bending_angle():

    df = pd.read_csv('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir)

    # Part 1: select double-peaked, single-contour sources with optical counterparts

    '''
    df2 = df[(df['nlobes'] == 2)]

    with open('%s/bending_angles/angles_multipeaked_singles.csv' % rgz_dir,'w') as f:

        print >> f,'zooniverse_id,bending_angle'

        # get optical counterpart for zooniverse_id (loop over eventually)

        df2_pair1 = df2[::2]
        df2_pair2 = df2[1::2]
        _zid = np.array(df2_pair1['zooniverse_id'])
        _x1 = np.array(df2_pair1['x'])
        _y1 = np.array(df2_pair1['y'])
        _x2 = np.array(df2_pair2['x'])
        _y2 = np.array(df2_pair2['y'])

        for zooniverse_id,cx1,cy1,cx2,cy2 in zip(_zid,_x1,_y1,_x2,_y2):

            c = consensus.checksum(zooniverse_id)
            if (len(c['answer']) == 1) and (c['n_users']/float(c['n_total']) >= 0.75):
                if c['answer'][c['answer'].keys()[0]].has_key('ir_peak'):
                    peak_x,peak_y = c['answer'][c['answer'].keys()[0]]['ir_peak']

                    ir_x = peak_x * first_ir_scale_x
                    ir_y = peak_y * first_ir_scale_y
                    alpha = bending_angle(ir_x,ir_y,cx1,cy1,cx2,cy2)
                    alpha_deg = alpha * 180./np.pi

                    print >> f,'%s,%.3f' % (zooniverse_id,alpha_deg)
                    #print ir_x,ir_y,cx1,cy1,cx2,cy2
                    #print alpha_deg

            else:
                print "Had more than 1 IR sources and/or less than 75 percent consensus for %s" % zooniverse_id

    '''
    # Part 2: select triple-peaked, single-contour sources with no optical counterparts

    df3 = df[(df['nlobes'] == 3)]
    imgcen_x,imgcen_y = FIRST_FITS_WIDTH,FIRST_FITS_HEIGHT 

    with open('%s/bending_angles/angles_multipeaked_singles_no_optical.csv' % rgz_dir,'w') as f:

        print >> f,'zooniverse_id,bending_angle,position_angle'

        # get optical counterpart for zooniverse_id (loop over eventually)

        df3_pair1 = df3[::3]
        df3_pair2 = df3[1::3]
        df3_pair3 = df3[2::3]
        _zid = np.array(df3_pair1['zooniverse_id'])
        _x1 = np.array(df3_pair1['x'])
        _y1 = np.array(df3_pair1['y'])
        _x2 = np.array(df3_pair2['x'])
        _y2 = np.array(df3_pair2['y'])
        _x3 = np.array(df3_pair3['x'])
        _y3 = np.array(df3_pair3['y'])


        for zooniverse_id,cx1,cy1,cx2,cy2,cx3,cy3 in zip(_zid,_x1,_y1,_x2,_y2,_x3,_y3):

            c = consensus.checksum(zooniverse_id)
            try:
                if (c['n_users']/float(c['n_total']) >= 0.75):

                    # Use the maximum of the three possible opening angles, since we don't know the "center" of the galaxy
                    
                    alpha1 = bending_angle(cx1,cy1,cx2,cy2,cx3,cy3)
                    alpha2 = bending_angle(cx2,cy2,cx1,cy1,cx3,cy3)
                    alpha3 = bending_angle(cx3,cy3,cx1,cy1,cx2,cy2)
                    alpha_max = max(alpha1,alpha2,alpha3)

                    # Compute position angle between angular bisector and north

                    if alpha_max == alpha1:
                        alpha_deg = alpha1 * 180/np.pi
                        phi_deg = position_angle(cx1,cy1,cx2,cy2,cx3,cy3) * 180/np.pi
                    elif alpha_max == alpha2:
                        alpha_deg = alpha2 * 180/np.pi
                        alpha2 = bending_angle(cx2,cy2,cx1,cy1,cx3,cy3) * 180/np.pi
                    elif alpha_max == alpha3:
                        alpha_deg = alpha3 * 180/np.pi
                        alpha2 = bending_angle(cx3,cy3,cx1,cy1,cx2,cy2) * 180/np.pi
                    else:
                        print 'Warning: no match between opening angle and max angle'
                    
                    print >> f,'%s,%.3f,%.3f' % (zooniverse_id,alpha_deg,phi_deg)
                else:
                    print "Had less than 75 percent consensus for %s" % zooniverse_id
            except TypeError:
                print "No 'answer' key for %s" % zooniverse_id,c

    return None

if __name__ == '__main__':

    mps_bending_angle()
