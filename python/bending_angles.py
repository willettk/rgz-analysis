# Measure bending angles for double radio sources in the RGZ data
#
# Kyle Willett, UMN

from __future__ import division

# Local RGZ modules

import rgz
import consensus
from load_contours import get_contours,make_pathdict

# Default packages

import json
import cStringIO
import urllib
import time
import random
import os
from ast import literal_eval
from collections import Counter

# Other packages

import pandas as pd

import numpy as np

from astropy.io import ascii,fits
from astropy import wcs

from PIL import Image

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.path import Path
import matplotlib.patches as patches

from scipy.interpolate import griddata
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

from astroML.plotting import hist as histML

# Local paths and files
rgz_dir = '/Users/willettk/Astronomy/Research/GalaxyZoo/rgz-analysis'
rgz_consensus_file = '%s/csv/consensus_rgz_first.csv' % rgz_dir

# Various image parameters

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

def get_doubles(consensus_level=0.50):

    # Find examples of RGZ subjects with exactly two radio components

    rgzconsensus = ascii.read(rgz_consensus_file,format='csv')

    dblidx = (rgzconsensus['n_radio'] == 2) & (rgzconsensus['consensus_level'] >= consensus_level)
    doubles = rgzconsensus[dblidx]

    return doubles

def get_triples(consensus_level=0.50):

    # Find examples of RGZ subjects with exactly three radio components

    rgzconsensus = ascii.read(rgz_consensus_file,format='csv')

    trpidx = (rgzconsensus['n_radio'] == 3) & (rgzconsensus['consensus_level'] >= consensus_level)
    triples = rgzconsensus[trpidx]

    return triples

def all_doubles_pixradio(doubles,pathdict):

    # Compute the coordinates of the optical ID and radio component centroids, bending angle, and position angle
    # for all consensus RGZ subjects with exactly two radio components

    with open('%s/bending_angles/angles_double_pixradio.csv' % rgz_dir,'w') as f:
        print >> f,'zooniverse_id,bending_angle,position_angle'
        for double in doubles:
            #irx,iry,radio_components = pix_convert(double,pathdict)
            xc,yc = literal_eval(double['ir_peak'])
            if xc is not None:
                subject = subjects.find_one({'zooniverse_id':double['zooniverse_id']})
                contours = get_contours(subject,pathdict)
                radio_components = contours['contours']

                radio_centroids = pix_radio(radio_components)
                alpha = bending_angle(xc,yc,radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                alpha_deg = alpha * 180./np.pi
                phi = position_angle(xc,yc,radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                phi_deg = phi * 180./np.pi
                if alpha is not None:
                    print >> f,'%s,%.3f,%.3f' % (double['zooniverse_id'],alpha_deg,phi_deg)

    return None

def dblid(doubles,zooniverse_id):

    # Retrieve subject for a single two-component radio source in the doubles list

    dbl = doubles[doubles['zooniverse_id'] == zooniverse_id][0]

    return dbl

def pix_convert(galaxy,pathdict,local=False):

    # Convert IR coordinates from RA/dec into pixel

    subject = subjects.find_one({'zooniverse_id':galaxy['zooniverse_id']})
    contours = get_contours(subject,pathdict)
    
    radio_components = contours['contours']

    try:
        assert len(radio_components) > 1, \
            'Radio data only has %i component for %s' % (len(radio_components),galaxy['zooniverse_id'])
    except AssertionError:
        return None,None,None

    # Keep everything in pixel coordinates. Reverse what's done in consensus.py; 
    # transform an RA/dec pair into the pixel x/y pair.

    # Convert the pixel coordinates into RA,dec using the WCS object from the header

    hdulist = fits.open(pathdict[galaxy['first_id']])

    w = wcs.WCS(hdulist[0].header)

    worldcrd = np.array([[galaxy['ra'],galaxy['dec']]],np.float_)
    pix = w.wcs_world2pix(worldcrd,0)

    irx,iry = pix[0]
    irx_first,iry_first = np.round(pix[0][0] / first_ir_scale_x), np.round(IMG_HEIGHT_NEW - pix[0][1] / first_ir_scale_y)

    return irx,iry,radio_components

def bending_angle(xc,yc,x1,y1,x2,y2):

    # Compute the bending angle (in radians) between three points in pixel space

    '''
    Points are:

        - xc,yc:  x,y center of IR counterpart
        - x1,y1:  x,y center of 1st radio lobe
        - x2,y2:  x,y center of 2nd radio lobe
    '''

    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
    center = np.array([xc,yc])

    r1diff = r1-center
    r2diff = r2-center

    r1len = np.hypot(r1diff[0],r1diff[1])
    r2len = np.hypot(r2diff[0],r2diff[1])

    alpha = np.arccos(np.dot(r1diff,r2diff) / (r1len*r2len))

    return alpha

def bending_angle_sdss(zid,x1,y1,x2,y2):

    # Compute the bending angle (in radians) between three points in pixel space

    '''

        - zid
        - x1,y1:  x,y center of 1st radio lobe
        - x2,y2:  x,y center of 2nd radio lobe
    '''

    # I'd love to do this purely in RA/dec, converting all positions in Astropy, but functionality doesn't seem to be there.

    # Convert SDSS optical position into radio-frame pixel coordinates
    hdulist = fits.open(pathdict[galaxy['first_id']])

    w = wcs.WCS(hdulist[0].header)

    worldcrd = np.array([[galaxy['ra'],galaxy['dec']]],np.float_)
    pix = w.wcs_world2pix(worldcrd,0)

    xc,yc = np.round(pix[0][0] / first_ir_scale_x), np.round(IMG_HEIGHT_NEW - pix[0][1] / first_ir_scale_y)


    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
    center = np.array([xc,yc])

    r1diff = r1-center
    r2diff = r2-center

    r1len = np.hypot(r1diff[0],r1diff[1])
    r2len = np.hypot(r2diff[0],r2diff[1])

    alpha = np.arccos(np.dot(r1diff,r2diff) / (r1len*r2len))

    return alpha

def position_angle(xc,yc,x1,y1,x2,y2):

    # Compute the position angle (in radians, with respect to north) between three points in pixel space

    '''
    Points are:

        - xc,yc:    x,y center of bending angle
        - x1,y1:  x,y center of 1st component
        - x2,y2:  x,y center of 2nd component
    '''

    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
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

def pix_radio(radio_components):

    # From list of bounding boxes in radio pixel coordinates,
    # return the centroids of the boxes in IR pixel coordinates

    radio_centroids = []
    for comp in radio_components:
        bbox = comp[0]['bbox']

        cxu = np.median((bbox[0],bbox[2]))
        cyu = np.median((bbox[1],bbox[3]))

        cx,cy = cxu/first_ir_scale_x,cyu/first_ir_scale_y
        radio_centroids.append((cx,cy))

    return radio_centroids

def bbox_radio_to_ir(bbox):

    # Convert the bbox in RGZ subject from radio to infrared pixel scale

    bbox_ir = [bbox[0]/first_ir_scale_x,bbox[1]/first_ir_scale_x,bbox[2]/first_ir_scale_x,bbox[3]/first_ir_scale_x]

    return bbox_ir

def load_angles(filename):

    # Load the CSV file of the computed bending angles for multi-peaked or multi-lobed sources
     
    with open('%s/bending_angles/%s.csv' % (rgz_dir,filename),'r') as f:
        angstr = f.readlines()

    ba = [float(x.split(',')[1]) for x in angstr[1:]]
    #pa = [float(x.split(',')[2]) for x in angstr[1:]]

    return ba#,pa

def plothist(savefig=False):

    # Plot distribution of the bending angles for RGZ sources

    '''
    angles_double_pixradio = load_angles('angles_double_pixradio')
    angles_triple_pixradio = load_angles('angles_triple_pixradio')
    angles_double_mps = load_angles('angles_multipeaked_singles')
    angles_triple_mps = load_angles('angles_multipeaked_singles_no_optical')
    '''

    data = ascii.read('{:}/csv/static_catalog3.csv'.format(rgz_dir),delimiter=' ')

    angles_radiodouble = data[data['angle_type'] == 'double_pixradio']['bending_angle']
    angles_mps = data[data['angle_type'] == 'multipeaked_singles']['bending_angle']

    # Set up figure

    fig = plt.figure(2,(15,8))

    c1 = '#377eb8'
    c2 = '#e41a1c'
    c3 = '#4daf4a'
    c4 = '#984ea3'

    # Panel 1 - histogram
    ax1 = fig.add_subplot(121)
    
    histML(angles_radiodouble, bins=15, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c1, range=(0,90),label='double lobed, multi-contour')
    histML(angles_mps, bins=15, ax=ax1, histtype='step', lw=3, alpha=1.0, color=c2, range=(0,90),label='double-peaked, single-contour')

    ax1.set_xlim(0,90)
    ax1.vlines(x=np.median(angles_radiodouble),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c1,linestyle='--')
    ax1.vlines(x=np.median(angles_mps),ymin=ax1.get_ylim()[0],ymax = ax1.get_ylim()[1],color=c2,linestyle='--')
    ax1.set_xlabel(r'bending angle [deg]',fontsize=24)
    ax1.set_ylabel('count',fontsize=20)

    plt.tick_params(axis='both', which='major', labelsize=20)
    # Panel 2 - cumulative

    ax2 = fig.add_subplot(122)
    
    histML(angles_radiodouble, bins=15, ax=ax2, histtype='step', lw=3, alpha=1.0, color=c1, range=(0,90),label='double lobed, multi-contour',cumulative=True)
    histML(angles_mps, bins=15, ax=ax2, histtype='step', lw=3, alpha=1.0, color=c2, range=(0,90),label='double-peaked, single-contour',cumulative=True)

    ax2.set_xlim(0,90)
    ax2.vlines(x=np.median(angles_radiodouble),ymin=ax2.get_ylim()[0],ymax = ax2.get_ylim()[1],color=c1,linestyle='--')
    ax2.vlines(x=np.median(angles_mps),ymin=ax2.get_ylim()[0],ymax = ax2.get_ylim()[1],color=c2,linestyle='--')
    ax2.set_xlabel(r'bending angle [deg]',fontsize=24)
    ax2.set_ylabel('count',fontsize=20)
    ax2.legend(loc='upper left')

    plt.tick_params(axis='both', which='major', labelsize=20)
    # Finish adjusting plot parameters

    fig.tight_layout()

    if savefig:
        fig.savefig('%s/bending_angles/plots/bending_angles_hist.pdf' % rgz_dir)
    else:
        plt.show()

    return None

def plot_one_double(zooniverse_id,pathdict,figno=1,save_fig=False,anglepath='',dbltype='radio'):

    # Make a four-panel plot of the consensus identification with marked bending angle and position angle for a double source

    cons = consensus.checksum(zooniverse_id)

    subject = subjects.find_one({'zooniverse_id':zooniverse_id})
    contours = get_contours(subject,pathdict)
    radio_components = contours['contours']

    # Plot image
    
    answer = cons['answer']

    # Download contour data
    
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
    
    url_standard = subject['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.set_title('WISE')

    url_radio = subject['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title(subject['metadata']['source'])
    ax2.get_yaxis().set_ticklabels([])

    ax3.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(answer) > 0:
        ax4.add_patch(patch_black)
        radio_centers = []
        for component in components:
            bbox = component[0]['bbox']

            # Draw centers of bounding boxes
            xradiocen = np.median((bbox[0],bbox[2])) / first_ir_scale_x
            yradiocen = np.median((bbox[1],bbox[3])) / first_ir_scale_y
            radio_centers.append((xradiocen,yradiocen))

            ax4.scatter(xradiocen,yradiocen, c='g', marker='s', s=15, alpha=1)

            # Draw edge of bounding box

            xradiomin = bbox[2] / first_ir_scale_x
            xradiomax = bbox[0] / first_ir_scale_x
            yradiomin = bbox[3] / first_ir_scale_y
            yradiomax = bbox[1] / first_ir_scale_y
            ax4.plot([xradiomin,xradiomin,xradiomax,xradiomax,xradiomin],[yradiomin,yradiomax,yradiomax,yradiomin,yradiomin],color='g')

        for ans in answer:
            if answer[ans].has_key('ir_peak'):
                # Optical counterpart position
                xc,yc = answer[ans]['ir_peak']     

                # Position of radio sources for multi-peaked, single-component subjects
                if dbltype == "mps":
                    local_maxima = mps_cc(subject,pathdict,plot=False,verbose=False)
                    suffix = '' if len(local_maxima) == 1 else 's'
                    assert len(local_maxima) >= 2, \
                        "%i peak%s in first radio component of %s; must have exactly 2 peaks to plot bending angle using mps method." % (len(local_maxima),suffix,zooniverse_id)
                    x1 = local_maxima[0][1][0] / first_ir_scale_x
                    y1 = local_maxima[0][1][1] / first_ir_scale_y
                    x2 = local_maxima[1][1][0] / first_ir_scale_x
                    y2 = local_maxima[1][1][1] / first_ir_scale_y
                    ax4.scatter(x1,y1, color='darkorange', marker='s', s=15, alpha=1)
                    ax4.scatter(x2,y2, color='darkorange', marker='s', s=15, alpha=1)

                # Position of radio sources for double-lobed, two-component subjects
                elif len(radio_centers) == 2:
                    x1,y1 = radio_centers[0]
                    x2,y2 = radio_centers[1]
                else:
                    raise ValueError("Centers of radio boxes not defined.")

                m1 = (y1 - yc) / (x1 - xc)
                b1 = yc - m1*xc
                m2 = (y2 - yc) / (x2 - xc)
                b2 = yc - m2*xc

                xedge1 = 0 if x1 < xc else 500
                yedge1 = y1 - (x1-xedge1)*(yc-y1)/(xc-x1)

                xedge2 = 0 if x2 < xc else 500
                yedge2 = y2 - (x2-xedge2)*(yc-y2)/(xc-x2)

                # Draw and annotate the the bending angle

                ax4.plot([xedge1,xc],[yedge1,yc],color='orange',linestyle='--')
                ax4.plot([xedge2,xc],[yedge2,yc],color='orange',linestyle='--')
                alpha_deg = bending_angle(xc,yc,x1,y1,x2,y2) * 180/np.pi
                ax4.text(550,0,r'$\alpha$ = %.1f deg' % alpha_deg,fontsize=11)

                # Draw vector pointing north

                # Draw the bisector vector
                '''
                yd = y_bisect(xc,yc,xedge1,yedge1,xedge2,yedge2)
                ax4.arrow(xc,yc,-xc,yd-yc,head_width=20, head_length=40, fc='blue', ec='blue')
                '''

                # Compute the position angle with respect to north
                phi_deg = position_angle(xc,500-yc,x1,500-y1,x2,500-y2) * 180/np.pi
                ax4.text(550,50,r'$\phi$ = %.1f deg' % phi_deg,fontsize=11)
                ax4.arrow(xc,yc,0,-yc,head_width=20, head_length=40, fc='grey', ec='grey',ls='dotted')

            else:
                print "No peak for %s" % zooniverse_id


    ax4.yaxis.tick_right()
    
    ax1.get_xaxis().set_ticks([0,100,200,300,400])
    ax2.get_xaxis().set_ticks([0,100,200,300,400])
    ax3.get_xaxis().set_ticks([0,100,200,300,400])
    ax4.get_xaxis().set_ticks([0,100,200,300,400,500])

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if save_fig:
        fig.savefig('%s/bending_angles/plots/individual/%sba_%s.pdf' % (rgz_dir,anglepath,zooniverse_id))
        plt.close()
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    
    return None

def plot_one_triple(zooniverse_id,pathdict,figno=1,save_fig=False,anglepath=''):

    # Make a four-panel plot of the consensus identification with marked bending angle and position angle for a triple source

    cons = consensus.checksum(zooniverse_id)

    subject = subjects.find_one({'zooniverse_id':zooniverse_id})
    contours = get_contours(subject,pathdict)
    radio_components = contours['contours']

    # Plot image
    
    answer = cons['answer']

    # Download contour data
    
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
    
    url_standard = subject['location']['standard']
    im_standard = Image.open(cStringIO.StringIO(urllib.urlopen(url_standard).read()))
    ax1 = fig.add_subplot(141)
    ax1.imshow(im_standard,origin='upper')
    ax1.set_title('WISE')

    url_radio = subject['location']['radio']
    im_radio = Image.open(cStringIO.StringIO(urllib.urlopen(url_radio).read()))
    ax2 = fig.add_subplot(142)
    ax2.imshow(im_radio,origin='upper')
    ax2.set_title(subject['metadata']['source'])
    ax2.get_yaxis().set_ticklabels([])

    ax3.get_yaxis().set_ticklabels([])

    # Plot contours identified as the consensus
    if len(answer) > 0:
        ax4.add_patch(patch_black)
        # Add centers of bounding boxes
        for comp in components:
            bbox_radio = comp[0]['bbox']
            bbox_ir = bbox_radio_to_ir(bbox_radio)
            xrad = np.median((bbox_ir[0],bbox_ir[2]))
            yrad = np.median((bbox_ir[1],bbox_ir[3]))
            ax4.scatter(xrad,yrad, c='g', marker='s', s=15, alpha=1)

            dbx = [bbox_ir[i] for i in (2,2,0,0,2)]
            dby = [bbox_ir[i] for i in (3,1,1,3,3)]
            ax4.plot(dbx,dby,color='g')

        radiobeamsize = 5. # arcsec
        imagesize = 3. # arcmin
        imagescale = IMG_HEIGHT_NEW/imagesize / 60. # pixel / arcsec
        radio_tol = radiobeamsize * imagescale

        for ans in answer:
            if answer[ans].has_key('ir_peak'):
                # Optical counterpart position
                xc,yc = answer[ans]['ir_peak']     

                # Measure all positions in radio pixel coordinates
                radio_centroids = pix_radio(components)
    
                maxdist = 0
                for centroid in radio_centroids:
                    d = pix_dist(xc,yc,centroid[0],centroid[1])
                    maxdist = d if d > maxdist else maxdist
                    if d <= radio_tol:
                        middle_radio = centroid
                        radio_centroids.remove(middle_radio)
    
                        x1 = radio_centroids[0][0]
                        y1 = radio_centroids[0][1]
                        x2 = radio_centroids[1][0]
                        y2 = radio_centroids[1][1]

                if len(radio_centroids) == 2:

                    m1 = (y1 - yc) / (x1 - xc)
                    b1 = yc - m1*xc
                    m2 = (y2 - yc) / (x2 - xc)
                    b2 = yc - m2*xc

                    xedge1 = 0 if x1 < xc else 500
                    yedge1 = y1 - (x1-xedge1)*(yc-y1)/(xc-x1)

                    xedge2 = 0 if x2 < xc else 500
                    yedge2 = y2 - (x2-xedge2)*(yc-y2)/(xc-x2)

                    # Draw and annotate the the bending angle

                    ax4.plot([xedge1,xc],[yedge1,yc],color='orange',linestyle='--')
                    ax4.plot([xedge2,xc],[yedge2,yc],color='orange',linestyle='--')
                    alpha_deg = bending_angle(xc,yc,x1,y1,x2,y2) * 180/np.pi
                    ax4.text(550,0,r'$\alpha$ = %.1f deg' % alpha_deg,fontsize=11)

                else:
                    print "\tDidn't find match to optical ID for triple radio source %s" % zooniverse_id

            else:
                print "\tNo IR peak for %s" % zooniverse_id

    ax4.yaxis.tick_right()
    
    ax1.get_xaxis().set_ticks([0,100,200,300,400])
    ax2.get_xaxis().set_ticks([0,100,200,300,400])
    ax3.get_xaxis().set_ticks([0,100,200,300,400])
    ax4.get_xaxis().set_ticks([0,100,200,300,400,500])

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if save_fig:
        fig.savefig('{0:}/bending_angles/plots/individual/triples/{1:}ba_{2:}.pdf'.format(rgz_dir,anglepath,zooniverse_id))
        plt.close()
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    
    return None

def y_bisect(xc,yc,xt,yt,xb,yb):

    # Finds the point yd such that the vector (xc,yc) -> (0,yd)
    # bisects the angle formed by the vectors (xb,yb) -> (xc,yc)
    #                                     and (xt,yt) -> (xc,yc)

    bc_length = np.hypot(xb - xc,yb - yc)
    tc_length = np.hypot(xt - xc,yt - yc)
    numerator = ((xb - xc)*xc + (yb - yc)*yc)/bc_length - (xc*(xt - xc) + yc*(yt - yc))/tc_length
    denominator = (yb - yc)/bc_length - (yt - yc)/tc_length

    return numerator/denominator

def plot_some(n,random_selection=False):

    # Plot a random selection of double and triple sources

    pathdict = make_pathdict()

    # Doubles

    doubles = get_doubles()
    somedoubles = random.sample(doubles,n) if random_selection else doubles[:n]

    for dbl in somedoubles:
        plot_one_double(dbl['zooniverse_id'],pathdict,save_fig=True)

    # Triples
    sometriples = get_triples()
    for triple in sometriples:
        plot_one_triple(triple['zooniverse_id'],pathdict,save_fig=True)

    return None

def pix_dist(x1,y1,x2,y2):

    # Find the distance between two sets of Cartesian points via the Pythagorean theorem

    dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    return dist

def all_triples_pixradio(triples,pathdict):

    # Compute the bending angle for RGZ subjects with three radio components and an optical ID within 1 beam size of the center component

    radiobeamsize = 5. # arcsec
    imagesize = 3. # arcmin
    imagescale = IMG_HEIGHT_NEW/imagesize / 60. # pixel / arcsec
    radio_tol = radiobeamsize * imagescale

    with open('%s/bending_angles/angles_triple_pixradio.csv' % rgz_dir,'w') as f:
        print >> f,'zooniverse_id,bending_angle,position_angle'
        for triple in triples:
            irx,iry = literal_eval(triple['ir_peak'])

            subject = subjects.find_one({'zooniverse_id':triple['zooniverse_id']})
            contours = get_contours(subject,pathdict)
            radio_components = contours['contours']

            # Measure all positions in radio pixel coordinates
            radio_centroids = pix_radio(radio_components)
    
            maxdist = 0
            for centroid in radio_centroids:
                d = pix_dist(irx,iry,centroid[0],centroid[1])
                maxdist = d if d > maxdist else maxdist
                if d <= radio_tol:
                    middle_radio = centroid
                    radio_centroids.remove(middle_radio)
    
                    alpha = bending_angle(middle_radio[0],middle_radio[1],radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                    alpha_deg = alpha * 180./np.pi
                    phi = position_angle(middle_radio[0],middle_radio[1],radio_centroids[0][0],radio_centroids[0][1],radio_centroids[1][0],radio_centroids[1][1])
                    phi_deg = phi * 180./np.pi
    
                    print >> f,'%s,%.3f,%.3f' % (triple['zooniverse_id'],alpha_deg,phi_deg)
                    break
                else:
                    "Couldn't match the optical ID within 1 beam size of center for %s" % triple['zooniverse_id']

    return None

def find_multipeaked_singles(subject,plot=False,verbose=True):

    # Deprecated in favor of mps_cc
    # Find multi-peaked single component sources via binary kernels. 

    # Download contour data
    
    contours = get_contours(subject,pathdict)
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

    # Find the centroid of a polygon defined by a list of (x,y) points

    x = [l['x'] for l in arr]
    y = [l['y'] for l in arr]

    xmean = np.mean(x)
    ymean = np.mean(y)

    return xmean,ymean

def point_in_poly(x,y,poly):

    # Determine whether a given point (x,y) is within a convex polygon defined by an array of points

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

    # Create a list of x,y pairs out of an array to draw a polygon

    x = [l['x'] for l in arr]
    y = [l['y'] for l in arr]

    polygon = [(xx,yy) for xx,yy in zip(x,y)]

    return polygon
    
def mps_cc(subject,pathdict,plot=True,verbose=True):

    # Find location of peaks within a single-component radio source via contour counting

    contours = get_contours(subject,pathdict)
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

    # Find location of peaks within all single-component radio sources via contour counting
    
    '''
    Time estimate:

    Contour data retrieved over network:
        5013.01 seconds (~83 minutes) for 38,750 images
        7.73 images per second

    Contour data stored locally:
        1559.67 seconds (~26 minutes) for 46,068 images
        29.54 images per second

    '''
    
    # Note - only enabled for FIRST now.
    
    tstart = time.time()
    mps = subjects.find({'state':'complete','metadata.contour_count':1,'metadata.survey':'first'},timeout=False)
    n = mps.count()

    with open('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,'w') as f:
        print >> f,"zooniverse_id,nlobe,ntotal,xc,yc"
        idx_s = 0
        for subject in mps:
            try:
                local_maxima = mps_cc(subject,pathdict,plot=False,verbose=False)
                if len(local_maxima) > 1:
                    for idx,lm in enumerate(local_maxima):
                            print >> f,"{0:},{1:d},{2:d},{3:.4f},{4:.4f}".format(subject['zooniverse_id'],idx+1,len(local_maxima),lm[1][0],lm[1][1])
            except ValueError:
                print "Error retrieving JSON object for {0:}".format(subject['zooniverse_id'])
            idx_s += 1
        if ~idx_s % 100:
            print "%i completed" % idx_s
        
    mps.close()
    tend = time.time()

    print '%.2f minutes for %i images' % ((tend - tstart)/60.,n)
    print '%.2f images per second' % (n/(tend - tstart))

    return None

def hist_mps_cc():

    # Plot the distribution of the number of peaks in single-component radio subjects

    data = ascii.read('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,delimiter=' ',data_start=1,header_start=0)

    c = Counter(data['zooniverse_id'])

    # Doesn't include npeaks = 1, so calculate that separately

    ntotal = subjects.find({'state':'complete','metadata.contour_count':1,'metadata.survey':'first'}).count()
    runningsum = 0
    for v in c.itervalues():
        runningsum += v
    c[1] = ntotal - runningsum

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

def batch_mps_kernel():

    # Deprecated in favor of batch_mps_cc
    # Find location of peaks within all single-component radio sources via binary kernels

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

def mps_bending_angle(consensus_level = 0.50):

    # Compute the bending and position angles for double-peaked, single-contour sources with optical counterparts

    tstart = time.time()
    
    # Load data
    df = pd.read_csv('%s/bending_angles/multipeaked_singles_cc.csv' % rgz_dir,delimiter=',')

    df2 = df[(df['ntotal'] == 2)]

    with open('%s/bending_angles/angles_multipeaked_singles.csv' % rgz_dir,'w') as f:

        print >> f,'zooniverse_id,bending_angle,position_angle'

        # get optical counterpart for zooniverse_id (loop over eventually)

        df2_pair1 = df2[::2]
        df2_pair2 = df2[1::2]
        _zid = np.array(df2_pair1['zooniverse_id'])
        _x1 = np.array(df2_pair1['xc'])
        _y1 = np.array(df2_pair1['yc'])
        _x2 = np.array(df2_pair2['xc'])
        _y2 = np.array(df2_pair2['yc'])

        for zooniverse_id,x1,y1,x2,y2 in zip(_zid,_x1,_y1,_x2,_y2):

            c = consensus.checksum(zooniverse_id)
            try:
                if (len(c['answer']) == 1) and (c['n_users']/float(c['n_total']) >= consensus_level):
                    if c['answer'][c['answer'].keys()[0]].has_key('ir_peak'):
                        peak_x,peak_y = c['answer'][c['answer'].keys()[0]]['ir_peak']

                        ir_x = peak_x * first_ir_scale_x
                        ir_y = peak_y * first_ir_scale_y
                        alpha = bending_angle(ir_x,ir_y,x1,y1,x2,y2)
                        alpha_deg = alpha * 180./np.pi
                        phi = position_angle(ir_x,ir_y,x1,y1,x2,y2)
                        phi_deg = phi * 180./np.pi

                        print >> f,'{0:s},{1:.4f},{2:.4f}'.format(zooniverse_id,alpha_deg,phi_deg)

                else:
                    print "Had more than 1 IR sources and/or less than {0:.2f} percent consensus for {1:s}".format(consensus_level,zooniverse_id)
            except TypeError:
                print "No 'answer' key for %s" % zooniverse_id

    # Timing the process
    tend = time.time()

    n = len(df2)/2
    print '%.2f minutes for %i subjects' % ((tend - tstart)/60.,n)
    print '%.2f subjects per second' % (n/(tend - tstart))


    return None

if __name__ == '__main__':

    # If run from command line, computes bending angles for all double and multi-peaked single component RGZ sources

    pathdict = make_pathdict()

    doubles = get_doubles()
    all_doubles_pixradio(doubles,pathdict)

    batch_mps_cc()
    mps_bending_angle()

    '''
    triples = get_triples()
    all_triples_pixradio(triples,pathdict)
    '''
