from consensus import *
import pandas as pd
import csv
import operator
from ast import literal_eval
import numpy as np
import matplotlib.pyplot as plt

'''

Some very quick and dirty analysis on the outputs of the RGZ consensus code and testing the effect
of various weighting schemes.  --- Kyle Willett (UMN), 25 May 2016

To do: check to see what difference changing the weights has had. I tried generating three data sets:

    1. Consensus by plurality vote (default method).
    2. Consensus by plurality, but people who saw at least 5 gold standard subjects and agreed with the science team >50% had their weights doubled.
    3. Consensus by plurality, but people who saw at least 5 gold standard subjects and agreed with the science team >50% had their weights increased by a factor of ten.

N_0  = 110277
N_2  = 104232
N_10 = 103725

Consensus numbers are slightly different, dropping when some users are highly upweighted. I guess that would occur if users who said "No Radio/No IR" were highly upweighted and thus the subject would be dropped. Born out by the screen messages on tabernacle.

What did the initial results show?
    
    1. Total consensus votes per subject goes up, as expected; scales roughly as factor of x2 or x10 for the weighted data.
    2. Consensus fraction goes up slightly for increasing weights. At the 75% level we've been using, the fraction of galaxies increases from ~60% to 70% of the sample. So they're more confident.
    3. The number of radio components per subjects doesn't appreciably change. Weighted votes are very slightly more likely to find more radio components in a particular image.
'''

def load_data():

    # Load the data sets for FIRST with different weighting schemes

    
    w01 = pd.read_csv('{0}/csv/weights/noweights.csv'.format(rgz_path))
    w02 = pd.read_csv('{0}/csv/weights/doubleweight.csv'.format(rgz_path))
    w10 = pd.read_csv('{0}/csv/weights/tenweight.csv'.format(rgz_path))
    
    return w01,w02,w10

def plot_hist():

    # Make a cumulative histogram for each of the various weighting schemes of:
    #   the total number of weighted votes being counted
    #   the consensus level in each of source
    #   the number of radio components in each source

    
    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,6))
    
    n_bins = 50
    
    # Axis 1
    
    ax1.hist(w01['n_total'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,200),label='w=1')
    ax1.hist(w02['n_total'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,200),label='w=2')
    ax1.hist(w10['n_total'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,200),label='w=10')
    
    ax1.set_xlabel('total votes per subject',fontsize=16)
    ax1.set_ylabel('cumulative fraction',fontsize=16)
    
    ax1.grid(True)
    ax1.set_ylim(0, 1.05)
    ax1.legend(loc='lower right',fontsize=12)
    
    n_bins = 25
    
    # Axis 2
    
    ax2.hist(w01['consensus_level'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,1),label='w=1')
    ax2.hist(w02['consensus_level'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,1),label='w=2')
    ax2.hist(w10['consensus_level'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,1),label='w=10')

    ax2.set_xlabel('consensus level',fontsize=16)
    ax2.set_ylabel('cumulative fraction',fontsize=16)
    
    ax2.grid(True)
    ax2.set_ylim(0, 1.05)
    ax2.legend(loc='lower right',fontsize=12)
    
    # Axis 3
    
    h01,b01,p01 = ax3.hist(w01['n_radio'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,10),label='w=1')
    h01,b02,p02 = ax3.hist(w02['n_radio'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,10),label='w=2')
    h10,b10,p10 = ax3.hist(w10['n_radio'], n_bins, normed=1, histtype='step', lw=2, cumulative=True, range=(0,10),label='w=10')
    
    ax3.set_xlabel('number of radio components',fontsize=16)
    ax3.set_ylabel('cumulative fraction',fontsize=16)
    
    ax3.grid(True)
    ax3.set_ylim(0, 1.05)
    ax3.legend(loc='lower right',fontsize=12)
    
    fig.tight_layout()
    plt.show()

    return None
    
def checksum_val(g):

    # Compute the unique checksum (product value of xmax) that identifies a unique set of radio components

    bbox_list = literal_eval(g)
    
    xmaxlist = [x[1] for x in bbox_list]
    product = reduce(operator.mul, xmaxlist, 1)

    return round(product,3)

def add_checksum():

    # Write checksum to its own column in the file of weights

    wf = open('{0}/csv/weights/tenweight.csv'.format(rgz_path),'w')
    writer = csv.writer(wf)
    
    with open('{0}/csv/tenweight/consensus_rgz_first.csv'.format(rgz_path)) as f:
        reader = csv.reader(f)
        header = reader.next()
        header.append('checksum')
        writer.writerow(header)
        for row in reader:
            checksum = checksum_val(row[7])
            newrow = row
            newrow.append(checksum)
            writer.writerow(newrow)
    
    wf.close()

def find_differences():

    # Find the number of galaxies that actually changed classifications with different weighting schemes

    w01,w02,w10 = load_data()
    grouped = w01.groupby('zooniverse_id')
    
    wf02 = open('{0}/csv/weights/single_double.csv'.format(rgz_path),'wb')
    writer02 = csv.writer(wf02)
    wf10 = open('{0}/csv/weights/single_tenfold.csv'.format(rgz_path),'wb')
    writer10 = csv.writer(wf10)
    
    for i,(name,group) in enumerate(grouped):
        n01 = len(group)
        n02 =(w02['zooniverse_id'] == name).sum() 
        n10 =(w10['zooniverse_id'] == name).sum() 
        if n02 > 0 and n02 != n01:
            writer02.writerow([name,n01,n02])
        if n10 > 0 and n10 != n01:
            writer10.writerow([name,n01,n10])
    
        # Roughly 10 seconds per
        if not i % 1000:
            print i

def plot_weighted_consensus(zid,savefig=False):

    # Plot 4-panel image of IR, radio, KDE estimate, and consensus

    consensus01 = checksum(zid)
    consensus02 = checksum(zid,weights=2)
    consensus10 = checksum(zid,weights=10)

    consenses = (consensus01,consensus02,consensus10)

    fig,axarr = plt.subplots(3,4,figsize=(15,13))

    for cons, axrow, w in zip(consenses,axarr,(1,2,10)):

        ax1,ax2,ax3,ax4 = axrow.ravel()
    
        zid = cons['zid']
        answer = cons['answer']
        sub = subjects.find_one({'zooniverse_id':zid})
        survey = sub['metadata']['survey']

        # Get contour data
        contours = get_contours(sub,pathdict)
        
        # Key bit that sets the difference between surveys.
        #   contours['width'] = img_params[survey]['FITS_WIDTH']
        #   contours['height'] = img_params[survey]['FITS_HEIGHT']

        sf_x = img_params[survey]['IMG_WIDTH_NEW'] * 1./contours['width']
        sf_y = img_params[survey]['IMG_HEIGHT_NEW'] * 1./contours['height']
        
        verts_all = []
        codes_all = []
        components = contours['contours']

        for comp in components:
        
            # Order of bounding box components is (xmax,ymax,xmin,ymin)
            comp_xmax,comp_ymax,comp_xmin,comp_ymin = comp[0]['bbox']
            
            # Only plot radio components identified by the users as most common answer
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
            print 'Users found no components for consensus match of {0}'.format(zid)
        
        # Plot the infrared results
        
        colormaparr = [cm.hot_r,cm.Blues,cm.RdPu,cm.Greens,cm.PuBu,cm.YlGn,cm.Greys][::-1]
        colorarr = ['r','b','m','g','c','y','k'][::-1]
        
        # If, in the rare case, that the consensus has more unique sources than the number of colors:
        if len(answer) > len(colorarr):
            colorarr *= int(len(answer)/len(colorarr))+1
            colormaparr *= int(len(answer)/len(colorarr))+1
        
        if len(answer) > 0: # At least one galaxy was identified
            for idx,ans in enumerate(answer.itervalues()):

                if ans.has_key('peak_data'):

                    xmin = 1.
                    xmax = img_params[survey]['IMG_HEIGHT_NEW']
                    ymin = 1.
                    ymax = img_params[survey]['IMG_WIDTH_NEW']

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
                    ax4.text(img_params[survey]['IMG_WIDTH_NEW']+50,idx*25,'#{0:d} - no IR host'.format(idx),fontsize=11)

        ax3.set_xlim([0, img_params[survey]['IMG_WIDTH_NEW']])
        ax3.set_ylim([img_params[survey]['IMG_HEIGHT_NEW'], 0])
        ax3.set_title(zid)
        ax3.set_aspect('equal')
        
        ax4.set_xlim([0, img_params[survey]['IMG_WIDTH_NEW']])
        ax4.set_ylim([img_params[survey]['IMG_HEIGHT_NEW'], 0])
        ax4.set_title('Consensus ({0:d}/{1:d} users)'.format(cons['n_users'],cons['n_total']))

        ax4.text(15,50,r'$n_\mathrm{src}=$'+'{0}'.format(len(answer)),fontsize=12)
        
        ax4.set_aspect('equal')
        
        # Display IR and radio images
        
        im_standard = grab_image(sub,imgtype='standard')

        ax1.imshow(im_standard,origin='upper')
        ax1.set_title('weights = {0}'.format(w))

        im_radio = grab_image(sub,imgtype='radio')

        ax2.imshow(im_radio,origin='upper')
        ax2.set_title(sub['metadata']['source'])
        ax2.get_yaxis().set_ticklabels([])

        ax3.get_yaxis().set_ticklabels([])

        # Plot contours identified as the consensus
        if len(answer) > 0:
            ax4.add_patch(patch_black)
        ax4.yaxis.tick_right()

        nticks = 5
        
        ax1.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
        ax2.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
        ax3.get_xaxis().set_ticks(np.arange(nticks)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)
        ax4.get_xaxis().set_ticks(np.arange(nticks+1)*img_params[survey]['IMG_WIDTH_NEW'] * 1./nticks)

    plt.subplots_adjust(wspace=0.02)
    
    # Save hard copy of the figure
    if savefig == True:
        fig.savefig('{0}/csv/weights/plots/{1}.pdf'.format(rgz_path,zid))
    else:
        plt.show()

    # Close figure after it's done; otherwise mpl complains about having thousands of stuff open
    plt.close()

    return None
