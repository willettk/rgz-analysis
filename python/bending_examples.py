# Make a gallery of images showing the RGZ consensus double sources, sorted by bending angle.

from astropy.io import ascii

path = '/Users/willettk/Astronomy/Research/GalaxyZoo'
data = ascii.read('{:}/rgz-analysis/csv/static_catalog3.csv'.format(path),delimiter=' ')

import bending_angles as ba
import numpy as np

pathdict = ba.make_pathdict()

def bending_examples():
    for a in np.linspace(0,80,9):
        bdata = data[(data['bending_angle'] >= a) & (data['bending_angle'] < a+10.)]
        count,errcount = 0,0
        if len(bdata) > 0:
            for b in bdata:
                zid = b['zooniverse_id']
                try:
                    if b['angle_type'] == 'multipeaked_singles':
                        angle_type = 'mps'
                    else:
                        angle_type = 'radio'
                    ba.plot_one_double(zid,pathdict,savefig=True,anglepath='{0:.0f}_{1:.0f}/'.format(a,a+10),dbltype=angle_type)
                    count += 1
                except ValueError as inst:
                    print "ValueError,",inst.args,zid
                    errcount += 1
        print '{:d} galaxies with bending angle, {:d} with errors for angles between {:.0f} and {:.0f} deg'.format(count,errcount,a,a+10)

def triples():

    # Triples

    data = ascii.read('{:}/rgz-analysis/bending_angles/angles_triple_pixradio.csv'.format(path),delimiter=',')

    print '{:d} triple galaxies with a bending angle'.format(len(data))

    for a in np.linspace(0,170,18):
        bdata = data[(data['bending_angle'] >= a) & (data['bending_angle'] < a+10.)]
        if len(bdata) > 0:
            for t in bdata:
                zid = t['zooniverse_id']
                ba.plot_one_triple(zid,pathdict,savefig=True,anglepath='{0:.0f}_{1:.0f}/'.format(a,a+10))
        print '{:d} triple galaxies with bending angle for angles between {:.0f} and {:.0f} deg'.format(len(bdata),a,a+10)

