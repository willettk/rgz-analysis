"""
2) PEAKS POSITIONS
I'm working with the ['radio']['peaks'] positions and I noticed an
interesting feature. For entries in the catalog with
['radio']['numberPeaks']=2, I plotted the difference in ra e in dec of
the two peaks. The resulting plot is pp.png and it clearly shows a
pattern. What confuses me is that the discretization happens only along
the declination axis. This is probably related to how the flux gradient
is defined, but I though you should know anyway. In fact, I'd like you
to confirm the behaviour so I know it's not a bug in my code.
"""

# MongoDB parameters

from pymongo import MongoClient
from matplotlib import pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

client = MongoClient('localhost', 27017)
db = client['radio'] 

#subjects = db['radio_subjects'] 		# subjects = images
#classifications = db['radio_classifications']	# classifications = classifications of each subject per user
#consensus = db['consensus']	# classifications = classifications of each subject per user
catalog = db['catalog']	# classifications = classifications of each subject per user


twopeaks = catalog.find({'radio.numberPeaks':2})

'''
diffarr,radiff,decdiff = [],[],[]
for source in twopeaks:
    peak1 = source['radio']['peaks'][0]
    peak2 = source['radio']['peaks'][1]

    c1 = SkyCoord(peak1['ra'],peak1['dec'],unit=u.deg,frame='icrs')
    c2 = SkyCoord(peak2['ra'],peak2['dec'],unit=u.deg,frame='icrs')

    diffarr.append(c1.separation(c2).arcsecond)
    radiff.append((c1.ra - c2.ra).arcsecond)
    decdiff.append((c1.dec - c2.dec).arcsecond)
'''

ra1,dec1 = [],[]
ra2,dec2 = [],[]
for source in twopeaks:
    peak1 = source['radio']['peaks'][0]
    peak2 = source['radio']['peaks'][1]

    ra1.append(peak1['ra'])
    dec1.append(peak1['dec'])
    ra2.append(peak2['ra'])
    dec2.append(peak2['dec'])

carr1 = SkyCoord(ra1,dec1,unit='deg')
carr2 = SkyCoord(ra2,dec2,unit='deg')

diffarr = carr2.separation(carr1).arcsec
radiff = (carr1.ra - carr2.ra).arcsec
decdiff = (carr1.dec - carr2.dec).arcsec

fig,axarr = plt.subplots(1,3,figsize=(15,6))
ax1,ax2,ax3 = axarr.ravel()

ax1.scatter(radiff,decdiff,s=3,color='r')
ax1.set_xlim(-30,30)
ax1.set_ylim(-30,30)
ax1.set_xlabel(r"$\Delta$ RA [arcsec]",fontsize=20)
ax1.set_ylabel(r"$\Delta$ dec [arcsec]",fontsize=20)

ax2.hist(radiff,bins=100,range=(-10,10))
ax2.set_xlabel(r"$\Delta$ RA [arcsec]",fontsize=20)
ax3.hist(decdiff,bins=100,range=(-10,10))
ax3.set_xlabel(r"$\Delta$ dec [arcsec]",fontsize=20)

plt.show()
