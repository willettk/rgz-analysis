'''
This file contains an implementation of a contour tree object. Each Node contains a contour and links to the immediately interior contours.
The total radio flux through a contour, as well as radio peaks, can be calculated from the FITS file plus the contour pixels.
'''

import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.special import erfinv
from matplotlib import path
import catalogFunctions as fn #contains custom functions

class Node(object):
    '''tree implementation for contours'''
    
    def __init__(self, value=None, contour=None, fits_loc=None, img=None, w=None, sigmaJyBeam=0):
        '''tree initializer'''
        self.value = value #contour curve and level data
        self.children = [] #next contour curves contained within this one
        if fits_loc is not None:
            self.getFITS(fits_loc)
        else:
            self.img = img #FITS data as an array
            self.imgSize = int(img.shape[0]) #size in pixels of FITS data
            self.w = w #WCS converter object
        dec = self.w.wcs_pix2world( np.array( [[self.imgSize/2, self.imgSize/2]] ), 1)[0][1] #dec of image center
        if dec > 4.5558: #northern region, above +4*33'21"
            self.beamAreaArcsec2 = 1.44*np.pi*5.4*5.4/4 #5.4" FWHM circle
        elif 4.5558 > dec > -2.5069: #middle region, between +4*33'21" and -2*30'25"
            self.beamAreaArcsec2 = 1.44*np.pi*6.4*5.4/4 #6.4"x5.4" FWHM ellipse
        else: #southern region, below -2*30'25"
            self.beamAreaArcsec2 = 1.44*np.pi*6.8*5.4/4 #6.8"x5.4" FWHM ellipse
        self.pixelAreaArcsec2 = wcs.utils.proj_plane_pixel_area(self.w)*3600*3600 #arcsecond^2
        if contour is not None:
            mad2sigma = np.sqrt(2)*erfinv(2*0.75-1) #conversion factor
            self.sigmaJyBeam = (contour[0]['level']/3) / mad2sigma #standard deviation of flux density measurements
            self.sigmamJy = self.sigmaJyBeam*1000*self.pixelAreaArcsec2/self.beamAreaArcsec2
            for i in contour:
                self.insert(Node(value=i, img=self.img, w=self.w, sigmaJyBeam=self.sigmaJyBeam))
            vertices = []
            for pos in contour[0]['arr']:
                vertices.append([pos['x'], pos['y']])
            self.pathOutline = path.Path(vertices) #self.pathOutline is a Path object tracing the contour
            self.getTotalFlux() #self.fluxmJy and self.fluxErrmJy are the total integrated flux and error, respectively; also sets self.areaArcsec2, in arcsec^2
            self.getPeaks() #self.peaks is list of dicts of peak fluxes and locations
        else:
            self.sigmaJyBeam = sigmaJyBeam
            self.pathOutline = None
            self.fluxmJy = 0
            self.fluxErrmJy = 0
            self.peaks = []
    
    def insert(self, newNode):
        '''insert a contour node'''
        if self.value is None: #initialize the root with the outermost contour
            self.value = newNode.value
        elif self.value == newNode.value: #no duplicate contours
            return
        else:
            if newNode.value['k'] == self.value['k'] + 1: #add a contour one level higher as a child
                self.children.append(newNode)
            elif newNode.value['k'] <= self.value['k']: #if a contour of lower level appears, something went wrong
                raise Exception('Inside-out contour')
            else: #otherwise, find the next level that has a bounding box enclosing the new contour
                inner = fn.findBox(newNode.value['arr'])
                for child in self.children:
                    outer = fn.findBox(child.value['arr'])
                    if outer[0]>inner[0] and outer[1]>inner[1] and outer[2]<inner[2] and outer[3]<inner[3]:
                        child.insert(newNode)
    
    def check(self):
        '''manually check the topology of the tree by printing level numbers and bboxes to screen (for testing only)'''
        if self.value is None:
            print 'Empty'
        else:
            print 'Level ' + str(self.value['k']) + ': ' + str(fn.findBox(self.value['arr']))
            if self.children == []:
                print 'End'
            else:
                for child in self.children:
                    child.check()
    
    def getFITS(self, fits_loc):
        '''read FITS data from file'''
        self.img = fits.getdata(fits_loc, 0) #imports data as array
        self.img[np.isnan(self.img)] = 0 #sets NANs to 0
        self.imgSize = int(self.img.shape[0])
        self.w = wcs.WCS(fits.open(fits_loc)[0].header) #gets pixel-to-WCS conversion from header
        return self.img
    
    def getTotalFlux(self):
        '''find the total integrated flux of the component and its error'''
        fluxDensityJyBeam = 0
        pixelCount = 0
        for i in range(self.imgSize-1):
            for j in range(self.imgSize-1):
                if self.contains([i+1, self.imgSize-j]):
                    fluxDensityJyBeam += self.img[j][i]
                    pixelCount += 1
        self.areaArcsec2 = pixelCount*self.pixelAreaArcsec2
        fluxDensityErrJyBeam = np.sqrt(pixelCount)*self.sigmaJyBeam
        self.fluxmJy = fluxDensityJyBeam*1000*self.pixelAreaArcsec2/self.beamAreaArcsec2
        self.fluxErrmJy = fluxDensityErrJyBeam*1000*self.pixelAreaArcsec2/self.beamAreaArcsec2
        return [self.fluxmJy, self.fluxErrmJy]
    
    def getPeaks(self, pList=None):
        '''finds the peak values (in mJy) and locations (in DS9 pixel space) and return as dict'''
        if pList is None:
            pList = []
        if self.children == []:
            bbox = fn.bboxToDS9(fn.findBox(self.value['arr']), self.imgSize)[0] #bbox of innermost contour
            fluxDensityJyBeam = self.img[ bbox[3]-1:bbox[1]+1, bbox[2]-1:bbox[0]+1 ].max() #peak flux in bbox, with 1 pixel padding
            locP = np.where(self.img == fluxDensityJyBeam) #location in pixels
            locRD = self.w.wcs_pix2world( np.array( [[locP[1][0]+1, locP[0][0]+1]] ), 1) #location in ra and dec
            peak = dict(ra=locRD[0][0], dec=locRD[0][1], flux=fluxDensityJyBeam*1000)
            pList.append(peak)
        else:
            for child in self.children:
                child.getPeaks(pList)
        self.peaks = pList
        return self.peaks
    
    def contains(self, point):
        '''returns 1 if point is within the contour, returns 0 if otherwise or if there is no contour data'''
        if self.pathOutline is not None:
            return self.pathOutline.contains_point(point)
        else:
            return 0
