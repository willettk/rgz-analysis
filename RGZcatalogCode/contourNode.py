import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.special import erfinv
from matplotlib import path

#tree implementation for contours
class Node(object):

    #initialize tree
    def __init__(self, value=None, contour=None, fits_loc=None, img=None, w=None, sigma=0):
        self.value = value #contour curve and level data
        self.children = [] #next contour curves contained within this one
        if fits_loc is not None:
            self.getFITS(fits_loc)
        else:
            self.img = img #FITS data as an array
            self.w = w #WCS converter object
        dec = self.w.wcs_pix2world( np.array( [[66, 66]] ), 1)[0][1] #dec of image center
        if dec > 4.5558: #northern region, above +4*33'21"
            self.beamArea = 1.44*np.pi*5.4*5.4/4 #5.4" FWHM circle
        elif 4.5558 > dec > -2.5069: #middle region, between +4*33'21" and -2*30'25"
            self.beamArea = 1.44*np.pi*6.4*5.4/4 #6.4"x5.4" FWHM ellipse
        else: #southern region, below -2*30'25"
            self.beamArea = 1.44*np.pi*6.8*5.4/4 #6.8"x5.4" FWHM ellipse
        self.pixelArea = wcs.utils.proj_plane_pixel_area(self.w)*3600*3600 #arcsecond^2
        if contour is not None:
            mad2sigma = np.sqrt(2)*erfinv(2*0.75-1) #conversion factor
            self.sigma = (contour[0]['level']/3) / mad2sigma #standard deviation of flux density measurements
            for i in contour:
                self.insert(Node(value=i, img=self.img, w=self.w, sigma=self.sigma))
            vertices = []
            for pos in contour[0]['arr']:
                vertices.append([pos['x'], pos['y']])
            self.pathOutline = path.Path(vertices) #self.pathOutline is a Path object tracing the contour
            self.getTotalFlux() #self.flux and self.fluxErr are the total integrated flux and error, respectively; also sets self.area, in arcsec^2
            self.getPeaks() #self.peaks is list of dicts of peak fluxes and locations
        else:
            self.sigma = sigma
            self.pathOutline = None
            self.flux = 0
            self.fluxErr = 0
            self.peaks = []

    #insert a contour node
    def insert(self, newNode):
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

    #manually check the topology of the tree by printing level numbers and bboxes to screen
    #for testing only
    def check(self):
        if self.value is None:
            print 'Empty'
        else:
            print 'Level ' + str(self.value['k']) + ': ' + str(fn.findBox(self.value['arr']))
            if self.children == []:
                print 'End'
            else:
                for child in self.children:
                    child.check()

    #get FITS data from file
    def getFITS(self, fits_loc):
        self.img = fits.getdata(fits_loc, 0) #imports data as array
        self.img[np.isnan(self.img)] = 0 #sets NANs to 0
        self.w = wcs.WCS(fits.open(fits_loc)[0].header) #gets pixel-to-WCS conversion from header
        return self.img

    #find the total integrated flux of the component and its error
    def getTotalFlux(self):
        fluxDensity = 0
        pixelCount = 0
        for i in range(131):
            for j in range(131):
                if self.contains([i+1, j+1]):
                    fluxDensity += self.img[j][i]
                    pixelCount += 1
        self.area = pixelCount*self.pixelArea
        fluxDensityErr = np.sqrt(pixelCount)*self.sigma
        self.flux = fluxDensity*self.pixelArea/self.beamArea*1000
        self.fluxErr = fluxDensityErr*self.pixelArea/self.beamArea*1000
        return [self.flux, self.fluxErr]
        
    #finds the peak values (in mJy) and locations (in DS9 pixel space) and return as dict
    def getPeaks(self, pList=None):
        if pList is None:
            pList = []
        if self.children == []:
            bbox = fn.bboxToDS9(fn.findBox(self.value['arr']))[0] #bbox of innermost contour
            flux = self.img[ bbox[3]-1:bbox[1]+1, bbox[2]-1:bbox[0]+1 ].max() #peak flux in bbox, with 1 pixel padding
            locP = np.where(self.img == flux) #location in pixels
            locRD = self.w.wcs_pix2world( np.array( [[locP[1][0]+1, locP[0][0]+1]] ), 1) #location in ra and dec
            peak = dict(ra=locRD[0][0], dec=locRD[0][1], peakFlux=flux*self.pixelArea/self.beamArea*1000)
            pList.append(peak)
        else:
            for child in self.children:
                child.getPeaks(pList)
        self.peaks = pList
        return self.peaks

    #returns 1 if point is within the contour, returns 0 if otherwise or if there is no contour data
    def contains(self, point):
        if self.pathOutline is not None:
            return self.pathOutline.contains_point(point)
        else:
            return 0

#creates a bounding box for a given contour path
#loop = data['contours'][0][0]['arr'] #outermost contour (for testing)
def findBox(loop):
   xmax = loop[0]['x']
   ymax = loop[0]['y']
   xmin = loop[0]['x']
   ymin = loop[0]['y']
   for i in loop:
      if i['x']>xmax:
         xmax = i['x']
      elif i['x']<xmin:
         xmin = i['x']
      if i['y']>ymax:
         ymax = i['y']
      elif i['y']<ymin:
         ymin = i['y']
   return [xmax, ymax, xmin, ymin]

#finds the coordinates of the bbox in DS9's system and the input values for drawing a box in DS9
#bbox = tree.value['bbox'] #outermost bbox (for testing)
def bboxToDS9(bbox):
    xmax = bbox[0]
    ymax = bbox[1]
    xmin = bbox[2]
    ymin = bbox[3]
    temp = 133-ymax
    ymax = 133-ymin
    ymin = temp
    newBbox = [xmax, ymax, xmin, ymin]
    ds9Box = [ (xmax+xmin)/2., (ymax+ymin)/2., xmax-xmin, ymax-ymin ]
    return [newBbox, ds9Box]
