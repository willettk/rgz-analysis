import numpy as np

#shoelace algorithm to find area of polygon
#loop is a list of dicts, each containing an 'x' and 'y' coordinate
def area(loop):
   area = 0
   for i in xrange(0, len(loop)-1):
      area += loop[i]['x']*loop[i+1]['y'] - loop[i+1]['x']*loop[i]['y']
   area = 0.5 * abs(area + loop[-1]['x']*loop[0]['y'] - loop[0]['x']*loop[-1]['y'])
   return area

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

#finds the coordinates of the bbox in DS9's system
#and the imput values for drawing a box in DS9
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

#find distance between two [x,y] points
def distance(a, b):
   return np.sqrt( pow(a[0]-b[0],2) + pow(a[1]-b[1],2) )

#determines if two floats are approximately equal
def approx(a, b, uncertainty=1e-5):
   return np.abs(a-b) < uncertainty
