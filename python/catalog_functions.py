'''This file contains miscellaneous functions for the catalog pipeline, namely determine_paths, findBox, bboxToDS9, and approx.'''

import logging, os
import numpy as np

class DataAccessError(Exception):
   '''Raised when we can't access data on an external server for any reason.'''
   def __init__(self, message='Unable to connect to server'):
      self.message = message
   def __str__(self):
      return repr(self.message)

def determinePaths(paths):
   '''Set up the local data paths. Currently works from UMN servers on tabernacle, plus Kyle Willett's laptop.'''
   found_path = False
   for path in paths:
      if os.path.exists(path):
         found_path = True
         return path
   if found_path == False:
      output = 'Unable to find a hardcoded local path in %s; exiting' % str(paths)
      print output
      logging.warning(output)
      exit()

def findBox(loop):
   '''
   creates a bounding box for a given contour path
   loop = data['contours'][0][0]['arr'] #outermost contour (for testing)
   '''
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

def bboxToDS9(bbox, imgSize):
   '''
   finds the coordinates of the bbox in DS9's system and the imput values for drawing a box in DS9
   bbox = tree.value['bbox'] #outermost bbox (for testing)
   '''
   xmax = bbox[0]
   ymax = bbox[1]
   xmin = bbox[2]
   ymin = bbox[3]
   temp = imgSize+1-ymax
   ymax = imgSize+1-ymin
   ymin = temp
   newBbox = [xmax, ymax, xmin, ymin]
   ds9Box = [ (xmax+xmin)/2., (ymax+ymin)/2., xmax-xmin, ymax-ymin ]
   return [newBbox, ds9Box]

def approx(a, b, uncertainty=1e-5):
   '''determines if two floats are approximately equal'''
   return np.abs(a-b) < uncertainty
