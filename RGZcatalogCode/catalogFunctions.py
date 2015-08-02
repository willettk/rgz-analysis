import numpy as np
import pandas as pd
from StringIO import StringIO
import mechanize

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

#determines if two floats are approximately equal
def approx(a, b, uncertainty=1e-5):
   return np.abs(a-b) < uncertainty

#pass an SQL query to SDSS and return a pandas dataframe
def SDSS_select(sql):
   br = mechanize.Browser()
   try:
      br.open('http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx', timeout=4)
   except (mechanize.HTTPError, mechanize.URLError): #try once more
      br.open('http://skyserver.sdss.org/dr12/en/tools/search/sql.aspx', timeout=4)
   br.select_form(name='sql')
   br['cmd'] = sql
   br['format'] = ['csv']
   response = br.submit()
   file_like = StringIO(response.get_data())
   return pd.read_csv(file_like, skiprows=1)
