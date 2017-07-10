'''
This file contains an implementation of a contour tree object. Each Node contains a contour (stored as a matplotlib Path) and links to the immediately interior contours.
'''

import numpy as np
from matplotlib import path
import matplotlib.pyplot as plt
import catalog_functions as fn # Contains custom functions

class Node(object):
	'''Tree implementation for contours'''
	
	def __init__(self, w, value=None, contour=None):
		'''Tree initializer'''
		self.y = w._naxis2
		self.value = value # Contour curve and level data
		self.children = [] # Next contour curves contained within this one
		if contour is not None:
			for c in contour:
				self.insert(Node(w, value=c))
	
	def insert(self, new_node):
		'''Insert a contour node'''
		if self.value is None: # Initialize the root with the outermost contour
			self.value = new_node.value
			vertices = []
			for pos in self.value['arr']:
				vertices.append([pos['x'], 1.+self.y-pos['y']])
			self.path = path.Path(vertices) # self.path is a Path object tracing the contour
		elif self.value == new_node.value: # No duplicate contours
			return
		else:
			vertices = []
			for pos in new_node.value['arr']:
				vertices.append([pos['x'], 1.+self.y-pos['y']])
			new_node.path = path.Path(vertices) # self.path is a Path object tracing the contour
			if new_node.value['k'] == self.value['k'] + 1: # Add a contour one level higher as a child
				self.children.append(new_node)
			elif new_node.value['k'] <= self.value['k']: # If a contour of lower level appears, something went wrong
				raise RuntimeError('Inside-out contour')
			else: # Otherwise, find the next level that has a bounding box enclosing the new contour
				inner = fn.findBox(new_node.value['arr'])
				for child in self.children:
					outer = fn.findBox(child.value['arr'])
					if outer[0]>inner[0] and outer[1]>inner[1] and outer[2]<inner[2] and outer[3]<inner[3]:
						child.insert(new_node)
	
	def check(self):
		'''Manually check the topology of the tree by printing level numbers and bboxes to screen (for testing only)'''
		if self.value is None:
			print 'Empty'
		else:
			print 'Level {}: {}'.format(self.value['k'], fn.findBox(self.value['arr']))
			if self.children == []:
				print 'End'
			else:
				for child in self.children:
					child.check()
	
	def contains(self, point):
		'''Returns boolean or array of booleans'''
		if self.path is not None:
			if len(point.shape) == 1:
				return np.array([self.path.contains_point(point)])
			else:
				return self.path.contains_points(point)
		else:
			return np.array([False])
	
	def print_contours(self, skip=0, **kwargs):
		'''Print every contour with matplotlib
		Skips denotes how many contours to ignore when plotting'''
		assert type(skip) == int, 'Skip parameter must be int'
		if self.value is None:
			print 'Empty'
		else:
			if skip:
				skip -= 1
			else:
				plt.plot(self.path.vertices.T[0], self.path.vertices.T[1], **kwargs)
			for child in self.children:
				child.print_contours(skip, **kwargs)
	
	def get_levels(self, max_k=0):
		'''Returns how many levels deep this tree goes'''
		if self.children == []:
			return self.value['k']
		else:
			for child in self.children:
				max_k = max(max_k, child.get_levels(max_k))
			return max_k
	
	def print_contour_levels(self):
		'''Print the values of the contours to screen'''
		if 'level' in self.value:
			print self.value['level']
		for child in self.children:
			child.print_contour_levels()
	
	def get_min_disjoint(self, point, roots):
		'''Travels down each branch until it reaches a child that doesn't contain the supplied point'''
		if not self.contains(point):
			roots.append(self)
		elif self.children != []:
			for child in self.children:
				child.get_min_disjoint(point, roots)
	
	def get_kth_level(self, k, roots):
		'''Returns all contours with level k'''
		if self.value['k'] == k:
			roots.append(self)
		elif self.value['k'] < k and self.children != []:
			for child in self.children:
				child.get_kth_level(k, roots)
	
	def get_equal_disjoint(self, point, roots):
		'''Finds the smallest k such that the supplied point is outside all contours and returns all roots at that level'''
		temp_roots = []
		self.get_min_disjoint(point, temp_roots)
		max_k = []
		for root in temp_roots:
			max_k.append(root.value['k'])
		if len(max_k) == 0:
			return
		self.get_kth_level(max(max_k), roots)
		
		# Check that one of the paths didn't include the point at all levels
		delete_roots = []
		new_max_k = max(max_k)
		for ix, root in enumerate(roots):
			if root.contains(point):
				delete_roots.append(ix)
				new_max_k = max(new_max_k, root.get_levels()+1)
		
		# If one or more did, update roots to exclude them
		if new_max_k > max(max_k):
			temp_roots = [root for ix, root in enumerate(roots) if ix not in delete_roots]
			roots[:] = []
			for root in temp_roots:
				subroots = []
				root.get_kth_level(new_max_k, subroots)
				roots[:] = roots + subroots
	
	def remove_triple_center(self, ir_pos, peak_pos):
		'''Finds the center component of a triple source (i.e. outermost contour that contains the IR host and a radio peak) and removes it from the tree'''
		for ix, child in enumerate(self.children):
			if child.contains(ir_pos):
				if sum(child.contains(peak_pos))==1:
					del self.children[ix]
				elif np.any(child.contains(peak_pos)):
					child.remove_triple_center(ir_pos, peak_pos)
