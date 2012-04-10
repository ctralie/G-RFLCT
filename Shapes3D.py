from Primitives3D import *
from OpenGL.GL import *
import sys


class BBox3D(object):
	def __init__(self, xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1):
		[self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax] = [xmin, xmax, ymin, ymax, zmin, zmax]
	
	def XLen(self):
		return self.xmax - self.xmin
	
	def YLen(self):
		return self.ymax - self.ymin
	
	def ZLen(self):
		return self.zmax - self.zmin
	
	def getCenter(self):
		return Point3D((self.xmax+self.xmin)/2.0, (self.ymax+self.ymin)/2.0, (self.zmax+self.zmin)/2.0)
	
	def addPoint(self, P):
		if P.x < self.xmin: self.xmin = P.x
		if P.x > self.xmax: self.xmax = P.x
		if P.y < self.ymin: self.ymin = P.y
		if P.y > self.ymax: self.ymax = P.y
		if P.z < self.zmin: self.zmin = P.z
		if P.z > self.zmax: self.zmax = P.z
