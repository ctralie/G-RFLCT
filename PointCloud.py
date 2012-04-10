import os
from Primitives3D import *

class PointCloud(object):
	def __init__(self):
		self.points = []
		self.colors = []
	
	def loadXYZFile(self, filename)
		fin = fopen(filename, 'r')
		lineNum = 1
		for line in fin:
			fields = line.split()
			if len(fields) < 3:
				print "ERROR: Too few fields detected on line %i of %s"%(lineNum, filename)
			else:
				self.points.append(Point3D(*fields[0:3]))
				if len(fields) > 3:
					self.colors.append(fields[3:])
			lineNum = lineNum+1
	
	def loadFile(self, filename):
		suffix = filename.rsplit('.')[-1]
		if suffix[-1] == '\n':
			suffix = suffix[0:-1]
		if suffix == "xyz":
			self.loadXYZFile(filename)
		else:
			print "ERROR: Unrecognized file extension %s"%(suffix)
	
	def Draw(self, filename):
		
