from Primitives3D import *
from Shapes3D import *
from OpenGL.GL import *
from Graphics3D import *

class PointCloud(object):
	def __init__(self):
		self.points = []
		self.colors = []
		self.needsDisplayUpdate = True
		self.DisplayList = -1
	
	def loadXYZFile(self, filename):
		fin = fopen(filename, 'r')
		lineNum = 1
		for line in fin:
			fields = line.split()
			if len(fields) < 3:
				print "ERROR: Too few fields detected on line %i of %s"%(lineNum, filename)
			else:
				self.points.append(Point3D(*fields[0:3]))
				if len(fields) > 3:
					self.colors.append([float(c) for c in fields[3:]])
			lineNum = lineNum+1
	
	def loadFile(self, filename):
		suffix = filename.rsplit('.')[-1]
		if suffix[-1] == '\n':
			suffix = suffix[0:-1]
		if suffix == "xyz":
			self.loadXYZFile(filename)
		else:
			print "ERROR: Unrecognized file extension %s"%(suffix)
	
	def renderGL(self):
		if self.needsDisplayUpdate:
			if self.DisplayList != -1: #Deallocate previous display list
				glDeleteLists(self.DisplayList, 1)
			self.DisplayList = glGenLists(1)
			glNewList(self.DisplayList, GL_COMPILE)
			glPointSize(5)
			glBegin(GL_POINTS)
			for P, C in zip(self.points, self.colors):
				glColor3f(C[0], C[1], C[2])
				glVertex3f(P.x, P.y, P.z)
			glEnd()
			glEndList()
			self.needsDisplayUpdate = False
		glCallList(self.DisplayList)

	def getBBox(self):
		if len(self.vertices) == 0:
			return BBox3D(0, 0, 0, 0, 0, 0)
		P0 = self.points[0]
		bbox = BBox3D(P0.x, P0.x, P0.y, P0.y, P0.z, P0.z)
		for P in self.points:
			bbox.addPoint(P)
		return bbox
