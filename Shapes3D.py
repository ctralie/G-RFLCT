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
	
	def getDiagLength(self):
		dX = self.XLen()/2
		dY = self.YLen()/2
		dZ = self.ZLen()/2
		return math.sqrt(dX*dX + dY*dY + dZ*dZ)
	
	def getCenter(self):
		return Point3D((self.xmax+self.xmin)/2.0, (self.ymax+self.ymin)/2.0, (self.zmax+self.zmin)/2.0)
	
	def addPoint(self, P):
		if P.x < self.xmin: self.xmin = P.x
		if P.x > self.xmax: self.xmax = P.x
		if P.y < self.ymin: self.ymin = P.y
		if P.y > self.ymax: self.ymax = P.y
		if P.z < self.zmin: self.zmin = P.z
		if P.z > self.zmax: self.zmax = P.z

class Circle3D(object):
	#center is a point
	#R is a vector whose magnitude is the radius length
	#normal is a vector to orient the circle
	def __init__(self, center, R, normal):
		self.center = center
		self.R = R
		self.normal = normal
	
	#This method assumes the circles are coplanar
	def intersectsOtherCircle(self, other):
		dCenters = self.center - other.center
		R1 = self.R.Length()
		R2 = other.R.Length()
		if dCenters.Length() < R1 + R2 + EPS:
			return True
		return False

#Implement a *slow* N^4 version (since N will be small most of the time)
def getMinimumEnclosingCircle(points, normal):
	minimum = None
	minRadiusSqr = float("inf")
	for i in range(0, len(points)):
		for j in range(i+1, len(points)):
			for k in range(j+1, len(points)):
				[P0, P1, P2] = [points[i], points[j], points[k]]
				#Construct circle by intersecting perpendicular bisectors
				#between pairs of points to find the center
				bisector1 = Line3D(P0 + 0.5*(P1-P0), (P1-P0)%normal)
				bisector2 = Line3D(P1 + 0.5*(P2-P1), (P2-P1)%normal)
				center = bisector1.intersectOtherLine(bisector2)
				#if not isinstance(center, Point3D):
				#	print "i = %i, j = %i, k = %i"%(i, j, k)
				#	print "P0: %s\nP1: %s\nP2: %s"%(P0, P1, P2)
				#	print "Normal: %s\n(P1 - P0): %s\n(P2 - P1): %s"%(normal, P1 - P0, P2 - P1)
				#	print "bisector1: %s\nbisector2: %s"%(bisector1, bisector2)
				#	print "center: %s\n\n"%center
				if isinstance(center, Point3D):
					#Check to make sure all other points lie inside the circle
					RVec = P0 - center
					RLenSqr = RVec.squaredMag()
					#print "\n\n\n"
					#print "%s\t%s\t%s"%(P0, P1, P2)
					#print "Center = %s, R = %g"%(center, RVec.Length())
					if RLenSqr < minRadiusSqr:
						allInside = True
						for P in points:
							if (P - center).squaredMag() > RLenSqr+EPS:
								#print "%s not inside"%P
								allInside = False
								break
						if allInside:
							minimum = Circle3D(center, RVec, normal)
	return minimum

if __name__ == '__main__':
	normal = Vector3D(-1, 0, -0)
	points = [Point3D(-2.5, 2.5, 2.5), Point3D(-2.5, 2.5, -2.5), Point3D(-2.5, -2.5, -2.5), Point3D(-2.5, -2.5, 2.5)]
	circle = getMinimumEnclosingCircle(points, normal)
	print circle
	#points = [Point3D(-1, 0, 0), Point3D(0, 1, 0), Point3D(1, 0, 0), Point3D(0, 0, 0), Point3D(0.1, 1.5, 0)];
	#PointRot = Point3D(0, 0, 0)
	#AxisRot = Vector3D(1, 1, 1)
	#Angle = 0.5
	#points = [rotateAroundAxis(PointRot, AxisRot, Angle, P) for P in points]
	#circle = getMinimumEnclosingCircle(points, (points[1] - points[0])%(points[2] - points[1]))
	#circle.center = rotateAroundAxis(PointRot, AxisRot, -Angle, P)
	#print circle.center
	#for P in points:
	#	print (P - circle.center).Length()
	#print circle.R.Length()
