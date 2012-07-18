from Primitives3D import *
from Shapes3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
import math

#TODO: Update image sources pruning to make use of this class
#Get the epsilon to be used for numerical precision
def getEPS(A, B, C):
	avgdx = (abs(A.x-B.x) + abs(A.x-C.x) + abs(B.x-C.x))/3.0
	avgdy = (abs(A.y-B.y) + abs(A.y-C.y) + abs(B.y-C.y))/3.0
	avgdz = (abs(A.z-B.z) + abs(A.z-C.z) + abs(B.z-C.z))/3.0
	return min(min(avgdx, avgdy), avgdz)*1e-7

def PointsEqual2D(P1, P2, eps):
	if (abs(P1.x-P2.x) < eps and abs(P1.y-P2.y) < eps):
		return True
	return False

def CCW2D(A, B, C):
	det = B.x*C.y - B.y*C.x - A.x*C.y + A.y*C.x + A.x*B.y - A.y*B.x
	eps = getEPS(A, B, C)
	if (det > eps):
		return -1
	elif (det < -eps):
		return 1
	#Are the points all equal?
	if (PointsEqual2D(A, B, eps) and PointsEqual2D(B, C, eps)):
		return 0
	if (PointsEqual2D(A, B, eps)):
		return 2
	#Is C in the closure of A and B?
	#Vectors must be in opposite directions or one of the vectors
	#must be zero (C is on one of the endpoints of A and B)
	vAC = C - A
	vBC = C - B
	vAC.z = 0
	vBC.z = 0
	if (vAC.Dot(vBC) <= 0):
		return 0;#This fires for C in the closure of A and B (including endpoints)
	vBA = A - B
	vBA.z = 0
	#C to the left of AB
	if (vBA.Dot(vBC) > 0):
		return -2
	#C to the right of AB
	else:
		return 2

#Each beam should have a coordinate system and a near distance
#The coordinate system should orient the "towards" direction
#perpendicular to the image plane of the beam
#"face" points to the MeshFace associated with the beam, if there is one
#Otherwise, the beam it is one of the initial cast beams and has not 
#reflected across a face yet. In that case the default neardist can be used
#Other member variables: 
#mvMatrix: Modelview matrix for putting in the coordinate system of a beam
#frustPoints: Points of the frustum on the 2D image plane
class Beam3D(object):
	def __init__(self, origin, frustVertices, face = None):
		self.origin = origin
		self.neardist = 0.01
		self.face = face
		if len(frustVertices) < 3:
			print "ERROR: Less than 3 frustVertices on beam projection face"
			return
		self.towards = getFaceNormal(frustVertices)
		dV = frustVertices[0] - origin
		if dV.Dot(self.towards) < 0:
			self.towards = -1*self.towards
		self.right = frustVertices[1] - frustVertices[0]
		self.right.normalize()
		self.up = self.right%self.towards
		self.up.normalize()
		#Now calculate and store the modelview matrix to get into the 
		#beam's coordinate system
		#(NOTE this function follows the OpenGL convention that things in
		#front of the beam are -z)
		self.mvMatrix = getCameraMatrix(self.towards, self.up, self.right, self.origin)
		#Now calculate the near distance if the beam is associated with a face
		if face:
			P = self.mvMatrix*frustVertices[0]
			self.neardist = -P.z
		#Now map the frustum points to 2D image plane coordinates
		self.frustPoints = self.projectPolygon(frustVertices, False)
	
	#Project a polygon onto the beam's image plane (before clipping to frustum)
	def projectPolygon(self, polygon, doClip = True):
		#First transform all points into the beam's field of view
		mvVerts = [self.mvMatrix*v for v in polygon]
		#Now clip the polygon to the near plane
		clippedVerts = []
		if doClip:
			for i in range(0, len(mvVerts)):
				v1 = mvVerts[i]
				v2 = mvVerts[(i+1)%len(mvVerts)]
				(d1, d2) = (-v1.z, -v2.z)
				#v1 is behind plane and v2 is in front of plane
				if d1 < self.neardist and d2 >= self.neardist:
					ratio = (self.neardist-d1)/(d2-d1)
					vNew = Vector3D(v1.x+(v2.x-v1.x)*ratio, v1.y+(v2.y-v1.y)*ratio, -self.neardist)
					clippedVerts.append(vNew)
				#v1 is in front of plane and v2 is behind plane
				elif d1 >= self.neardist and d2 < self.neardist:
					ratio = (self.neardist-d2)/(d1 - d2)
					vNew = Vector3D(v2.x+(v1.x-v2.x)*ratio, v2.y+(v1.y-v2.y)*ratio, -self.neardist)
					clippedVerts.append(v1)
					clippedVerts.append(vNew)
				#Both vertices are in front of the plane
				elif d1 >= self.neardist and d2 >= self.neardist:
					#Just add the left one
					clippedVerts.append(v1)
				#Otherwise two vertices are behind
		else:
			clippedVerts = mvVerts
		#FOR DEBUGGING
		#for v in clippedVerts:
		#	print self.mvMatrix.Inverse()*v
		#Now do the perspective projection
		for i in range(0, len(clippedVerts)):
			clippedVerts[i].x = -self.neardist*clippedVerts[i].x/clippedVerts[i].z
			clippedVerts[i].y = -self.neardist*clippedVerts[i].y/clippedVerts[i].z
			clippedVerts[i].z = -self.neardist
		#Make sure the points are in counter clockwise order
		if len(clippedVerts) >= 3:
			ccw = CCW2D(clippedVerts[0], clippedVerts[1], clippedVerts[2])
			if ccw == 1:
				clippedVerts.reverse()
		return clippedVerts
	
	
	#Perform Sutherland Hodgman Clipping to clip a projected polygon to
	#the inside of the beam
	#Following pseudocode on http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
	#The function assumes that polygon2D has been clipped to the near plane
	#and put in image plane coordinates (with the use of the function projectPolygon)
	#TODO: Make sure this function can handle a 1 Point polygon, so I can use that
	#to test whether a receiver position is within a beam
	def clipToFrustum(self, polygon2D):
		outputList = polygon2D[:]
		for i in range(0, len(self.frustPoints)):
			clipEdge = [self.frustPoints[i], self.frustPoints[(i+1)%len(self.frustPoints)]]
			inputList = outputList
			outputList = []
			S = inputList[-1]
			for E in inputList:
				#E is inside the clip edge
				if CCW2D(clipEdge[0], clipEdge[1], E) != 1:
					#S is not inside the clip edge
					if CCW2D(clipEdge[0], clipEdge[1], S) == 1:
						#Polygon going from outside to inside
						if CCW2D(clipEdge[0], clipEdge[1], E) != 0:
							#Only add the intersection if E is not on the clip edge
							#(otherwise E gets added twice)
							line1 = Line3D(clipEdge[0], clipEdge[1]-clipEdge[0])
							line2 = Line3D(S, E-S)
							intersection = line1.intersectOtherLine(line2)
							if not intersection:
								print "1: Clip intersection not found: %s, %s, %s, %s"%(self.mvMatrix.Inverse()*clipEdge[0], self.mvMatrix.Inverse()*clipEdge[1], self.mvMatrix.Inverse()*S, self.mvMatrix.Inverse()*E)
							outputList.append(intersection)
					outputList.append(E)
				elif CCW2D(clipEdge[0], clipEdge[1], S) != 1:
					#Polygon going from inside to outside
					if CCW2D(clipEdge[0], clipEdge[1], S) != 0:
						#Only add intersection if S is not on the clip edge
						#(otherwise it gets added twice since it's already been added)
						line1 = Line3D(clipEdge[0], clipEdge[1]-clipEdge[0])
						line2 = Line3D(S, E-S)
						intersection = line1.intersectOtherLine(line2)
						if not intersection:
							print "2: Clip intersection not found: %s, %s, %s, %s"%(self.mvMatrix.Inverse()*clipEdge[0], self.mvMatrix.Inverse()*clipEdge[1], self.mvMatrix.Inverse()*S, self.mvMatrix.Inverse()*E)
						outputList.append(intersection)
				S = E
		#for v in outputList:
		#	print self.mvMatrix.Inverse()*v
		return outputList
	
	def clipPolygon(self, poly):
		ret = self.projectPolygon(poly)
		return self.clipToFrustum(ret)
	
	def clipMeshFace(self, face):
		return self.clipPolygon([v.pos for v in face.getVertices()])

if __name__ == '__main__':
	mesh = getRectMesh(Point3D(-1, -1, -1), Point3D(-1, 1, -1), Point3D(0, 1, -1), Point3D(0, -1, -1))
	beam = Beam3D(Point3D(0, 0, 0), [v.pos for v in mesh.faces[0].getVertices()], mesh.faces[0])
	#poly = [Vector3D(-3, 0, 0), Vector3D(0, 0, 1), Vector3D(2, 0, 1), Vector3D(3, 0, 0), Vector3D(7, 0, -4), Vector3D(5, 0, -5), Vector3D(3, 0, -5), Vector3D(0, 0, -3)]
	poly = [Vector3D(-1, 0, -1), Vector3D(-0.5, 0.5, -1), Vector3D(1.5, 0.5, -1), Vector3D(2, 0, -1), Vector3D(1.5, -0.5, -1), Vector3D(0, -1, -1)]
	poly = beam.projectPolygon(poly)
	poly = beam.clipToFrustum(poly)
