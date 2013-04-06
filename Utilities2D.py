#I deal with a lot of special cases in my projects that involve 2D
#polygons and 2D line segments (especially in the beam tracer when
#looking in perspective space).  This file holds a lot of the 2D
#utility functions that I need to deal with 2D points.  Note that I 
#still use all of the 3D primitive classes to describe objects, they
#are just degenerate cases (e.g. Point3D(x, y, 0) describes a 2D point)
from Primitives3D import *

#This helper function is used to print 2D polygons
#as parallel lists of x and y coordinates
#so that they can be used with the "patch" command
#in Matlab
def printMatlabPoly(poly, suffix = ""):
	print "x%s = ["%suffix,
	for i in range(0, len(poly)):
		print poly[i].x,
		if i < len(poly)-1:
			print ",",
	print "]"
	print "y%s = ["%suffix,
	for i in range(0, len(poly)):
		print poly[i].y,
		if i < len(poly)-1:
			print ",",
	print "]"

#TODO: Update image sources pruning to make use of this class
#Get the epsilon to be used for numerical precision
def getEPS(A, B, C):
	return EPS
	avgdx = (abs(A.x-B.x) + abs(A.x-C.x) + abs(B.x-C.x))/3.0
	avgdy = (abs(A.y-B.y) + abs(A.y-C.y) + abs(B.y-C.y))/3.0
	avgdz = (abs(A.z-B.z) + abs(A.z-C.z) + abs(B.z-C.z))/3.0
	mins = [avgdx, avgdy, avgdz]
	mins = [x for x in mins if x > 0]
	if len(mins) > 0:
		return mins[0]*1e-4
	return 1e-7

def PointsEqual2D(P1, P2, eps):
	if (abs(P1.x-P2.x) < eps and abs(P1.y-P2.y) < eps):
		return True
	#print "P1 = %s, P2 = %s, abs(P1.x-P2.x) = %g, abs(P1.y - P2.y) = %g"%(P1, P2, abs(P1.x-P2.x), abs(P1.y-P2.y))
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
	if (vAC.Dot(vBC) < eps):
		return 0;#This fires for C in the closure of A and B (including endpoints)
	vBA = A - B
	vBA.z = 0
	#C to the left of AB
	if (vBA.Dot(vBC) > eps):
		return -2
	#C to the right of AB
	else:
		return 2

#Find the intersection of two lines segments in a numerically stable
#way by looking at them parametrically
def intersectSegments2D(A, B, C, D, countEndpoints = True):
	denomDet = (D.x-C.x)*(A.y-B.y) - (D.y-C.y)*(A.x-B.x)
	if (denomDet == 0): #Segments are parallel
		return None
	num_t = (A.x-C.x)*(A.y-B.y) - (A.y-C.y)*(A.x-B.x);
	num_s = (D.x-C.x)*(A.y-C.y) - (D.y-C.y)*(A.x-C.x);
	t = float(num_t) / float(denomDet)
	s = float(num_s) / float(denomDet)
	if (s < 0 or s > 1):
		return None #Intersection not within the bounds of segment 1
	if (t < 0 or t > 1):
		return None #Intersection not within the bounds of segment 2

	#Don't count intersections that occur at the endpoints of both segments
	#if the user so chooses
	if ((t == 0 or t == 1) and (s == 0 or s == 1) and (not countEndpoints)):
		return None

	ret = Point3D(A.x, A.y, 0)
	ret.x = ret.x + (B.x-A.x)*s;
	ret.y = ret.y + (B.y-A.y)*s;
	return ret

#Perform Sutherland Hodgman Clipping to clip the polygon "toClipPoly"
#to the inside of the polygon "boundaryPoly"
#Following pseudocode on http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
def clipSutherlandHodgman(boundaryPoly, toClipPoly):
	outputList = toClipPoly[:]
	for v in outputList:
		v.clippedVertex = False
	for i in range(0, len(boundaryPoly)):
		if len(outputList) == 0: #Special case: No Points left
			break
		clipEdge = [boundaryPoly[i], boundaryPoly[(i+1)%len(boundaryPoly)]]
		inputList = outputList
		outputList = []
		S = inputList[-1]
		for E in inputList:
			CCWS = CCW2D(clipEdge[0], clipEdge[1], S)
			CCWE = CCW2D(clipEdge[0], clipEdge[1], E)
			if CCWE != 1: #E is inside the clip edge
				if CCWS == 1: #S is not inside the clip edge
					#Polygon going from outside to inside
					if CCWE != 0:
						#Only add the intersection if E is not on the clip edge
						#(otherwise E gets added twice)
						#NOTE: Make line1 the line from which the intersection point is calculated
						#since it is a line created from the beam, which is known to have certain
						#numerical precision guarantees
						line1 = Line3D(clipEdge[0], clipEdge[1]-clipEdge[0])
						line2 = Line3D(S, E-S)
						ret = line1.intersectOtherLineRet_t(line2)
						if not ret:
							print "CCWE = %i, CCWS = %i"%(CCWE, CCWS)
							print "EPS_S = %g"%getEPS(clipEdge[0], clipEdge[1], S)
							print "EPS_E = %g"%getEPS(clipEdge[0], clipEdge[1], E)
							print "1: Clip intersection not found: ClipEdge = [%s, %s], S = %s, E = %s"%(self.mvMatrixInverse*clipEdge[0], self.mvMatrixInverse*clipEdge[1], self.mvMatrixInverse*S, self.mvMatrixInverse*E)
						else:
							(t, intersection) = ret
							intersection.clippedVertex = True
							outputList.append(intersection)
				outputList.append(E)
			elif CCWS != 1:
				#Polygon going from inside to outside
				if CCWS != 0:
					#Only add intersection if S is not on the clip edge
					#(otherwise it gets added twice since it's already been added)
					line1 = Line3D(clipEdge[0], clipEdge[1]-clipEdge[0])
					line2 = Line3D(S, E-S)
					ret = line1.intersectOtherLineRet_t(line2)
					if not ret:
						print "CCWE = %i, CCWS = %i"%(CCWE, CCWS)
						print "EPS_S = %g"%getEPS(clipEdge[0], clipEdge[1], S)
						print "EPS_E = %g"%getEPS(clipEdge[0], clipEdge[1], E)
						print "2: Clip intersection not found: ClipEdge = [%s, %s], S = %s, E = %s"%(self.mvMatrixInverse*clipEdge[0], self.mvMatrixInverse*clipEdge[1], self.mvMatrixInverse*S, self.mvMatrixInverse*E)
					else:
						(t, intersection) = ret
						intersection.clippedVertex = True
						outputList.append(intersection)
			S = E
	#for v in outputList:
	#	print self.mvMatrixInverse*v
	#Check outputList to make sure no points are overlapping
	ret = []
	for i in range(0, len(outputList)):
		if not PointsEqual2D(outputList[i], outputList[(i+1)%len(outputList)], EPS):
			ret.append(outputList[i])
	return ret



#Helper function for "pointInsideTriangle()"
def pointOnRightSideOfEdge2D(A, B, P, CLOSENESS_EPS = EPS):
	CCWABP = CCW2D(A, B, P)
	if CCWABP != 1 and CCWABP != 0:
		if CCWABP == -1:
			#Do a perpendicular projection onto the segment
			#to make sure it isn't a super close call
			vAB = B - A
			vAP = P - A
			proj = vAB.projPerp(vAP)
			if proj.squaredMag() < CLOSENESS_EPS:
				return True
			return False
		#Check endpoints
		elif CCWABP == -2:
			vPA = A - P
			if vPA.squaredMag() < CLOSENESS_EPS:
				return True
			return False
		elif CCWABP == 2:
			vPB = B - P
			if vPB.squaredMag() < CLOSENESS_EPS:
				return True
			return False
		else:
			print "ERROR in pointOnRightSideOfEdge2D: Shouldn't have gotten here"
	return True

#This is a helper function for "getCutsInsideTriangle()" in the Equidecomposability project
#and also a helper function for ear cutting triangulation
def pointInsideTriangle2D(A, B, C, P, CLOSENESS_EPS = EPS):
	isInside = pointOnRightSideOfEdge2D(A, B, P, CLOSENESS_EPS)
	isInside = isInside and (pointOnRightSideOfEdge2D(B, C, P, CLOSENESS_EPS))
	isInside = isInside and (pointOnRightSideOfEdge2D(C, A, P, CLOSENESS_EPS))
	return isInside

if __name__ ==  '__main__':
	A = Point3D(-1, 0, 0)
	B = Point3D(0, 1, 0)
	C = Point3D(1, 0, 0)
	P = Point3D(0, 1, 0)
	print pointInsideTriangle2D(A, B, C, P)
	print a
