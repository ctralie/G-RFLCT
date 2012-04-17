from Primitives3D import *
from Shapes3D import *
from PolyMesh import *

class RGB3D(object):
	def __init__(self, R = 1.0, G = 1.0, B = 1.0, A = 0.0):
		[self.R, self.G, self.B, self.A] = [R, G, B, A]

class OpticalMaterial(object):
	def __init__(self, ka = RGB3D(0.2, 0.2, 0.2), kd = RGB3D(0.5, 0.5, 0.5), ks = RGB3D(0.5, 0.5, 0.5), kt = RGB3D(0.0, 0.0, 0.0), emission = RGB3D(0.0, 0.0, 0.0)):
		[self.ka, self.kd, self.ks] = [ka, kd, ks]
		[self.kt, self.emission] = [kt, emission]

#fP0 is the apex of the frustum
#frustPoints are the points that define the frustum face 
#(assumed to be convex and specified in order)
#facePoints are the points that define the face that's being checked
#for being inside the frustum (assumed to be convex and specified in order)
def meshFaceInFrustum(fP0, frustPoints, facePoints):
	if len(frustPoints) < 3:
		print "Error: Less than three points in face defining frustum"
		return False
	if len(facePoints) < 3:
		print "Error: Less than three points in face that's being checked against frustum"
		return False
	#Calculate the frustum normal (pointing away from fP0)
	frustNormal = getFaceNormal(frustPoints)
	if frustNormal.Dot(frustPoints[0] - fP0) < 0:
		frustNormal = (-1)*frustNormal
	#Calculate the face normal (pointing away from P0)
	faceNormal = getFaceNormal(facePoints)
	if faceNormal.Dot(frustPoints[0] - fP0) < 0:
		faceNormal = (-1)*faceNormal
	facePlane = Plane3D(facePoints[0], faceNormal)
	#Check that at least part of the face is in front of the frustum face
	allBehind = True
	for P in facePoints:
		if frustNormal.Dot(P - frustPoints[0]) > 0:
			allBehind = False
			break
	if allBehind:
		return False
	#Find the polygon frustum cross section that is
	#in the plane of the face by intersecting a ray through each
	#point on the frustum face with the face plane
	frustSlicePoints = []
	for P in frustPoints:
		ray = Ray3D(fP0, P-fP0)
		intersection = ray.intersectPlane(facePlane)
		if intersection == None: #TODO: Need to handle this case better!
		#(this is the case where projection of the face's plane to the frustum may be unbounded)
			return True
		if intersection != None:
			frustSlicePoints.append(intersection[1])
	#Find the minimum enclosing circle for the plane cross section
	#and for the face
	faceCircle = getMinimumEnclosingCircle(facePoints, faceNormal)
	frustSliceCircle = getMinimumEnclosingCircle(frustSlicePoints, faceNormal)
	if faceCircle == None:
		print "Warning: No minimum enclosing circle found for face points"
		print faceNormal
		for P in facePoints:
			print P
		print "\n\n"
		return True #Return true to be safe
	if frustSliceCircle == None:
		print "Warning: No minimum enclosing circle found for frustum slice"
		for P in frustSlicePoints:
			print P
		print "\n\n"
		return True
	return faceCircle.intersectsOtherCircle(frustSliceCircle)
