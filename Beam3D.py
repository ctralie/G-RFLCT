from Primitives3D import *
from Shapes3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from OpenGL.GL import *
import math
import pickle

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
	return 1e-7
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

#Each beam should have a coordinate system and a near distance
#The coordinate system should orient the "towards" direction
#perpendicular to the image plane of the beam
#"face" points to the MeshFace associated with the beam, if there is one
#Otherwise, the beam it is one of the initial cast beams and has not 
#reflected across a face yet. In that case the default neardist can be used
#Other member variables: 
#mvMatrix: Modelview matrix for putting a point into the coordinate system of a beam
#mvMatrixInverse: Takes a point from the beam's coordinates back into world coordinates
#frustPoints: Points of the frustum on the 2D image plane
#parent: The parent beam
#order: The "depth" of the beam; e.g. order 0 means this beam started at
#the origin and is heading towards its first boundary; order 1 means the beam
#has reflected off of one face already
#dummy: Used to construct the root of the beam tree
#freeBeam: True if this beam is not associated with a face, false otherwise
class Beam3D(object):
	def __init__(self, origin, frustVertices, parent = None, order = 0, face = None, dummy = False):
		self.dummy = dummy
		self.children = []
		self.parent = parent
		self.order = order
		self.origin = origin
		if dummy:
			return
		self.frustVertices = frustVertices
		self.neardist = 0.01
		self.face = face
		if len(frustVertices) < 3:
			print "ERROR: Only %i frustVertices on beam projection face"%len(frustVertices)
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
		self.mvMatrixInverse = self.mvMatrix.Inverse()
		#Now calculate the near distance if the beam is associated with a face
		self.freeBeam = True
		if face:
			P = self.mvMatrix*frustVertices[0]
			self.neardist = -P.z
			self.freeBeam = False
		#Now map the frustum points to 2D image plane coordinates
		self.frustPoints = self.projectPolygon(frustVertices, False)
	
	#Project a polygon onto the beam's image plane (before clipping to frustum)
	#Side effect: The field "origZ" is added to all point objects to store
	#the original z coordinates before projection
	def projectPolygon(self, polygon, doClip = True):
		#First transform all points into the beam's field of view
		mvVerts = [self.mvMatrix*v for v in polygon]
		#Now clip the polygon to the near plane
		clippedVerts = mvVerts
		clipDist = self.neardist
		if doClip:
			clippedVerts = []
			for i in range(0, len(mvVerts)):
				v1 = mvVerts[i]
				v2 = mvVerts[(i+1)%len(mvVerts)]
				(d1, d2) = (-v1.z, -v2.z)
				#v1 is behind the beam and v2 is in front of the beam
				if d1 < clipDist and d2 >= clipDist:
					ratio = (clipDist-d1)/(d2-d1)
					vNew = Point3D(v1.x+(v2.x-v1.x)*ratio, v1.y+(v2.y-v1.y)*ratio, -clipDist)
					clippedVerts.append(vNew)
				#v1 is in front of beam and v2 is behind beam
				elif d1 >= clipDist and d2 < clipDist:
					ratio = (clipDist-d2)/(d1 - d2)
					vNew = Point3D(v2.x+(v1.x-v2.x)*ratio, v2.y+(v1.y-v2.y)*ratio, -clipDist)
					clippedVerts.append(v1)
					clippedVerts.append(vNew)
				#Both vertices are in front of the beam
				elif d1 >= clipDist and d2 >= clipDist:
					#Just add the left one
					clippedVerts.append(v1)
				#Otherwise two vertices are behind
		#FOR DEBUGGING
		#for v in clippedVerts:
		#	print self.mvMatrixInverse*v
		#Now do the perspective projection
		for i in range(0, len(clippedVerts)):
			if clippedVerts[i].z >= 0:
				if doClip:
					print "ERROR: Near clipping did not work properly"
				else:
					print "WARNING: Beam focal point is in the plane of the defining face!!"
					clippedVerts[i].z = 1e-5
			#Store the original z position before projection onto image plane
			clippedVerts[i].x = -self.neardist*clippedVerts[i].x/clippedVerts[i].z
			clippedVerts[i].y = -self.neardist*clippedVerts[i].y/clippedVerts[i].z
			clippedVerts[i].z = -clipDist
		#Make sure the points are in counter clockwise order
		for i in range(0, len(clippedVerts)-2):
			P0 = clippedVerts[i]
			P1 = clippedVerts[i+1]
			P2 = clippedVerts[(i+2)%len(clippedVerts)]
			ccw = CCW2D(P0, P1, P2)
			if ccw == 1:
				clippedVerts.reverse()
				break
		return clippedVerts
	
	
	#Perform Sutherland Hodgman Clipping to clip a projected polygon to
	#the inside of the beam
	#Following pseudocode on http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
	#The function assumes that polygon2D has been clipped to the near plane
	#and put in image plane coordinates (with the use of the function projectPolygon)
	#Vertices that result from clipping are marked with the field "clippedVertex" as True
	#so that they can be distinguished from vertices that are unchanged
	#TODO: Make sure this function can handle a 1 Point polygon, so I can use that
	#to test whether a receiver position is within a beam
	def clipToFrustum(self, polygon2D):
		outputList = polygon2D[:]
		for v in outputList:
			v.clippedVertex = False
		for i in range(0, len(self.frustPoints)):
			if len(outputList) == 0: #Special case: No Points left
				break
			clipEdge = [self.frustPoints[i], self.frustPoints[(i+1)%len(self.frustPoints)]]
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
	
	def clipPolygon(self, poly):
		ret = self.projectPolygon(poly)
		return self.clipToFrustum(ret)
	
	def projectMeshFace(self, face):
		return self.projectPolygon([v.pos for v in face.getVertices()])
	
	def clipMeshFace(self, face):
		return self.clipPolygon([v.pos for v in face.getVertices()])
	
	#Put a face that has been clipped to the beam back into world coordinates
	#by casting rays through the clipped vertices in the image plane
	#and intersecting them with the face in world coordinates
	#clipped: a list of clipped points in 2D coordinates
	#face: a MeshFace object representing the face in world coordinates
	def projectBackClippedFace(self, clipped, face):
		#NOTE: Can do ray intersect plane instead of ray intersect face
		#because the rays are known to reside inside of the bounds of the 
		#face due to the nature of the clipping algorithm
		plane = face.getPlane()
		vertices = []
		for vOrig in clipped:
			v = vOrig.Copy()
			v = self.mvMatrixInverse*v
			line = Line3D(self.origin, v-self.origin)
			intersection = line.intersectPlane(plane)
			if intersection:
				vertices.append(intersection[1])
				if intersection[1].z > 100 and False:
					print line
					print plane
					print intersection[1]
					print v
					print vOrig
					print self.mvMatrixInverse
					print self.mvMatrix
					print self.origin
					print "\n\n"
			else:
				print "ERROR: Unable to find intersected point on face"
				print vOrig
				print line
		return vertices	
	
	#Purpose: Return the largest face that is completely unobstructed from
	#the point of view of this beam
	#Faces contains a list of MeshFace objects to be clipped against this beam
	#Returns a tuple (projected and clipped vertices, clipped vertices back on face plane, face object)
	#backProjectedFaces is a by-reference array for storing the clipped faces back
	#in world coordinates
	#Side effects: adds the field "visibleFaces" to the beam
	def findLargestUnobstructedFace(self, faces, backProjectedFaces = []):
		validFaces = []
		#First find the faces that are actually within the beam
		for face in faces:
			clipped = self.clipMeshFace(face)
			if len(clipped) > 2:
				#Only consider the face if it is within the beam
				validFaces.append((clipped, face))
		if len(validFaces) == 0:
			return None
		#Stores the visible faces to save work for split beams
		self.visibleFaces = [f[1] for f in validFaces]
		retFace = None
		faceArea = 0.0
		for i in range(0, len(validFaces)):
			face = validFaces[i]
			#Put the clipped coordinates of this face back into world
			#coordinates
			vertices = self.projectBackClippedFace(face[0], face[1])
			backProjectedFaces.append(vertices)
			#Check "face" against the plane of every other valid face
			faceInFront = True
			for j in range(0, len(validFaces)):
				if i == j:
					continue
				otherFace = validFaces[j]
				normal = otherFace[1].getNormal()
				P0 = otherFace[1].startV.pos
				dV = P0 - self.origin
				#Make sure plane normal is pointing in the direction of this beam
				if dV.Dot(normal) < 0:
					normal = (-1)*normal
				plane = Plane3D(P0, normal)
				#Check every clipped point of face against the plane of otherFace
				for v in vertices:
					if plane.distFromPlane(v) > EPS:
						#One of the points of this face is behind the plane
						#of another face
						faceInFront = False
						break
			if faceInFront:
				area = getPolygonArea(vertices)
				if area > faceArea:
					faceArea = area
					retFace = (face[0], vertices, face[1])
		if not retFace:
			print "Warning: %i faces visible in beam but no face found in front"%len(validFaces)
		return retFace
			
	def drawPolygon(self, poly, minX, maxX, minY, maxY, dim):
		border = 10
		scaleX = float(dim - border*2)/(maxX - minX)
		scaleY = float(dim - border*2)/(maxY - minY)
		glBegin(GL_POLYGON)
		for P in poly:
			glVertex2f(border+(P.x-minX)*scaleX, border+(P.y-minY)*scaleY)
		glEnd()
		glPointSize(6)
		glBegin(GL_POINTS)
		for P in poly:
			glVertex2f(border+(P.x-minX)*scaleX, border+(P.y-minY)*scaleY)
		glEnd()
	
	def drawProjectedMeshFaces(self, faces, dim, toggleDrawSplits):
		glLineWidth(2)
		glDisable(GL_LIGHTING)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(0, dim, dim, 0, 0, 1)
		glDisable(GL_DEPTH_TEST)
		glMatrixMode (GL_MODELVIEW)
		glLoadIdentity()
		glColor3f(1, 1, 1)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
		glBegin(GL_POLYGON)
		glVertex2f(0, 0)
		glVertex2f(dim, 0)
		glVertex2f(dim, dim)
		glVertex2f(0, dim)
		glEnd()
		#Calculate bounding box for beam
		minX = min([P.x for P in self.frustPoints])
		maxX = max([P.x for P in self.frustPoints])
		minY = min([P.y for P in self.frustPoints])
		maxY = max([P.y for P in self.frustPoints])
		glColor3f(1, 0, 0)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
		faceInFront = self.findLargestUnobstructedFace(faces)
		if faceInFront:
			self.drawPolygon(faceInFront[0], minX, maxX, minY, maxY, dim)
		glColor3f(0, 0, 0)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
		self.drawPolygon(self.frustPoints, minX, maxX, minY, maxY, dim)
		if toggleDrawSplits:
			if faceInFront:
				glColor3f(0, 1, 0)
				print "BEAM POINTS: "
				for P in self.frustPoints:
					print P
				print "FACEINFRONT POINTS:"
				for P in faceInFront[0]:
					print P
				subBeams = splitBeam(self, faceInFront[0])
				print "There are %i subbeams"%len(subBeams)
				for subBeam in subBeams:
					self.drawPolygon(subBeam, minX, maxX, minY, maxY, dim)
				fout = open('beam.dat', 'w')
				pickle.dump(self, fout)
				fout.close()
				fout = open('faceInFront.dat', 'w')
				pickle.dump(faceInFront, fout)
				fout.close()
		else:
			glColor3f(0, 0, 1)
			clippedFaces = [self.clipMeshFace(face) for face in faces]
			clippedFaces = [face for face in clippedFaces if len(face) > 2]
			for i in range(0, len(clippedFaces)):
				face = clippedFaces[i]
				val = float(i)/float(len(clippedFaces))
				glColor3f(val, val, val)
				glColor3f(0, 0, 0)
				self.drawPolygon(face, minX, maxX, minY, maxY, dim)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_LIGHTING)
		
	
	def __str__(self):
		ret = "Beam3D: origin = %s, neardist = %s, Points = "%(self.origin, self.neardist)
		for v in self.frustVertices:
			ret = ret + "%s "%v
		return ret

#split the beam around face into multiple convex regions
#beam and face are both assumed to be on the image plane so
#that convex splitting can occur in 2D
#it is also assumed that "face" has already been clipped to the beam
#The split beams are returned as a list of point lists, where
#the frustum points are still in the beam's image plane coordinates
def splitBeam(beam, face):
	newBeams = []
	beamPoints = [P.Copy() for P in beam.frustPoints]
	for i in range(0, len(face)):
		#Step 1: Figure out where the line constructed from
		#this face segment intersects the remaining part of the beam
		#The two points of "face" under consideration
		fP1 = face[i]
		fP2 = face[(i+1)%len(face)]
		faceLine = Line3D(fP1, fP2-fP1)
		#Left and right beam indexes are the starting indexes of
		#the beam segment that was intersected (i.e. the index before
		#the intersection)
		leftIntersection = None
		leftBeamIndex = -1 
		rightIntersection = None
		rightBeamIndex = -1
		for j in range(0, len(beamPoints)):
			bP1 = beamPoints[j]
			bP2 = beamPoints[(j+1)%len(beamPoints)]
			if CCW2D(bP1, bP2, fP1) == 0:
				leftIntersection = fP1
				leftBeamIndex = j
			if CCW2D(bP1, bP2, fP2) == 0:
				rightIntersection = fP2
				rightBeamIndex = j
			if leftIntersection and rightIntersection:
				break
			beamLine = Line3D(bP1, bP2-bP1)
			ret = beamLine.intersectOtherLineRet_t(faceLine)
			if ret:
				(t, intersection) = ret
				#If the intersection is actually within the segment
				if CCW2D(bP1, bP2, intersection) == 0:
					#Left intersection
					ccw = CCW2D(fP1, fP2, intersection)
					if ccw == -2:
						leftIntersection = intersection
						leftBeamIndex = j
					elif ccw == 2:
						rightIntersection = intersection
						rightBeamIndex = j
					#else:
					#	print "ERROR: Unexpected beam segment face segment intersection case, ccw = %i"%ccw
		#Step 2: Figure out the polygon that's formed between the line
		#from this face segment and the convex section of the beam it cuts off
		if leftBeamIndex == -1 or rightBeamIndex == -1:
			print "ERROR: Not all interesections were found while trying to split beam"
			print "fP1 = %s, fP2 = %s"%(fP1, fP2)
			if leftBeamIndex == -1:
				print "ERROR: No left intersection found while trying to split beam"
			if rightBeamIndex == -1:
				print "ERROR: No right intersection found while trying to split beam"
			continue
		if leftBeamIndex == rightBeamIndex:
			continue
		eps = 1e-4
		newBeam = [leftIntersection]
		index = (leftBeamIndex+1)%len(beamPoints)
		while True:
			if not PointsEqual2D(newBeam[-1], beamPoints[index], eps):
				newBeam.append(beamPoints[index])
			if index == rightBeamIndex:
				break
			index = (index+1)%len(beamPoints)
		if not PointsEqual2D(newBeam[-1], rightIntersection, eps):
			newBeam.append(rightIntersection)
		#This takes care of the corner case where only an edge or 
		#point is inside the beam
		if len(newBeam) >= 3:
			newBeams.append(newBeam)
		#Step 3: Cut the polygon we just found off of the remaining part of the beam
		#As long as the intersection was not one of the beam endpoints
		#add it to the appropriate place within the beam polygon
		cutBeamPoints = [rightIntersection]
		index = (rightBeamIndex+1)%len(beamPoints)
		while True:
			if not PointsEqual2D(cutBeamPoints[-1], beamPoints[index], eps):
				cutBeamPoints.append(beamPoints[index])
			if index == leftBeamIndex:
				break
			index = (index + 1)%len(beamPoints)
		if not PointsEqual2D(cutBeamPoints[-1], leftIntersection, eps):
			cutBeamPoints.append(leftIntersection)
		beamPoints = cutBeamPoints
	return newBeams
		

class BeamTree(object):
	#Construct all beams up to a maximum order of "maxOrder"
	#allFaces is a list of mesh face objects in the scene in world coordinates
	def __init__(self, origin, allFaces, maxOrder = 0):
		self.allFaces = allFaces
		self.origin = origin
		#"root" is where the beam tree starts; it is a dummy beam
		#There will be 6 roots of each beam, for the six faces of a cube
		#that encompasses the full surface area possible
		self.root = Beam3D(origin, [], None, 0, None, True)
		roots = []
		#Front Face
		verts = [v+origin for v in [Point3D(-1, -1, 1), Point3D(1, -1, 1), Point3D(1, 1, 1), Point3D(-1, 1, 1)]]
		roots.append(Beam3D(origin, verts, self.root))
		if False:
			#Back Face
			verts.reverse()
			verts = [v+Point3D(0, 0, -2) for v in verts]
			roots.append(Beam3D(origin, verts, self.root))
			#Left Face
			verts = [v+origin for v in [Point3D(-1, -1, -1), Point3D(-1, -1, 1), Point3D(-1, 1, 1), Point3D(-1, 1, -1)]]
			roots.append(Beam3D(origin, verts, self.root))
			#Right Face
			verts.reverse()
			verts = [v+Point3D(2, 0, 0) for v in verts]
			roots.append(Beam3D(origin, verts, self.root))
			#Top Face
			verts = [v+origin for v in [Point3D(-1, 1, 1), Point3D(1, 1, 1), Point3D(1, 1, -1), Point3D(-1, 1, -1)]]
			roots.append(Beam3D(origin, verts, self.root))
			#Bottom face
			verts.reverse()
			verts = [v+Point3D(0, -2, 0) for v in verts]
			roots.append(Beam3D(origin, verts, self.root))
		#Start the recursion on each one of these sub beams
		for beam in roots:
			self.root.children.append(beam)
			self.intersectBeam(beam, self.allFaces, maxOrder)
	
	#A recursive function that intersects "beam" with "faces" and splits/reflects the beam
	#until the maximum order is reached
	def intersectBeam(self, beam, faces, maxOrder):
		if beam.order == maxOrder:
			return
		if len(faces) == 0:
			return
		faceInFront = beam.findLargestUnobstructedFace(faces)
		if not faceInFront: #If there's just free space
			return
		
		#Split the beam around the visible face into sub-parts of the same
		#order and recursively split those sub-parts
		#NOTE: faceInFront[0] is the clipped and projected version of
		#the face in front so it should match 
		subBeams = splitBeam(beam, faceInFront[0])
		print "There are %i beams split around the front face"%len(subBeams)
		for subBeamPoints in subBeams:
			subBeamPoints = [beam.mvMatrixInverse*P for P in subBeamPoints]
			#print "[",
			#for P in subBeamPoints:
			#	print "%g,"%P.z,
			#print "]"
			subBeam = Beam3D(beam.origin, subBeamPoints, beam.parent, beam.order, beam.face)
			#Since this is a split and not a reflection the split part has the same
			#parent as "beam"
			beam.parent.children.append(subBeam)
			self.intersectBeam(subBeam, beam.visibleFaces, maxOrder)
		
		#Now Reflect the beam across this face and recursively intersect
		#that beam which has an order of +1
		#(this is second so that beam reflection occurs breadth first)
		#faceInFront is a tuple (projected and clipped vertices, clipped vertices back in 3D, face object)
		facePoints = faceInFront[1]
		face = faceInFront[2]
		#First generate split beam and add it to the list of children
		frontBeam = Beam3D(beam.origin, facePoints, beam.parent, beam.order, beam.face)
		beam.parent.children.append(frontBeam)
		#Next remove "beam" from its parent's children list since it's 
		#no longer relevant and has been split into a bunch of pieces
		#TODO: Make removing more efficient
		beam.parent.children.remove(beam)
		
		return
		#Now add the reflected beam as a child of the split front beam
		faceNormal = getFaceNormal(facePoints)
		dV = beam.origin - facePoints[0]
		perpFaceV = faceNormal.proj(dV)
		parFaceV = faceNormal.projPerp(dV)
		mirrorP0 = facePoints[0] + parFaceV - perpFaceV
		reflectedBeam = Beam3D(mirrorP0, facePoints, frontBeam, beam.order+1, face)
		frontBeam.children.append(reflectedBeam)
		self.intersectBeam(reflectedBeam, self.allFaces, maxOrder)

######################################################

#Testing a case where beam splitting did not work properly
if __name__ == '__main__':
	beamPoints = [Point3D(0.004, -0.004, -0.01), Point3D(0.01, -0.004, -0.01), Point3D(0.01, 0.01, -0.01), Point3D(0.004, 0.01, -0.01)]
	faceInFrontPoints2 = [Point3D(0.004, -0.004, -0.01), Point3D(0.01, -0.004, -0.01), Point3D(0.01, 0.01, -0.01), Point3D(0.004, 0.004, -0.01)]
	beam2 = Beam3D(Point3D(0, 0, 0), beamPoints)
	beam2.frustPoints = beamPoints
	print "SPLITTING TYPED BEAM"
	splitBeams2 = splitBeam(beam2, faceInFrontPoints2)
	print "DONE SPLITTING TYPED BEAM"
	fin = open("beam.dat", 'r')
	beam = pickle.load(fin)
	fin.close()
	fin = open("faceInFront.dat", 'r')
	faceInFront = pickle.load(fin)
	fin.close()
	print "SPLITTING LOADED BEAM"
	splitBeams = splitBeam(beam, faceInFront[0])
	print "FINISHED SPLITTING LOADED BEAM"
	print "DIFFERENCE BETWEEN LOADED AND TYPED:"
	for i in range(0, len(beam.frustPoints)):
		print beam.frustPoints[i] - beamPoints[i]
	for i in range(0, len(faceInFront[0])):
		print (faceInFront[0])[i] - faceInFrontPoints2[i]
	print "\n\nLOADED BEAM SPLIT:"
	for SB in splitBeams:
		for P in SB:
			print P
		print "\n"
	print "\n\nTYPED BEAM SPLIT:"
	for SB in splitBeams2:
		for P in SB:
			print P
		print "\n"	

if __name__ == '__main__6':
	fin = open("beam.dat", 'r')
	beam = pickle.load(fin)
	fin.close()
	fin = open("faceInFront.dat", 'r')
	faceInFront = pickle.load(fin)
	fin.close()
	[bP1, bP2, fP2] = [beam.frustPoints[0], beam.frustPoints[1], (faceInFront[0])[1]]
	print "bP1 = %s, bP2 = %s, fP2 = %s"%(bP1, bP2, fP2)
	print "%i"%CCW2D(bP1, bP2, fP2)

#Test a simple beam projection and split case
if __name__ == '__main__5':
	#Simple beam projection and splitting test for debugging
	origin = Point3D(0, 0, -2)
	frustVertices = [v+origin for v in [Point3D(-1, -1, 1), Point3D(1, -1, 1), Point3D(1, 1, 1), Point3D(-1, 1, 1)]]
	beam = Beam3D(origin, frustVertices)
	faceVertices = [Point3D(-2, -2, -4), Point3D(-2, -2, 4), Point3D(-2, 2, 4), Point3D(-2, 2, -4)]
	mesh = PolyMesh()
	meshVertices = []
	for V in faceVertices:
		meshVertices.append(mesh.addVertex(V))
	mesh.addFace(meshVertices)
	projected = beam.projectPolygon(faceVertices)
	print "PROJECTED POLYGON POINTS"
	for P in projected:
		print "%s"%(P)
	print "\nFRUSTUM POINTS:"
	for P in beam.frustPoints:
		print P
	clipped = beam.clipToFrustum(projected)
	print "\nCLIPPED POINTS"
	for P in clipped:
		print "%s"%(P)
	print "\nPOINTS PROJECTED BACK"
	projectedBack = beam.projectBackClippedFace(clipped, mesh.faces[0])
	for P in projectedBack:
		print P
		
#This main was used to debug numerical precision errors in Beam3D.clipToFrustum()
if __name__ == '__main__4':
	facePoints = [Point3D(0.8, -1.25, -1.3), Point3D(0.8, -1.25, -1.2), Point3D(0.8, -0.27, -1.2), Point3D(0.8, -0.27, -1.3)]
	origin = Point3D(0, 0, -2)
	beamVerts = [v+origin for v in [Point3D(-1, -1, 1), Point3D(1, -1, 1), Point3D(1, 1, 1), Point3D(-1, 1, 1)]]
	beam = Beam3D(origin, beamVerts)
	printMatlabPoly(beam.frustPoints, "beam")
	print "\n\n"
	projected = beam.projectPolygon(facePoints)
	printMatlabPoly(projected, "proj")
	print "\n\n"
	clipped = beam.clipPolygon(facePoints)
	printMatlabPoly(clipped, "clip")
	

if __name__== '__main__1':
	unitBox = getBoxMesh()
	tree = BeamTree(Point3D(0, 0, 0), unitBox.faces[:], 1)


if __name__ == '__main_2':
	origin =  Point3D(-2.5, 2.5, -2)
	frustPoints = [Point3D(1.25, 2.5, 2.5), Point3D(1.25, 2.5, -2.5), Point3D(3.75, 2.5, -2.5), Point3D(3.75, 2.5, 2.5)]
	points = [Point3D(-4.5, -3.75, 0), Point3D(0.5, -3.75, 0), Point3D(0.5, -6.25, 0), Point3D(-4.5, -6.25, 0)]
	beam = Beam3D(origin, frustPoints)
	beam.clipPolygon(points)

if __name__ == '__main__3':
	mesh = getRectMesh(Point3D(-1, -1, -1), Point3D(-1, 1, -1), Point3D(0, 1, -1), Point3D(0, -1, -1))
	beam = Beam3D(Point3D(0, 0, 0), [v.pos for v in mesh.faces[0].getVertices()], mesh.faces[0])
	#poly = [Vector3D(-3, 0, 0), Vector3D(0, 0, 1), Vector3D(2, 0, 1), Vector3D(3, 0, 0), Vector3D(7, 0, -4), Vector3D(5, 0, -5), Vector3D(3, 0, -5), Vector3D(0, 0, -3)]
	#poly = [Vector3D(-1, 0, -1)]#, Vector3D(-0.5, 0.5, -1), Vector3D(1.5, 0.5, -1), Vector3D(2, 0, -1), Vector3D(1.5, -0.5, -1), Vector3D(0, -1, -1)]
	#poly = [v+Vector3D(-100, 0, 0) for v in poly]
	poly = [Point3D(-1, -1, -1), Point3D(-1, 1, -1), Point3D(0, 1, -1), Point3D(0, -1, -2)]
	poly = beam.projectPolygon(poly)
	poly = beam.clipToFrustum(poly)
	for v in poly:
		print beam.mvMatrixInverse*v
