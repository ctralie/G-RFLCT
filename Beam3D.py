from Primitives3D import *
from Utilities2D import *
from Shapes3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from OpenGL.GL import *
import numpy as np
import math
import pickle


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
#face: The face of the image plane (None if this is a root beam)
#dummy: Used to construct the root of the beam tree
#terminalFace: The face on which this beam terminates
class Beam3D(object):
	def initMVMatrix(self, origin, points):
		self.towards = getFaceNormal(points)
		dV = points[0] - origin
		if dV.Dot(self.towards) < 0:
			self.towards = -1*self.towards
		self.right = points[1] - points[0]
		self.right.normalize()
		self.up = self.right%self.towards
		self.up.normalize()
		#Now calculate and store the modelview matrix to get into the 
		#beam's coordinate system
		#(NOTE this function follows the OpenGL convention that things in
		#front of the beam are -z)
		self.mvMatrix = getCameraMatrix(self.towards, self.up, self.right, self.origin)
		self.mvMatrixInverse = self.mvMatrix.Inverse()

	def __init__(self, origin, frustVertices, parent = None, order = 0, face = None, dummy = False, terminalFace = None):
		self.dummy = dummy
		self.children = []
		self.parent = parent
		self.order = order
		self.origin = origin
		if dummy:
			return
		self.frustVertices = frustVertices
		self.face = face
		self.terminalFace = terminalFace
		if face:
			#If the image plane is a face
			points = [V.pos for V in face.getVertices()]
			if len(points) < 3:
				print "ERROR: Only %i points on beam projection face"%len(points)
				return
			self.initMVMatrix(origin, points)
			P = self.mvMatrix*points[0]
			self.neardist = -P.z
		else:
			#If the beam is a root beam and there is no associated face
			if len(frustVertices) < 3:
				print "ERROR: Only %i frustVertices on beam projection face"%len(frustVertices)
				return
			self.initMVMatrix(origin, frustVertices)
			self.neardist = 0.01
		#Now map the frustum points to 2D image plane coordinates
		self.frustPoints = self.projectPolygon(frustVertices, False)
	
	#Project a polygon onto the beam's image plane (before clipping to frustum)
	#Also clip the polygon to the beam's image plane
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
	#The function assumes that polygon2D has been clipped to the near plane
	#and put in image plane coordinates (with the use of the function projectPolygon)
	#Vertices that result from clipping are marked with the field "clippedVertex" as True
	#so that they can be distinguished from vertices that are unchanged
	#TODO: Make sure this function can handle a 1 Point polygon, so I can use that
	#to test whether a receiver position is within a beam
	def clipToFrustum(self, polygon2D):
		return clipSutherlandHodgman(self.frustPoints, polygon2D)
	
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
			if face is self.face:
				#Don't consider the imaging plane
				continue
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
						print "plane.distFromPlane(v) = %g"%plane.distFromPlane(v)
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
				subBeams = splitBeam(self, faceInFront[0])
				for subBeam in subBeams:
					self.drawPolygon(subBeam, minX, maxX, minY, maxY, dim)
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
	
	#Draw the beam in space as a transparent object
	def drawBeam(self, color = None, beamTrans = 0.3):
		colormap = [(0, 1, 0), (1, 1, 0), (0, 1, 1), (0, 0, 1)]
		glDisable(GL_LIGHTING)
		glEnable(GL_BLEND)
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
		
		c = color
		#If the color is not specified, use a colormap which color codes the 
		#sections of the beam based on what order they are
		if not c:
			c = colormap[self.order%len(colormap)]
		P0 = self.origin
		Points = self.frustVertices

		if self.face:
			#The beam is a reflected beam and starts at a face
			#PointsStart are the points on the imaging face, and Points
			#are the points at the end of the beam
			PointsStart = [self.mvMatrixInverse*P for P in self.frustPoints]
			glColor4f(c[0], c[1], c[2], beamTrans)
			for i1 in range(0, len(Points)):
				glBegin(GL_POLYGON)
				i2 = (i1+1)%len(Points)
				toDraw = [PointsStart[i1], Points[i1], Points[i2], PointsStart[i2]]
				for P in toDraw:
					glVertex3f(P.x, P.y, P.z)
				glEnd()
			glColor3f(0, 0, 0.3)
			glLineWidth(5)
			glBegin(GL_LINES)
			for i in range(0, len(Points)):
				[P1, P2] = [PointsStart[i], Points[i]]
				glVertex3f(P1.x, P1.y, P1.z)
				glVertex3f(P2.x, P2.y, P2.z)
			glEnd()
		else:
			#The beam is a root beam and starts at a point
			glColor4f(c[0], c[1], c[2], beamTrans)
			glBegin(GL_TRIANGLES)
			for i in range(0, len(Points)):
				P1 = Points[i]
				P2 = Points[(i+1)%len(Points)]
				glVertex3f(P0.x, P0.y, P0.z)
				glVertex3f(P1.x, P1.y, P1.z)
				glVertex3f(P2.x, P2.y, P2.z)
			glEnd()
			glColor3f(0, 0, 0.3)
			glLineWidth(5)
			glBegin(GL_LINES)
			for P in Points:
				glVertex3f(P0.x, P0.y, P0.z)
				glVertex3f(P.x, P.y, P.z)
			glEnd()

		glDisable(GL_BLEND)
		glEnable(GL_LIGHTING)
		if self.order > 0:
			self.parent.drawBeam(color, beamTrans)
	
	#Draw outlines of shapes clipped to this beam in 3D
	def drawBackProjected(self, faces):
		glDisable(GL_LIGHTING)
		backProjectedFaces = []
		self.findLargestUnobstructedFace(faces, backProjectedFaces)
		glColor3f(1, 0, 0)
		glLineWidth(3)
		for face in backProjectedFaces:
			glBegin(GL_LINES)
			for i in range(0, len(face)):
				P1 = face[i]
				P2 = face[(i+1)%len(face)]
				glVertex3f(P1.x, P1.y, P1.z)
				glVertex3f(P2.x, P2.y, P2.z)
			glEnd()
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
	#"beams" is a list of lists of vertices which represent beams that start
	#at the origin.  If it is None,
	#make an omnidirectional set of beams that start at the beam origin by
	#making the beams go out each of 6 faces of a cube
	def __init__(self, origin, allFaces, maxOrder = 0, beams = None):
		self.allFaces = allFaces
		self.origin = origin
		#"root" is where the beam tree starts; it is a dummy beam
		#There will be 6 roots of each beam, for the six faces of a cube
		#that encompasses the full surface area possible
		self.root = Beam3D(origin, [], None, 0, None, True)
		roots = []
		if beams:
			for verts in beams:
				roots.append(Beam3D(origin, verts, self.root))
		else:
			#omnidirectional
			#Front Face
			verts = [v+origin for v in [Point3D(-1, -1, 1), Point3D(1, -1, 1), Point3D(1, 1, 1), Point3D(-1, 1, 1)]]
			roots.append(Beam3D(origin, verts, self.root))
			if False:#For now just look at front face for debugging
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
		if beam.order > maxOrder:
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
		#print "There are %i beams split around the front face"%len(subBeams)
		for subBeamPoints in subBeams:
			subBeamPoints = [beam.mvMatrixInverse*P for P in subBeamPoints]
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
		frontBeam.terminalFace = face
		beam.parent.children.append(frontBeam)
		#Next remove "beam" from its parent's children list since it's 
		#no longer relevant and has been split into a bunch of pieces
		#TODO: Make removing more efficient
		beam.parent.children.remove(beam)
		
		#Now add the reflected beam as a child of the split front beam
		faceNormal = getFaceNormal(facePoints)
		dV = beam.origin - facePoints[0]
		perpFaceV = faceNormal.proj(dV)
		parFaceV = faceNormal.projPerp(dV)
		mirrorP0 = facePoints[0] + parFaceV - perpFaceV
		reflectedBeam = Beam3D(mirrorP0, facePoints, frontBeam, beam.order+1, face)
		frontBeam.children.append(reflectedBeam)
		self.intersectBeam(reflectedBeam, self.allFaces, maxOrder)
	
	#faces: The faces this beam should intersect in order
	#If an element of faces is "None", then treat that as a wildcard
	def getBeamsIntersectingFaces(self, faces):
		matchedBeams = []
		beamsToCheck = self.root.children[:]
		while len(beamsToCheck) > 0:
			thisbeam = beamsToCheck.pop()
			beamIsValid = False
			if not faces[thisbeam.order]:
				beamIsValid = True
			elif faces[thisbeam.order] == thisbeam.terminalFace:
				beamIsValid = True
			if beamIsValid:
				if thisbeam.order == len(faces) - 1:
					matchedBeams.append(thisbeam)
				else:
					beamsToCheck = beamsToCheck + thisbeam.children
		return matchedBeams

	#Return an array of arrays of points that represent the faces
	#of the beam path 
	def extractBeamToOrigin(self, beam):
		print "TODO"

	#"point" is a point on a face of the beam.  This function will
	#return a set of line segments that represent a specular path
	#from this point on the terminal face of the beam back to the
	#origin
	def extractPathToOrigin(self, beam, point):
		#First generate all virtual images
		beamParts = [beam]
		while beamParts[-1].order > 0:
			beamParts.append(beamParts[-1].parent)
		beamParts.reverse()
		images = [None]*(beam.order+1)
		images[0] = self.origin
		for i in range(1, len(images)):
			P0 = images[i-1]
			face = beamParts[i-1].terminalFace
			faceNormal = face.getNormal()
			facePoint = face.startV.pos
			dV = P0 - facePoint
			perpFaceV = faceNormal.proj(dV)
			parFaceV = faceNormal.projPerp(dV)
			mirrorP0 = facePoint + parFaceV - perpFaceV
			images[i] = mirrorP0
		images = images[1:]
		path = [point]
		i = len(images) - 1
		while i >= 0:
			V = path[-1] - images[i]
			ray = Ray3D(images[i], V)
			plane = beamParts[i].terminalFace.getPlane()
			[t, P] = ray.intersectPlane(plane)
			path.append(P)
			i = i - 1
		path.reverse()
		path = [self.origin] + path
		return path
	
	#Return the image of the source of the beam tree reflected across
	#all faces until it gets to "beam"'s face
	def getImageSource(self, beam):
		#First generate all virtual images
		beamParts = [beam]
		while beamParts[-1].order > 0:
			beamParts.append(beamParts[-1].parent)
		beamParts.reverse()
		images = [None]*(beam.order+1)
		images[0] = self.origin
		for i in range(1, len(images)):
			P0 = images[i-1]
			face = beamParts[i-1].terminalFace
			faceNormal = face.getNormal()
			facePoint = face.startV.pos
			dV = P0 - facePoint
			perpFaceV = faceNormal.proj(dV)
			parFaceV = faceNormal.projPerp(dV)
			mirrorP0 = facePoint + parFaceV - perpFaceV
			images[i] = mirrorP0
		return images[-1]

	#Return an image that holds the interference pattern on the face "face"
	#due to beams
	#resx and resy are the number of samples in x and y, respectively, to take
	#freq is the frequency of the propagating wave (used to compute phase)
	def getInterferencePatternOnFace(self, face, beams, resx, resy, freq):
		wavelen = 3.0e8/freq
		#First construct a modelview matrix which transforms points on the plane
		#of the face into a coordinate system where Z = 0 on the plane
		towards = face.getNormal()
		verts = [v.pos for v in face.getVertices()]
		right = verts[1] - verts[0]
		towards.normalize()
		right.normalize()
		up = right % towards
		mvMatrix = getCameraMatrix(towards, up, right, verts[0])
		mvMatrixInverse = mvMatrix.Inverse()
		#Now calculate a 2D bounding box based on the vertices of the face
		#which will determine the area where the pattern is calculated
		facePoints2D = [mvMatrix*P for P in verts]
		for P in facePoints2D:
			print P
		bbox = BBox3D(0, 0, 0, 0, 0, 0)
		for P in facePoints2D:
			bbox.addPoint(P)
		if bbox.zmin < -EPS or bbox.zmax > EPS:
			print "ERROR: Plane transformation for interference pattern zmin = %g, zmax = %g"%(bbox.zmin, bbox.zmax)
		print bbox
		xstart = float(bbox.xmin)
		xend = float(bbox.xmax)
		ystart = float(bbox.ymin)
		yend = float(bbox.ymax)
		dx = (xend - xstart) / (resx - 1)
		dy = (yend - ystart) / (resy - 1)
		#Transform each of the beams' terminating face points into the plane coordinate system
		beamPoints2D = [[None]]*len(beams)
		for k in range(0, len(beams)):
			beamPoints2D[k] = [mvMatrix*P for P in beams[k].frustVertices]
			for i in range(0, len(beamPoints2D)-2):
				P0 = (beamPoints2D[k])[i]
				P1 = (beamPoints2D[k])[i+1]
				P2 = (beamPoints2D[k])[(i+2)%len(beamPoints2D)]
				ccw = CCW2D(P0, P1, P2)
				if ccw == 1:
					beamPoints2D[k].reverse()
					break
		#Precompute the image sources for each face
		imageSources = [self.getImageSource(beam) for beam in beams]
		#Now loop through and fill in the pattern
		pattern = np.zeros((resy, resx))
		for yi in range(0, resy):
			y = ystart + yi*dy
			for xi in range(0, resx):
				#TODO: Implement phasor addition from antenna response here
				RealResp = 0
				ImagResp = 0
				x = xstart + xi*dx
				P = Point3D(x, y, 0)
				PWorld = mvMatrixInverse*P
				#For every beam, check to see if the point (x, y) is inside of
				#that beam's terminating face
				for k in range(0, len(beams)):
					beam = beams[k]
					imageSource = imageSources[k]
					points = beamPoints2D[k]
					inside = True
					for i in range(0, len(points)):
						P0 = points[i]
						P1 = points[(i+1)%len(points)]
						if CCW2D(P0, P1, P) == 1:
							inside = False
							break
					if inside:
						#TODO: Maybe the exact paths will be more important later
						#instead of just the virtual image paths?
						dV = PWorld - imageSource
						pathLength = dV.Length()
						theta = 2*math.pi*pathLength/wavelen
						RealResp = RealResp + math.cos(theta)
						ImagResp = ImagResp + math.sin(theta)
				phase = math.atan2(ImagResp, RealResp)
				pattern[yi, xi] = phase
		return pattern

######################################################
