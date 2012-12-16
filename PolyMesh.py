from Primitives3D import *
from Shapes3D import *
from OpenGL.GL import *
from Graphics3D import *
import sys
import re

class MeshVertex(object):
	def __init__(self, P, ID):
		self.pos = P
		self.ID = ID
		self.edges = set() #Store reference to all emanating edges
		#NOTE: Edges are not guaranteed to be in any particular order
		self.component = -1 #Which connected component it's in
	
	def getVertexNeighbors(self):
		ret = [0]*len(self.edges)
		i = 0
		for edge in self.edges:
			ret[i] = edge.vertexAcross(self)
			i = i+1
		return ret
	
	#Return a set of all faces attached to this vertex
	def getAttachedFaces(self):
		ret = set()
		i = 0
		for edge in self.edges:
			if (edge.f1 != None):
				ret.add(edge.f1)
			if (edge.f2 != None):
				ret.add(edge.f2)
		return ret
	
	#Get an estimate of the vertex normal by taking a weighted
	#average of normals of attached faces
	def getNormal(self):
		faces = self.getAttachedFaces()
		totalArea = 0.0;
		normal = Vector3D(0, 0, 0)
		for f in faces:
			w = f.getArea()
			totalArea = totalArea + w
			normal = normal + w*f.getNormal()
		return (1.0/totalArea)*normal
	
	#Sort the edges so that they are in CCW order when projected onto
	#the plane defined by the vertex's normal
	def getCCWSortedEdges(self):
		n = self.getNormal()
		#Store an extra field, vec, in each edge that represents
		#the normalized vector along the edge starting at this vertex
		e0 = None
		for e in self.edges:
			e.vec = e.vertexAcross(self).pos - self.pos
			e.vec = n.projPerp(e.vec)
			e.vec.normalize()
			e0 = e
		#Put all vectors in the xy plane
		v = e0.vec
		w = v%n
		B = Matrix4([v.x, v.y, v.z, 0, w.x, w.y, w.z, 0, n.x, n.y, n.z, 0, 0, 0, 0, 1])
		B = B.Inverse()
		for e in self.edges:
			e.vec = B*e.vec
		#Perform the sorting by CCW angle made with (1, 0) in the xy plane
		def getAngle(e):
			return math.atan2(e.vec.y, e.vec.x)
		return sorted(self.edges, key=getAngle)

class MeshFace(object):
	def __init__(self, ID):
		self.ID = ID
		self.edges = [] #Store edges in CCW order
		self.startV = 0 #Vertex that starts it off
		#Cache area, normal, and centroid
		self.area = None
		self.normal = None
		self.centroid = None

	def flipNormal(self):
		#Reverse the specification of the edges to make the normal
		#point in the opposite direction
		self.edges.reverse()
		self.normal = None

	#Return a list of vertices on the face in CCW order
	def getVertices(self):
		ret = [0]*len(self.edges)
		v = self.startV
		for i in range(0, len(self.edges)):
			ret[i] = v
			v = self.edges[i].vertexAcross(v)
		return ret
	
	def getNormal(self):
		if not self.normal:
			verts = [v.pos for v in self.getVertices()]
			self.normal = getFaceNormal(verts)
		return self.normal
	
	def getArea(self):
		if not self.area:
			verts = self.getVertices()
			self.area = getPolygonArea(verts)
		return self.area
	
	def getCentroid(self):
		if not self.centroid:
			ret = Point3D(0, 0, 0)
			verts = self.getVertices()
			if len(verts) == 0:
				return ret
			for v in verts:
				ret = ret + v.pos
			ret = ret*(1.0/float(len(verts)))
			self.centroid = ret
		return self.centroid
	
	def drawFilled(self, drawNormal = True):
		glBegin(GL_POLYGON)
		if drawNormal:
			normal = self.getNormal()
			if isinstance(normal, Vector3D):
				glNormal3f(normal.x, normal.y, normal.z)
		verts = self.getVertices()
		for v in verts:
			P = v.pos
			glVertex3f(P.x, P.y, P.z)
		glEnd()
	
	def drawBorder(self):
		glLineWidth(3)
		glBegin(GL_LINES)
		for e in self.edges:
			P1 = e.v1.pos
			P2 = e.v2.pos
			glVertex3f(P1.x, P1.y, P1.z)
			glVertex3f(P2.x, P2.y, P2.z)
		glEnd()
	
	def getPlane(self):
		return Plane3D(self.startV.pos, self.getNormal())

class MeshEdge(object):
	def __init__(self, v1, v2, ID):
		self.ID = ID
		[self.v1, self.v2] = [v1, v2]
		[self.f1, self.f2] = [None, None]
	
	def vertexAcross(self, startV):
		if startV == self.v1:
			return self.v2
		if startV == self.v2:
			return self.v1
		sys.stderr.write("Warning (vertexAcross): Vertex not member of edge\n")
		return None
	
	def addFace(self, face, v1):
		if self.f1 == None:
			self.f1 = face
		elif self.f2 == None:
			self.f2 = face
		else:
			sys.stderr.write("Cannot add face to edge; already 2 there\n")
	
	#Remove pointer to face
	def removeFace(self, face):
		if self.f1 == face:
			self.f1 = None
		elif self.f2 == face:
			self.f2 = None
		else:
			sys.stderr.write("Cannot remove edge pointer to face that was never part of edge\n")
	
	def faceAcross(self, startF):
		if startF == self.f1:
			return self.f2
		if startF == self.f2:
			return self.f1
		sys.stderr.write("Warning (faceAcross): Face not member of edge\n")
		return None
	
	def getCenter(self):
		[P1, P2] = [self.v1.pos, self.v2.pos]
		return 0.5*(P1 + P2)

def getFaceInCommon(e1, e2):
	e2faces = []
	if e2.f1 != None:
		e2faces.append(e2.f1)
	if e2.f2 != None:
		e2faces.append(e2.f2)
	if e1.f1 in e2faces:
		return e1.f1
	if e1.f2 in e2faces:
		return e1.f2
	return None

#############################################################
####                	POLYGON MESH                    #####
#############################################################

class PolyMesh(object):
	def __init__(self):
		self.DisplayList = -1
		self.needsDisplayUpdate = True
		self.drawEdges = 0
		self.drawVerts = 0
		self.drawNormals = 0
		self.vertices = []
		self.edges = []
		self.faces = []
		#Pointers to a representative vertex in different 
		#connected components
		self.components = []
	
	#Return the edge between v1 and v2 if it exists, or
	#return None if an edge has not yet been created 
	#between them
	def getEdge(self, v1, v2):
		edge = v1.edges & v2.edges
		if len(edge) > 1:
			sys.stderr.write("Warning: More than one edge found on vertex list intersection\n")
		for e in edge:
			return e
		return None

	#############################################################
	####                ADD/REMOVE METHODS                  #####
	#############################################################

	def addVertex(self, P):
		vertex = MeshVertex(P, len(self.vertices))
		self.vertices.append(vertex)
		return vertex
	
	#Create an edge between v1 and v2 and return it
	#This function assumes v1 and v2 are valid vertices in the mesh
	def addEdge(self, v1, v2):
		edge = MeshEdge(v1, v2, len(self.edges))
		self.edges.append(edge)
		v1.edges.add(edge)
		v2.edges.add(edge)
		return edge
	
	#Given a list of pointers to mesh vertices in CCW order
	#create a face object from them
	def addFace(self, meshVerts):
		verts = [v.pos for v in meshVerts]
		if not arePlanar(verts):
			sys.stderr.write("Error: Trying to add mesh face that is not planar\n")
			for v in verts:
				print v
			return None
		if not are2DConvex(verts):
			sys.stderr.write("Error: Trying to add mesh face that is not convex\n")
			return None
		face = MeshFace(len(self.faces))
		face.startV = meshVerts[0]
		for i in range(0, len(meshVerts)):
			v1 = meshVerts[i]
			v2 = meshVerts[(i+1)%len(meshVerts)]
			edge = self.getEdge(v1, v2)
			if edge == None:
				edge = self.addEdge(v1, v2)
			face.edges.append(edge)
			edge.addFace(face, v1) #Add pointer to face from edge
		self.faces.append(face)
		return face
	
	#Remove the face from the list of faces and remove the pointers
	#from all edges to this face
	def removeFace(self, face):
		#Swap the face to remove with the last face (O(1) removal)
		self.faces[face.ID] = self.faces[-1]
		self.faces[face.ID].ID = face.ID #Update ID of swapped face
		face.ID = -1
		self.faces.pop()
		#Remove pointers from all of the face's edges
		for edge in face.edges:
			edge.removeFace(face)
	
	#Remove this edge from the list of edges and remove 
	#references to the edge from both of its vertices
	#(NOTE: This function is not responsible for cleaning up
	#faces that may have used this edge; that is up to the client)
	def removeEdge(self, edge):
		#Swap the edge to remove with the last edge
		self.edges[edge.ID] = self.edges[-1]
		self.edges[edge.ID].ID = edge.ID #Update ID of swapped face
		edge.ID = -1
		self.edges.pop()
		#Remove pointers from the two vertices that make up this edge
		edge.v1.edges.remove(edge)
		edge.v2.edges.remove(edge)
	
	#Remove this vertex from the list of vertices
	#NOTE: This function is not responsible for cleaning up any of
	#the edges or faces that may have used this vertex
	def removeVertex(self, vertex):
		self.vertices[vertex.ID] = self.vertices[-1]
		self.vertices[vertex.ID].ID = vertex.ID
		vertex.ID = -1
		self.vertices.pop()
	
	#############################################################
	####    TOPOLOGY SUBDIVISION AND REMESHING METHODS      #####
	#############################################################
	
	#Split every face into K+1 faces by creating a new vertex
	#at the midpoint of every edge, removing the original face, 
	#creating a new face connecting all the new vertices, and
	#creating triangular faces to fill in the gaps
	def splitFaces(self):
		#First subdivide each edge
		for e in self.edges:
			P = e.getCenter()
			e.centerVertex = self.addVertex(P)
		faces = list(self.faces)
		#Now create the new faces
		for f in faces:
			self.removeFace(f)
			#Add the inner face
			fInnerVerts = [e.centerVertex for e in f.edges]
			self.addFace(fInnerVerts)
			#Add the triangles around the border
			innerVerts = [0]*len(f.edges)
			outerVerts = [0]*len(f.edges)
			outerV = f.startV
			for i in range(0, len(f.edges)):
				outerV = f.edges[i].vertexAcross(outerV)
				outerVerts[i] = outerV
				innerVerts[i] = f.edges[i].centerVertex
			for i in range(0, len(f.edges)):
				triVerts = [innerVerts[i], outerVerts[i], innerVerts[(i+1)%len(innerVerts)]]
				self.addFace(triVerts)
		#Remove all edges that were on the original faces
		for f in faces:
			for e in f.edges:
				if e.ID != -1: #If the edge has not already been removed
					self.removeEdge(e)
		self.needsDisplayUpdate = True
			
	
	#Split every face into N triangles by creating a vertex at
	#the centroid of each face
	def starRemeshFaces(self, faces):
		for f in faces:
			#TODO: Implement normal meshes this way (move centroid along normal)??
			centroidP = f.getCentroid()
			centroid = self.addVertex(centroidP)
			verts = f.getVertices()
			#Remove face and replace with N triangles
			self.removeFace(f)
			for i in range(0, len(verts)):
				v1 = verts[i]
				v2 = verts[(i+1)%len(verts)]
				newVerts = [v1, v2, centroid]
				newFace = self.addFace(newVerts)
		self.needsDisplayUpdate = True
	
	def starRemesh(self):
		#Copy over the face list since it's about to be modified
		faces = list(self.faces)
		self.starRemeshFaces(faces)
	
	#Triangulate all faces that are not triangular by using
	#the star scheme
	def starTriangulate(self):
		faces = []
		for f in self.faces:
			if len(f.edges) > 3:
				faces.append(f)
		self.starRemeshFaces(faces)

	#Divide every polygonal face with N sides into (N-2) triangles
	#This assumes the polygonal faces are convex
	def minTrianglesRemesh(self):
		faces = list(self.faces)
		for f in faces:
			verts = f.getVertices()
			if len(verts) <= 3:
				continue
			#Remove face and replace with (N-2) triangles
			self.removeFace(f)
			v0 = verts[0]
			for i in range(1, len(verts)-1):
				v1 = verts[i]
				v2 = verts[i+1]
				newFace = self.addFace([v0, v1, v2])
		self.needsDisplayUpdate = True

	#For every vertex, create a new vertex a parameter t [0-0.5] of
	#the way along each of its N attached edges, and then "chop off"
	#the pyramid whose base is formed b the new vertices and whose apex
	#is the original vertex, creating a new planar face to cover the hole
	def truncate(self, t):
		choppedVertices = list(self.vertices)
		choppedFaces = list(self.faces)
		for f in choppedFaces:
			f.vertexMap = {}
		#Create all of the new vertices that result from chopping off
		#the top of the pyramid, and make a list of new faces
		for v in choppedVertices:
			edges = v.getCCWSortedEdges()
			newVerts = [0]*len(edges)
			#First make the new vertices
			for i in range(0, len(edges)):
				edge = edges[i]
				v2 = edge.vertexAcross(v)
				newPos = v.pos + t*(v2.pos-v.pos)
				newVerts[i] = self.addVertex(newPos)
			#Now create the faces around the base of the pyramid
			for i in range(0, len(edges)):
				[e1, e2] = [edges[i], edges[(i+1)%len(edges)]]
				[v1, v2] = [newVerts[i], newVerts[(i+1)%len(newVerts)]]
				face = getFaceInCommon(e1, e2)
				if face == None:
					#If there is a hole in between the edges, skip this
					continue
				face.vertexMap[v] = [v2, v1] #CCW is the opposite way across an edge
			#Add the face that fills the chopped off region
			self.addFace(newVerts)
		#Now update all of the topology with the following steps:
		#1) Remove the outdated faces
		for f in choppedFaces:
			self.removeFace(f)
		#2) Remove the outdated edges coming out of the chopped vertices
		#and then remove the chopped vertices themselves
		for v in choppedVertices:
			for e in list(v.edges):
				self.removeEdge(e)
			self.removeVertex(v)
		#3) Add back the original faces around the base of each pyramid
		#with the appropriate substitutions
		for f in choppedFaces:
			newFace = []
			thisV = f.startV
			for i in range(0, len(f.edges)):
				newVerts = f.vertexMap[thisV]
				newFace.append(newVerts[0])
				newFace.append(newVerts[1])
				thisV = f.edges[i].vertexAcross(thisV)
			if len(newFace) > 0:
				self.addFace(newFace)
	
	def flipNormals(self):
		for f in self.faces:
			f.flipNormal()
	
	def getConnectedComponents(self):
		for v in self.vertices:
			if v.component == -1:
				self.components.append(v)
				stack = [v]
				while len(stack) > 0:
					vcurr = stack[-1]
					if vcurr.component == -1:
						vcurr.component = len(self.components)-1
						stack.pop()
						for vNeighbor in vcurr.getVertexNeighbors():
							stack.append(vNeighbor)
					else:
						stack.pop()
	
	def getConnectedComponentCounts(self):
		counts = [0]*len(self.components)
		for v in self.vertices:
			if v.component > -1:
				counts[v.component] = counts[v.component] + 1
		return counts
	
	def deleteAllButLargestConnectedComponent(self):
		if len(self.components) == 0:
			self.getConnectedComponents()
		counts = self.getConnectedComponentCounts()
		largestComponent = 0
		largestCount = 0
		for i in range(0, len(counts)):
			if counts[i] > largestCount:
				largestCount = counts[i]
				largestComponent = i
		facesToDel = []
		edgesToDel = []
		#Delete faces first, then edges, then vertices
		for f in self.faces:
			if f.startV.component != largestComponent:
				edgesToDel = edgesToDel + f.edges
				facesToDel.append(f)
		for f in facesToDel:
			self.removeFace(f)
		for e in edgesToDel:
			if e.ID != -1:
				self.removeEdge(e)
		verticesToDel = []
		for v in self.vertices:
			if v.component != largestComponent:
				verticesToDel.append(v)
		for v in verticesToDel:
			self.removeVertex(v)
		#Now update the connected components list
		if len(self.vertices) > 0:
			self.components = [self.vertices[0]]
			for v in self.vertices:
				v.component = 0
		self.needsDisplayUpdate = True
	
	def fillHoles(self):
		print "TODO"

	#############################################################
	####                 GEOMETRY METHODS                   #####
	#############################################################
	def getCentroid(self):
		center = Vector3D(0.0, 0.0, 0.0)
		for v in self.vertices:
			center = center + v.pos
		center = center * (1.0/float(len(self.vertices)))
		return center
	
	def getBBox(self):
		if len(self.vertices) == 0:
			return BBox3D(0, 0, 0, 0, 0, 0)
		P0 = self.vertices[0].pos
		bbox = BBox3D(P0.x, P0.x, P0.y, P0.y, P0.z, P0.z)
		for v in self.vertices:
			bbox.addPoint(v.pos)
		return bbox
	
	def sliceBelowPlane(self, plane):
		facesToDel = []
		edgesToDel = []
		verticesToDel = []
		facesToAdd = []
		newVertices = []
		for e in self.edges:
			#Keep track of edges which intersect the plane
			#and the vertex that represents that intersection
			e.centerVertex = None
		for f in self.faces:
			v1 = f.startV
			deleteFace = False
			splitFaceStartE = None
			splitFaceStartEIndex = -1
			splitFaceStartV = None
			splitFaceEndE = None
			for i in range(0, len(f.edges)):
				e = f.edges[i]
				v2 = e.vertexAcross(v1)
				distV1 = plane.distFromPlane(v1.pos)
				distV2 = plane.distFromPlane(v2.pos)
				#Mark all vertices below the plane for deletion
				#(this will count vertices more than once but 
				#because mesh is manifold the amortized size will
				#still be linear in the number of vertices)
				if distV1 < 0:
					verticesToDel.append(v1)
				#Mark all edges that are fully or partially below
				#the plane for deletion.  This list will hold
				#every such edge at most twice since 2 faces meet
				#at an edge in a manifold mesh
				if distV1 < 0 or distV2 < 0:
					edgesToDel.append(e)
					deleteFace = True
				#Edge goes from negative side to plus side of plane
				if distV1 < 0 and distV2 >= 0:
					if not e.centerVertex:
						line = Line3D(v1.pos, v2.pos - v1.pos)
						P = line.intersectPlane(plane)[1]
						if P:
							e.centerVertex = self.addVertex(P)
					if e.centerVertex:
						splitFaceStartEIndex = i
						splitFaceStartV = v1
						splitFaceStartE = e
				#Edge goes from plus side to negative side of plane
				if distV1 >= 0 and distV2 <= 0:
					if not e.centerVertex:
						line = Line3D(v1.pos, v2.pos - v1.pos)
						P = line.intersectPlane(plane)[1]
						if P:
							e.centerVertex = self.addVertex(P)
					if e.centerVertex:
						splitFaceEndE = e					
				v1 = v2
			if deleteFace:
				facesToDel.append(f)
			#Walk along the split part of the face on the positive
			#side of the plane
			if splitFaceStartE and splitFaceEndE:
				newFace = [splitFaceStartE.centerVertex]
				i = splitFaceStartEIndex
				e = splitFaceStartE
				v1 = splitFaceStartV
				while e != splitFaceEndE:
					v1 = e.vertexAcross(v1)
					newFace.append(v1)
					i = (i+1)%len(f.edges)
					e = f.edges[i]
				newFace.append(splitFaceEndE.centerVertex)
				facesToAdd.append(newFace)
		#First remove all faces that are no longer relevant
		for f in facesToDel:
			self.removeFace(f)	
		#Now remove edges that are no longer relevant
		for e in edgesToDel:
			if e.ID != -1:
				self.removeEdge(e)
		#Now remove vertices that are no longer relevant
		for v in verticesToDel:
			if v.ID != -1:
				self.removeVertex(v)
		#Add new faces
		for f in facesToAdd:
			self.addFace(f)
		self.needsDisplayUpdate = True
	
	#############################################################
	####                INPUT/OUTPUT METHODS                #####
	#############################################################
	def loadFile(self, filename):
		suffix = re.split("\.", filename)[-1]
		if suffix == "off":
			self.loadOffFile(filename)
		elif suffix == "obj":
			self.loadObjFile(filename)
		else:
			print "Unsupported file suffix (%s) for loading mesh"%(suffix, filename)
	
	def saveFile(self, filename, verbose = False):
		suffix = re.split("\.", filename)[-1]
		if suffix == "off":
			self.saveOffFile(filename, verbose)
		elif suffix == "obj":
			self.saveObjFile(filename, verbose)
		else:
			print "Unsupported file suffix (%s) for saving mesh %s"%(suffix, filename)		
	
	def loadOffFile(self, filename):
		fin = open(filename, 'r')
		nVertices = 0
		nFaces = 0
		nEdges = 0
		lineCount = 0
		face = 0
		vertex = 0
		for line in fin:
			lineCount = lineCount+1
			fields = line.split() #Splits whitespace by default
			if len(fields) == 0: #Blank line
				continue
			if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
				continue
			#Check section
			if nVertices == 0:
				if fields[0] == "OFF":
					if len(fields) > 2:
						fields[1:4] = [int(field) for field in fields]
						[nVertices, nFaces, nEdges] = fields[1:4]
				else:
					fields[0:3] = [int(field) for field in fields]
					[nVertices, nFaces, nEdges] = fields[0:3]
			elif vertex < nVertices:
				fields = [float(i) for i in fields]
				P = Point3D(fields[0], fields[1], fields[2])
				self.addVertex(P)
				vertex = vertex+1
			elif face < nFaces:
				#Assume the vertices are specified in CCW order
				fields = [int(i) for i in fields]
				meshVerts = fields[1:fields[0]+1]
				verts = [self.vertices[i] for i in meshVerts]
				self.addFace(verts)
				face = face+1
		fin.close()
			
	
	def saveOffFile(self, filename, verbose = False):
		nV = len(self.vertices)
		nE = len(self.edges)
		nF = len(self.faces)
		fout = open(filename, "w")
		fout.write("#Generated with Chris Tralie's G-RFLCT Library\n")
		fout.write("OFF\n%i %i %i\n"%(nV, nF, 0))
		for v in self.vertices:
			P = v.pos
			fout.write("%g %g %g\n"%(P.x, P.y, P.z))
		for f in self.faces:
			verts = f.getVertices()
			fout.write("%i "%(len(verts)))
			for v in verts:
				fout.write("%i "%(v.ID))
			fout.write("\n")
		fout.close()
		if verbose:
			print "Saved file to %s"%filename
		
	def loadObjFile(self, filename):
		#TODO: Right now vertex normals, face normals, and texture coordinates are ignored
		#Later incorporate them??
		fin = open(filename, 'r')
		for line in fin:
			fields = line.split()
			if len(fields) == 0: #Blank line
				continue
			if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
				continue
			if fields[0] == "v":
				coords = [float(i) for i in fields[1:4]]
				self.addVertex(Point3D(coords[0], coords[1], coords[2]))
			if fields[0] == "f":
				#Indices are numbered starting at 1 (so need to subtract that off)
				indices = [int(re.split("/",s)[0])-1 for s in fields[1:]]
				verts = [self.vertices[i] for i in indices]
				self.addFace(verts)
		fin.close()
	
	def saveObjFile(self, filename, verbose = False):
		fout = open(filename, "w")
		fout.write("#Generated with Chris Tralie's G-RFLCT Library\n")
		for v in self.vertices:
			P = v.pos
			fout.write("v %g %g %g\n"%(P.x, P.y, P.z))
		for f in self.faces:
			verts = f.getVertices()
			fout.write("f ")
			i = 0
			for i in range(0, len(verts)):
				v = verts[i]
				fout.write("%i"%(v.ID+1))#Indices are numbered starting at 1
				if i < len(verts) - 1:
					fout.write(" ")
			fout.write("\n")
		fout.close()
		if verbose:
			print "Saved file to %s"%filename
	
	def renderGL(self, drawEdges = 0, drawVerts = 0, drawNormals = 0):
		if self.drawEdges != drawEdges:
			self.drawEdges = drawEdges
			self.needsDisplayUpdate = True
		if self.drawVerts != drawVerts:
			self.drawVerts = drawVerts
			self.needsDisplayUpdate = True
		if self.drawNormals != drawNormals:
			self.drawNormals = drawNormals
			self.needsDisplayUpdate = True
		if self.needsDisplayUpdate:
			if self.DisplayList != -1: #Deallocate previous display list
				glDeleteLists(self.DisplayList, 1)
			self.DisplayList = glGenLists(1)
			glNewList(self.DisplayList, GL_COMPILE)
			for f in self.faces:
				f.drawFilled()
			if self.drawEdges:
				glDisable(GL_LIGHTING)
				glColor3f(0, 0, 1)
				glLineWidth(3)
				glBegin(GL_LINES)
				for e in self.edges:
					P1 = e.v1.pos
					P2 = e.v2.pos
					glVertex3f(P1.x, P1.y, P1.z)
					glVertex3f(P2.x, P2.y, P2.z)
				glEnd()
			if self.drawVerts:
				glDisable(GL_LIGHTING)
				glColor3f(1, 0, 0)
				glPointSize(5)
				glBegin(GL_POINTS)
				for v in self.vertices:
					P = v.pos
					glVertex3f(P.x, P.y, P.z)
				glEnd()
			if self.drawNormals:
				#For now just draw the vertex normals
				glDisable(GL_LIGHTING)
				glColor3f(0, 1, 0)
				glLineWidth(3)
				glBegin(GL_LINES)
				for v in self.vertices:
					P1 = v.pos
					P2 = P1 + 0.05*v.getNormal()
					glVertex3f(P1.x, P1.y, P1.z)
					glVertex3f(P2.x, P2.y, P2.z)
				glEnd()
			glEndList()
			self.needsDisplayUpdate = False
		glCallList(self.DisplayList)
	
	def renderCCWEdgesDebug(self):
		for v in self.vertices:
			edges = v.getCCWSortedEdges()
			for i in range(0, len(edges)):
				v1 = edges[i].vertexAcross(v)
				v2 = edges[(i+1)%len(edges)].vertexAcross(v)
				(P0, P1, P2) = (v.pos, v1.pos, v2.pos)
				normal = (P1 - P0) % (P2 - P0)
				normal.normalize()
				glBegin(GL_POLYGON)
				glNormal3f(normal.x, normal.y, normal.z)
				glVertex3f(P0.x, P0.y, P0.z)
				glVertex3f(P1.x, P1.y, P1.z)
				glVertex3f(P2.x, P2.y, P2.z)
				glEnd()
				Pc = P0 + P1 + P2
				Pc = (1.0/3.0)*Pc
				glDisable(GL_LIGHTING)
				glColor3f(1, 0, 0)
				glLineWidth(3)
				glBegin(GL_LINES)
				glVertex3f(Pc.x, Pc.y, Pc.z)
				Pc = Pc + normal
				glVertex3f(Pc.x, Pc.y, Pc.z)
				glEnd()
				glEnable(GL_LIGHTING)				
	
	#Slow version with no spatial subdivision
	def getRayIntersection(self, ray):
		t = float("inf")
		Point = None
		Face = None
		for f in self.faces:
			intersection = ray.intersectMeshFace(f)
			if intersection != None:
				if intersection[0] < t:
					t = intersection[0]
					Point = intersection[1]
					Face = f
		if isinstance(Point, Point3D):
			return [t, Point, Face]
		return None
	
	#Transformations are simple because geometry information is only
	#stored in the vertices
	def Transform(self, matrix):
		for v in self.vertices:
			v.pos = matrix*v.pos
	
	def __str__(self):
		nV = len(self.vertices)
		nE = len(self.edges)
		nF = len(self.faces)
		topology = nV-nE+nF
		return "PolyMesh Object: NVertices = %i, NEdges = %i, NFaces = %i, topology=%i"%(nV, nE, nF, topology)	

#Helper function for getBoxMesh and addFaceTiles
def makeBoxEdge(mesh, v1, v2, stepSize):
	if stepSize < 0:
		return [v1, v2]
	verts = [v1]
	direc = v2.pos - v1.pos
	frac = stepSize/direc.Length()
	#Round to the nearest integer number of tiles
	N = int(math.floor(1.0/frac+0.5))
	if N == 0:
		N = 1
	frac = 1.0/float(N)
	for i in range(1, N):
		newVert = mesh.addVertex(v1.pos+direc*frac*i)
		verts.append(newVert)
	verts.append(v2)
	return verts

#Helper function for getBoxMesh
def addFaceTiles(mesh, stepSize, ebott, eright, etop, eleft):
	topRow = etop
	index = 1
	for index in range(1, len(eleft)):
		bottomRow = None
		if index == len(eleft)-1:
			bottomRow = ebott
		else:
			bottomRow = makeBoxEdge(mesh, eleft[index], eright[index], stepSize)
		#Now add the square faces on this part
		for i in range(0, len(topRow)-1):
			mesh.addFace([bottomRow[i], bottomRow[i+1], topRow[i+1], topRow[i]])
		topRow = bottomRow

#L is length along z
#W is width along x
#H is height along y
#stepSize is the length of each square tile.  By default there are no tiles
#(stepSize = -1).  If one of the sides is not an integer multiple of the step size,
#then round to the nearest step size that would make it an integer multiple along
#that dimension
def getBoxMesh(L = 1.0, W = 1.0, H = 1.0, C = Point3D(0, 0, 0), stepSize = -1):
	mesh = PolyMesh()
	endpoints = []
	for dZ in [L/2.0, -L/2.0]:
		for dH in [-H/2.0, H/2.0]:
			for dW in [-W/2.0, W/2.0]:
				endpoints.append(mesh.addVertex(C+Point3D(dW, dH, dZ)))
	edgeIndices = [[0, 1], [1, 3], [3, 2], [2, 0], [1, 5], [5, 7], [7, 3], [7, 6], [6, 2], [0, 4], [4, 6], [4, 5]]
	edges = []
	edgesRev = []
	for edgePointers in edgeIndices:
		[v1, v2] = [endpoints[edgePointers[0]], endpoints[edgePointers[1]]]
		edges.append(makeBoxEdge(mesh, v1, v2, stepSize))
	for edge in edges:
		edgeRev = edge[:]
		edgeRev.reverse()
		edgesRev.append(edgeRev)
	#def addFaceTiles(mesh, stepSize, ebott, eright, etop, eleft):
	#Front Face
	addFaceTiles(mesh, stepSize, edges[0], edgesRev[1], edgesRev[2], edges[3])
	#Back Face
	addFaceTiles(mesh, stepSize, edgesRev[11], edgesRev[10], edges[7], edgesRev[5])
	#Left Face
	addFaceTiles(mesh, stepSize, edgesRev[9], edges[3], edges[8], edgesRev[10])
	#Right Face
	addFaceTiles(mesh, stepSize, edges[4], edgesRev[5], edgesRev[6], edgesRev[1])
	#Top Face
	addFaceTiles(mesh, stepSize, edgesRev[2], edges[6], edgesRev[7], edges[8])
	#Bottom Face
	addFaceTiles(mesh, stepSize, edges[11], edges[4], edges[0], edges[9])
	return mesh

def getRectMesh(P0, P1, P2, P3, stepSize = -1):
	mesh = PolyMesh()
	endpoints = [P0, P1, P2, P3]
	for i in range(0, len(endpoints)):
		endpoints[i] = mesh.addVertex(endpoints[i])
	edgeIndices = [[0, 1], [2, 1], [3, 2], [3, 0]]
	edges = []
	for edgePointers in edgeIndices:
		[v1, v2] = [endpoints[edgePointers[0]], endpoints[edgePointers[1]]]
		edges.append(makeBoxEdge(mesh, v1, v2, stepSize))
	addFaceTiles(mesh, stepSize, edges[0], edges[1], edges[2], edges[3])
	return mesh

if __name__ == '__main__':
	mesh = PolyMesh()
	mesh.loadFile("meshes/cube.obj")
	#mesh.removeFace(mesh.faces[10])
	#for face in mesh.faces:
	#s	print face," %i"%(face.ID)
	for v in mesh.vertices:
		print v.ID, " ",
	print "\n"
	mesh.removeVertex(mesh.vertices[2])	
	for v in mesh.vertices:
		print v.ID, " ",
