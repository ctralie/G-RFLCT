from Primitives3D import *
from Graphics3D import *
from Shapes3D import *
from PolyMesh import *
from OpenGL.GL import *
import math
import elementtree.ElementTree as ET
#TODO: Switch to lxml so I can get line numbers with my errors

class EMMaterial(object):
	def __init__(self, R = 1.0, T = 0.0):
		self.R = R;#Reflection coefficient
		self.T = T;#Transmission coefficient

class EMNode(object):
	#Make identity matrix by default
	#An optical material is simply used for visualizing in OpenGL
	#it probably won't have an effect on any RF simulations
	#If mesh is "None" then this node is a transformation node that takes us
	#one level further down in the tree
	#If 
	def __init__(self, parent = None, mesh = None, transformation = Matrix4(), EMMat = EMMaterial(), OpticalMat = OpticalMaterial()):
		self.parent = parent
		self.mesh = mesh
		self.children = []
		self.transformation = transformation
		self.EMMat = EMMat
		self.OpticalMat = OpticalMat
		#TODO: Add bounding boxes

class EMScene(object):
	def __init__(self):
		self.rootEMNode = EMNode()
		self.meshes = [] #Keep a list of meshes in the scene
		self.Source = None
		self.Receiver = None
		self.vSources = [] #Virtual sources
		self.paths = [] #Paths discovered from source to receiver
		self.rays = [] #For debugging

	def Read(self, filename):
		parentNode = self.rootEMNode
		self.ReadRecurse(filename, {}, {}, parentNode, None)
		self.getMeshList()
		self.buildVirtualSourceTree(1)
		self.getPathsToReceiver()

	def ReadRecurse(self, filename = '', EMMaterials = {}, OpticalMaterials = {},  EMParentNode = None, XMLNode = None):
		if (len(filename) > 0):
			tree = ET.parse(filename)
			XMLNode = tree.getroot()
			#Load in materials
			if XMLNode.find('OpticalMaterials') != None:
				node = XMLNode.find('OpticalMaterials')
				for mat in node.getchildren():
					name = mat.get("name")
					args = ["ka", "kd", "ks", "kt", "emission"]
					for arg in args:
						if mat.get(arg) == None:
							print "Error: No %s defined in optical material %s"%(arg, name)
							return
					ka = [float(i) for i in mat.get("ka").split()]
					kd = [float(i) for i in mat.get("ks").split()]
					ks = [float(i) for i in mat.get("ka").split()]
					kt = [float(i) for i in mat.get("kt").split()]
					emission = [float(i) for i in mat.get("emission").split()]
					OpticalMaterials[name] = OpticalMaterial(RGB3D(ka[0], ka[1], ka[2]), RGB3D(kd[0], kd[1], kd[2]), RGB3D(ks[0], ks[1], ks[2]), RGB3D(kt[0], kt[1], kt[2]), RGB3D(emission[0], emission[1], emission[2]))
			
			if XMLNode.find('EMMaterials') != None:
				node = XMLNode.find('EMMaterials')
				for mat in node.getchildren():
					name = mat.get("name")
					args = ["R", "T"]
					for arg in args:
						if mat.get(arg) == None:
							print "Error: No %s defined in EM material %s"%(arg, name)
							return
					R = float(mat.get("R"))
					T = float(mat.get("T"))
					EMMaterials[name] = EMMaterial(R, T)
		for currNode in XMLNode.getchildren():
			EMMat = EMMaterial()
			OpticalMat = OpticalMaterial()
			#Get materials (if specified...otherwise use defaults)
			if currNode.tag in ["tri", "box", "mesh"]:
				if currNode.get("em") != None:
					EMMat = EMMaterials[currNode.get("em")]
				if currNode.get("om") != None:
					OpticalMat = OpticalMaterials[currNode.get("om")]
			#Get transformation matrix (if it exists...otherwise use identity matrix)
			matrix = Matrix4()
			if currNode.tag in ["tri", "box", "mesh", "node", "include"]:
				if currNode.text != None:
					m = [float(i) for i in currNode.text.split()]
					if len(m) == 16:
						matrix = Matrix4(m)
					elif len(m) != 0:
						print "Invalid transformation matrix specification %s"%(currNode.text)
						return
			#Now process tag-specific elements	
			if currNode.tag == "tri":
				args = ["P0", "P1", "P2"]
				for arg in args:
					if currNode.get(arg) == None:
						print "Error: No %s defined in triangle"%(arg)
						return
				P0args = [float(i) for i in currNode.get("P0").split()]
				P1args = [float(i) for i in currNode.get("P1").split()]
				P2args = [float(i) for i in currNode.get("P2").split()]
				P0 = Point3D(P0args[0], P0args[1], P0args[2])
				P1 = Point3D(P1args[0], P1args[1], P1args[2])
				P2 = Point3D(P2args[0], P2args[1], P2args[2])
				mesh = PolyMesh()
				[V0, V1, V2] = [mesh.addVertex(P0), mesh.addVertex(P1), mesh.addVertex(P2)]
				mesh.addFace([V0, V1, V2])
				sceneNode = EMNode(EMParentNode, mesh, matrix, EMMat, OpticalMat)
				EMParentNode.children.append(sceneNode)
			elif currNode.tag == "box":
				args = ["length", "width", "height"]
				for arg in args:
					if currNode.get(arg) == None:
						print "Error: No %s defined in box"%(arg)
						return
				L = float(currNode.get("length")) #Length in z
				W = float(currNode.get("width")) #Width in x
				H = float(currNode.get("height")) #Height in y
				C = Point3D(0, 0, 0)
				#Center is an optional argument
				if currNode.get("center") != None:
					CPoints = [float(i) for i in currNode.get("center").split()] #Center of box
					C = Point3D(CPoints[0], CPoints[1], CPoints[2])
				mesh = getBoxMesh(L, W, H, C)
				sceneNode = EMNode(EMParentNode, mesh, matrix, EMMat, OpticalMat)
				EMParentNode.children.append(sceneNode)			
			elif currNode.tag == "mesh":
				args = ["filename"]
				for arg in args:
					if currNode.get(arg) == None:
						print "Error: No %s defined in mesh"%(arg)
						return
				meshfilename = currNode.get("filename")
				mesh = PolyMesh()
				mesh.loadFile(meshfilename)
				sceneNode = EMNode(EMParentNode, mesh, matrix, EMMat, OpticalMat)
				EMParentNode.children.append(sceneNode)
			elif currNode.tag == "node":
				#This is a transformation node in the graph
				sceneNode = EMNode(parent = EMParentNode, transformation = matrix)
				EMParentNode.children.append(sceneNode)
				#Now recursively add the branch in this node
				self.ReadRecurse('', EMMaterials, OpticalMaterials, sceneNode, currNode)
			elif currNode.tag == "include":
				#Recursively include another scene file as a branch in this tree
				if currNode.get("filename") == None:
					print "Error: No filename defined for included scene"
					return
				nextfilename = currNode.get("filename")
				sceneNode = EMNode(parent = EMParentNode)
				EMParentNode.children.append(sceneNode)
				self.ReadRecurse(nextfilename, {}, {}, matrix, sceneNode, None)
			elif currNode.tag == "Source":
				if currNode.get("pos") == None:
					print "Error: No position defined for EM Source"
					return
				coords = [float(i) for i in currNode.get("pos").split()]
				self.Source = Point3D(coords[0], coords[1], coords[2])
			elif currNode.tag == "Receiver":
				if currNode.get("pos") == None:
					print "Error: No position defined for EM Receiver"
					return
				coords = [float(i) for i in currNode.get("pos").split()]
				self.Receiver = Point3D(coords[0], coords[1], coords[2])
			else:
				if not (currNode.tag in ["EMMaterials", "OpticalMaterials"]):
					print "Unrecognized tag %s"%currNode.tag

	#Get a list of meshes in the scene
	#Transform all of the meshes into world coordinates but keep
	#track of the matrix that was used to do that
	def getMeshListRecurse(self, currEMNode, meshes, matrix):
		for child in currEMNode.children:
			transform = matrix*child.transformation
			if len(child.children) > 0:
				self.getMeshListRecurse(child, meshes, transform)
			if child.mesh != None:
				meshes.append(child.mesh)
				child.mesh.transform = transform
				child.mesh.Transform(transform)#Put mesh into world coordinates
				child.mesh.EMNode = child #Store this so I can lookup the material later

	def getMeshList(self):
		self.meshes = []
		self.getMeshListRecurse(self.rootEMNode, self.meshes, self.rootEMNode.transformation)
		

	def renderGLRecurse(self, currEMNode, matrix):
		for child in currEMNode.children:
			transform = matrix*child.transformation
			if len(child.children) > 0:
				self.renderGLRecurse(child, transform)
			if child.mesh != None:
				#Set optical material
				glEnable(GL_LIGHTING)
				ka = child.OpticalMat.ka
				kd = child.OpticalMat.kd
				ks = child.OpticalMat.ks
				emission = child.OpticalMat.emission
				glMaterialfv(GL_FRONT, GL_AMBIENT, [ka.R, ka.G, ka.B, ka.A])
				glMaterialfv(GL_FRONT, GL_DIFFUSE, [kd.R, kd.G, kd.B, kd.A])
				glMaterialfv(GL_FRONT, GL_SPECULAR, [ks.R, ks.G, ks.B, ks.A])
				#TODO: Include shininess in Optical Material Spec??
				glMaterialf(GL_FRONT, GL_SHININESS, 80)
				#Now draw
				#NOTE: I'm row major and GL is column major
				#so need to do transpose
				#NOTE: No longer need to do this step because the 
				#meshes have been pre-transformed
				#glPushMatrix()
				#glMultMatrixd(transform.Transpose().m)
				child.mesh.renderGL(drawEdges = 1)
				#glPopMatrix()
	
	def renderGL(self):
		self.renderGLRecurse(self.rootEMNode, self.rootEMNode.transformation)
	
	#TODO: Add better recursive intersect method with bounding box tests
	#NOTE: All meshes are in world coordinates now so no longer need to
	#transform the ray
	def getRayIntersection(self, ray):
		t = float("inf")
		Point = None
		face = None
		for m in self.meshes:
			thisRay = ray.Copy()
			intersection = m.getRayIntersection(thisRay)
			if intersection != None:
				thisPoint = intersection[1]
				thisFace = intersection[2]
				dVec = thisPoint - ray.P0
				this_t = dVec.Length()
				if this_t < t:
					t = this_t
					Point = thisPoint
					face = thisFace
		if isinstance(Point, Point3D):
			#Get the transformed face normal
			verts = [v.pos for v in face.getVertices()]
			normal = getFaceNormal(verts)
			#Make sure the normal is pointing in the right direction
			if normal.Dot(ray.V) > 0:
				normal = (-1)*normal;
			return (t, Point, normal, face)
		return None

	def buildVirtualSourceTreeRecurse(self, currNode, level, maxLevel):
		self.vSources.append(currNode)
		if level == maxLevel:
			return
		#Try to mirror this source around every face in the scene
		for m in self.meshes:
			EMNode = m.EMNode
			for f in m.faces:
				if f != currNode.meshFace:
					fP0 = currNode.pos
					facePoints = [P.pos for P in f.getVertices()]
					if level > 0: #Do pruning test to avoid making
						#virtual sources that cannot be reached from this one
						#meshFaceInFrustum(fP0, frustPoints, facePoints)
						frustPoints = [P.pos for P in currNode.meshFace.getVertices()]
						if not meshFaceInFrustum(fP0, frustPoints, facePoints):
							continue
					#Tests passed: need to make a virtual image here
					#First mirror the current virtual source point
					#across the plane of the new face
					faceNormal = getFaceNormal(facePoints)
					dV = fP0 - facePoints[0]
					perpFaceV = faceNormal.proj(dV)
					parFaceV = faceNormal.projPerp(dV)
					mirrorP0 = facePoints[0] + parFaceV - perpFaceV
					newSourceNode = EMVSourceNode(mirrorP0, currNode, f, m.EMNode)
					currNode.children.append(newSourceNode)
					self.buildVirtualSourceTreeRecurse(newSourceNode, level+1, maxLevel)

	def buildVirtualSourceTree(self, maxLevel):
		if not isinstance(self.Source, Point3D):
			print "Error: Trying to initialize virtual sources but no initial source specified"
			return
		self.vSources = []
		rootSource = EMVSourceNode(self.Source, None, None, None)
		self.buildVirtualSourceTreeRecurse(rootSource, 0, maxLevel)
	
	#This assumes the source tree has been built
	def getPathsToReceiver(self):
		self.paths = []
		for source in self.vSources:
			path = [self.Receiver]
			currSource = source
			validPath = True
			while currSource != None:
				target = path[-1]
				P0 = target
				V = currSource.pos - target
				ray = Ray3D(P0, V)
				intersection = self.getRayIntersection(ray)
				if intersection == None:
					validPath = False
					break
				#(t, Point, normal, face)
				if currSource.parent == None:
					#The last source
					if intersection[0] < (target-currSource.pos).Length() - EPS:
						#Something is blocking the final source
						validPath = False
					else:
						path.append(currSource.pos)
				else:
					if intersection[3] != currSource.meshFace:
						validPath = False
						break
					else:
						#Move the intersection point slightly away from the face
						path.append(intersection[1]+EPS*intersection[2])
				currSource = currSource.parent
			if validPath:
				self.paths.append(path)

#Used to build virtual source trees
class EMVSourceNode(object):
	#pos is the position of the virtual source in world coordinates
	#meshFace is the face around which this source was reflected
	def __init__(self, pos, parent, meshFace, EMNode):
		self.pos = pos
		self.parent = parent
		self.meshFace = meshFace
		self.EMNode = EMNode
		self.children = []
			
if __name__ == '__main__':
	scene = EMScene()
	scene.Read('test.xml')
