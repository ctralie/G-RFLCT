from Primitives3D import *
from Graphics3D import *
from Shapes3D import *
import math
import elementtree.ElementTree as ET

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

class EMScene(object):
	def EMScene(self):
		self.root = EMNode

	def Read(self, filename = '', EMMaterials = {}, OpticalMaterials = {},  EMNode = self.root, XMLNode = None):
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
			if currNode.tag == "tri":
				args = ["om", "em", "P0", "P1", "P2"]
				for arg in args:
					if currNode.get(arg) == None:
						print "Error: No %s defined in triangle"%(arg)
						return
				#Get material
				EMMat = EMMaterials[currNode.get("om")]
				OpticalMat = EMMaterials[currNode.get("em")]
				P0 = [float(i) for i in currNode.get("P0").split()]
				P1 = [float(i) for i in currNode.get("P1").split()]
				P2 = [float(i) for i in currNode.get("P2").split()]

