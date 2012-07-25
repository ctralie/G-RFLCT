from Primitives3D import *
from Shapes3D import *
from PolyMesh import *

class RGB3D(object):
	def __init__(self, R = 1.0, G = 1.0, B = 1.0, A = 0.0):
		[self.R, self.G, self.B, self.A] = [R, G, B, A]
	def clone(self):
		return RGB3D(self.R, self.G, self.B, self.A)
	def __add__(self, other):
		return RGB3D(self.R + other.R, self.G + other.G, self.B + other.B, self.A + other.A)
	def __mul__(self, a):
		return RGB3D(a*self.R, a*self.G, a*self.B, a*self.A)
	def __rmul__(self, a):
		return RGB3D(a*self.R, a*self.G, a*self.B, a*self.A)
	def Scale(self, other):
		self.R = self.R*other.R
		self.G = self.G*other.G
		self.B = self.B*other.B
		self.A = self.A*other.A
	def squaredMag(self):
		return self.R*self.R + self.G*self.G + self.B*self.B

class OpticalMaterial(object):
	def __init__(self, ka = RGB3D(0.2, 0.2, 0.2), kd = RGB3D(0.5, 0.5, 0.5), ks = RGB3D(0.5, 0.5, 0.5), kt = RGB3D(0.0, 0.0, 0.0), emission = RGB3D(0.0, 0.0, 0.0)):
		[self.ka, self.kd, self.ks] = [ka, kd, ks]
		[self.kt, self.emission] = [kt, emission]

class RadiosityMaterial(object):
	#em is the emission (RGB per unit area)
	#p is the reflectance (RGB)
	#LIncident is the incident light gathered in at one iteration (absolute)
	#BExcident is the best estimate of the radiosity so far (light per unit area)
	#BUnshot is the unshot radiosity used for progressive radiosity (unshot light per unit area)
	#NOTE: All radiosity values are per unit area
	def __init__(self, em = RGB3D(0, 0, 0), p = RGB3D(0.5, 0.5, 0.5)):
		self.em = em
		self.p = p
		self.LIncident = RGB3D(0, 0, 0)
		self.BExcident = em
		self.BUnshot = em.clone()
	
	#This function is important for when the material gets subdivided
	def clone(self):
		ret = RadiosityMaterial(self.em, self.p)
		ret.LIncident = self.LIncident.clone()
		ret.BExcident = self.BExcident.clone()
		ret.BUnshot = self.BUnshot.clone()
		return ret
