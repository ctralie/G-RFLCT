from Primitives3D import *
from Shapes3D import *
from PolyMesh import *

def splitIntoRGBA(val):
	A = (0xff000000&val)>>24
	R = (0x00ff0000&val)>>16
	G = (0x0000ff00&val)>>8
	B = (0x000000ff&val)
	return [R, G, B, A]

def extractFromRGBA(R, G, B, A):
	return ((A<<24)&0xff000000) | ((R<<16)&0x00ff0000) | ((G<<8)&0x0000ff00) | (B&0x000000ff)

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

if __name__ == '__main__':
	val = 0xdeefacc
	[R, G, B, A] = splitIntoRGBA(val)
	val2 = extractFromRGBA(R, G, B, A)
	print "%x, %x"%(val, val)
	print "[R, G, B, A] = [%u, %u, %u, %u]"%(R, G, B, A)
