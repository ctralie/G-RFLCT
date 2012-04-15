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

class Ray3D(object):
	def __init__(self, P0, V):
		self.P0 = P0
		self.V = V
	
	def Copy(self):
		return Ray3D(self.P0.Copy(), self.V.Copy())
	
	def Transform(self, matrix):
		self.P0 = matrix*self.P0
		self.V = matrix.getUpperLeft3x3() * self.V
		self.V.normalize()


	def intersectMeshFace(self, face):
		facePlane = face.getPlane()
		P0 = facePlane.P0
		N = facePlane.N
		P = self.P0
		V = self.V
		if abs(N.Dot(V)) < EPS:
			return None
		t = (P0.getVector().Dot(N) - N.Dot(P.getVector())) / (N.Dot(V))
		if t < 0:
			return None
		intersectP = P + t*V
		intersectN = N
		#Now check to see if the intersection is within the polygon
		#Do this by verifying that intersectP is on the same side
		#of each segment of the polygon
		verts = [v.pos for v in face.getVertices()]
		if len(verts) < 3:
			return None
		lastCross = (verts[1] - verts[0]) % (intersectP - verts[1])
		for i in range(1, len(verts)):
			v0 = verts[i]
			v1 = verts[(i+1)%len(verts)]
			cross = (v1 - v0) % (intersectP - v1)
			if cross.Dot(lastCross) < 0: #The point must be on the outside
				return None
			lastCross = cross
		return [t, intersectP, intersectN]
