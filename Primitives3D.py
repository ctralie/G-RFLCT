EPS = 1e-7
M_PI = 3.1415925
import math

class Vector3D(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	
	def Dot(self, other):
		return self.x*other.x + self.y*other.y + self.z*other.z
		
	def squaredMag(self):
		return self.x**2.0 + self.y**2.0 + self.z**2.0
	
	def Length(self):
		return self.squaredMag()**0.5
	
	def normalize(self):
		mag = self.squaredMag()**0.5
		mag = float(mag)
		if (mag > EPS):
			self.x = self.x/mag
			self.y = self.y/mag
			self.z = self.z/mag
	
	def __add__(self, other):
		if isinstance(other, Point3D):#A Point plus a vector is still a point
			return other + self
		return Vector3D(self.x+other.x, self.y+other.y, self.z+other.z)

	def __sub__(self, other):
		return Vector3D(self.x-other.x, self.y-other.y, self.z-other.z)
	
	def __rmul__(self, a):
		if isinstance(a, float) or isinstance(a, int):
			return Vector3D(a*self.x, a*self.y, a*self.z)
	
	def __mul__(self, a):
		if isinstance(a, float) or isinstance(a, int):
			return Vector3D(a*self.x, a*self.y, a*self.z)
	
	def __neg__(self):
		return Vector3D(-self.x, -self.y, -self.z)
	
	def __eq__(self, other):
		return (self.x == other.x and self.y == other.y and self.z == other.z)
		
	def __ne__(self, other):
		return not (self == other)
	
	#Dot product
	def Dot(self, other):
		return self.x*other.x + self.y*other.y + self.z*other.z
	
	#Use this to implement cross product
	def __mod__(self, other):
		newX = self.y*other.z - other.y*self.z
		newY = -self.x*other.z + other.x*self.z
		newZ = self.x*other.y - other.x*self.y
		return Vector3D(newX, newY, newZ)
		
	#Project V onto this vector
	def proj(self, V):
		thisSqrMag = self.squaredMag()
		scale = 0
		if (thisSqrMag > EPS):
			scale = float(self.Dot(V)) / float(thisSqrMag)
		return Vector3D(self.x*scale, self.y*scale, self.z*scale)
		
	#Do the perpendicular projection of V onto this vector
	def projPerp(self, V):
		ret = Vector3D(V.x, V.y, V.z)
		parProj = self.proj(V)
		return ret - parProj
	
	def Copy(self):
		return Vector3D(self.x, self.y, self.z)
	
	def getPoint(self):
		return Point3D(self.x, self.y, self.z)
	
	def __str__(self):
		return "Vector3D(%g, %g, %g)"%(self.x, self.y, self.z)

class Point3D(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def squaredMag(self):
		return self.x**2.0 + self.y**2.0 + self.z**2.0
	
	def Length(self):
		return self.squaredMag()**0.5
	
	def __add__(self, other):
		return Point3D(self.x+other.x, self.y+other.y, self.z+other.z)

	def __radd__(self, other):
		return self + other
	
	def __sub__(self, other):
		if isinstance(other, Vector3D): #A point minus a vector is still a point
			return Point3D(self.x-other.x, self.y - other.y, self.z - other.z)
		return Vector3D(self.x-other.x, self.y-other.y, self.z-other.z)	
	
	def __rmul__(self, a):
		if isinstance(a, float) or isinstance(a, int):
			return Point3D(a*self.x, a*self.y, a*self.z)
	
	def __mul__(self, a):
		if isinstance(a, float) or isinstance(a, int):
			return Point3D(a*self.x, a*self.y, a*self.z)
	
	def __neg__(self):
		return Point3D(-self.x, -self.y, -self.z)

	def __eq__(self, other):
		return (self.x == other.x and self.y == other.y and self.z == other.z)
		
	def __ne__(self, other):
		return not (self == other)

	def Copy(self):
		return Point3D(self.x, self.y, self.z)

	def getVector(self):
		return Vector3D(self.x, self.y, self.z)

	def __str__(self):
		return "Point3D(%g, %g, %g)"%(self.x, self.y, self.z)

class Matrix4(object):
	def __init__(self, args = []):
		if len(args) > 0: #Matrix specified in row-major order in argument
			self.m = args[:]
		else: #Make an identity matrix by default
			self.m = [0]*16 #Zeros
			self.m[0:16:5] = [1, 1, 1, 1] #Ones along diagonal
	
	def Transpose(self):
		return Matrix4([self.m[4*(i%4)+i/4] for i in range(0, 16)])
	
	def getColVector3(self, col):
		if col < 0 or col > 3:
			print "Column %i out of range\n"%(col)
			return
		return Vector3D(self.m[col], self.m[col+4], self.m[col+8])
	
	def __eq__(self, other):
		return sum([abs(self.m[i]-other.m[i]) for i in range(0, 16)]) == 0
	
	def __ne__(self, other):
		return not (self == other)
	
	def __add__(self, other):
		return Matrix4([self.m[i]+other.m[i] for i in range(0, 16)])
	
	def __sub__(self, other):
		return Matrix4([self.m[i]-other.m[i] for i in range(0, 16)])
	
	def __mul__(self, other):
		if isinstance(other, Matrix4):
			#Matrix-Matrix multiplication
			retm = [0]*16
			for i in range(0, 16):
				(row, col) = (i/4, i%4)
				retm[i] = sum([self.m[row*4+k]*other.m[k*4+col] for k in range(0, 4)])
			return Matrix4(retm)
		if isinstance(other, Point3D) or isinstance(other, Vector3D):
			#Point/vector transformation by matrix
			retv = [0]*4
			otherm = [other.x, other.y, other.z, 1] #4D homogenous point/vector
			for row in range(0, 4):
				retv[row] = sum([self.m[row*4+k]*otherm[k] for k in range(0, 4)])
			#Divide by homogenous coordinate
			retv = [float(retv[k]) / float(retv[3]) for k in range(0, 3)] 
			if isinstance(other, Point3D):
				return Point3D(retv[0], retv[1], retv[2])
			else:
				return Vector3D(retv[0], retv[1], retv[2])
	
	#http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	def Inverse(self):
		m = self.m
		inv = [0]*16
		inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10]
		inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
		inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]	+ m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
		inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]	- m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
		inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
		inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
		inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]	- m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
		inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]	+ m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
		inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]	+ m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
		inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
		inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
		inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
		inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
		inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
		inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
		inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10] + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];
		det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12]
		if det == 0:
			print "ERROR: Determinant of matrix is zero; cannot invert"
			return
		det = 1.0/det
		return Matrix4([inv[i]*det for i in range(0, 16)])

	def getUpperLeft3x3(self):
		m = self.m[:]
		m[3] = 0
		m[7] = 0
		m[11] = 0
		return Matrix4(m)

	def __str__(self):
		fmt = "[%g, %g, %g, %g]\n"*4
		return fmt%(tuple(self.m))

#Rotate a vector (or point) V by angle "theta" around a line through P0 
#whose direction is specified by V
def rotateAroundAxis(P0, axis, theta, V):
	#print "P0 = %s, axis = %s, V = %s"%(P0, axis, V)
	diffV = V - P0
	parV = axis.proj(diffV) #Part of v along axis unaffected by rotation
	perpV = axis.projPerp(diffV)
	if perpV.Dot(perpV) < EPS: #Hardly any perpendicular component
		return V
	u = perpV.Copy()
	u.normalize()
	w = axis
	w.normalize()
	v = w%u
	#Put perpV into a frame where the rotation is about (0, 0, 1)
	fromFrame = Matrix4([u.x, v.x, w.x, 0, u.y, v.y, w.y, 0, u.z, v.z, w.z, 0, 0, 0, 0, 1])
	toFrame = fromFrame.Inverse()
	perpV = toFrame*perpV
	#Rotate perpV by theta in that frame
	(cosTheta, sinTheta) = (math.cos(theta), math.sin(theta))
	rotMat = Matrix4([cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
	perpV = rotMat*perpV
	#Bring perpV back into world coordinates and compute result
	perpV = fromFrame*perpV
	res = P0 + perpV + parV
	if isinstance(V, Vector3D):
		return Vector3D(res.x, res.y, res.z)
	if isinstance(V, Point3D):
		return Point3D(res.x, res.y, res.z)

#Return True if the vertices in the list "verts" all lie
#in the same plane and False otherwise
def arePlanar(verts):
	if len(verts) <= 3:
		return True
	v0 = verts[1] - verts[0]
	v1 = verts[2] - verts[0]
	n = v0 % v1
	n.normalize()
	for i in range(3, len(verts)):
		v = verts[i] - verts[0]
		v.normalize()
		if abs(v.Dot(n)) > EPS:
			return False
	return True

#If the vertices in "verts" form a convex 2D polygon 
#(in the order specified) return True.  Return False otherwise
def are2DConvex(verts):
	if len(verts) <= 3:
		return True
	if not arePlanar(verts):
		return False
	v0 = verts[0]
	v1 = verts[1]
	v2 = verts[2]
	lastCross = (v1-v0) % (v2-v1)
	for i in range(3, len(verts)+1):
		v0 = v1
		v1 = v2
		v2 = verts[i%len(verts)]
		cross = (v1-v0) % (v2-v1)
		if cross.Dot(lastCross) < 0:
			return False
		lastCross = cross
	return True

class Plane3D(object):
	#P0 is some point on the plane, N is the normal
	#Also store A, B, C, and D, the coefficients of the implicit plane equation
	def __init__(self, P0 = Point3D(0, 0, 0), N = Vector3D(0, 0, 1)):
		self.P0 = P0
		self.N = N
		self.N.normalize()
		self.resetEquation()

	def resetEquation(self):
		[self.A, self.B, self.C] = [self.N.x, self.N.y, self.N.z]
		self.D = -self.P0.getVector().Dot(self.N)

	def initFromEquation(self, A, B, C, D):
		self.N = Vector3D(A, B, C)
		self.P0 = Point3D(A, B, D)
		self.P0 = (-D/self.N.MagSquared())*self.P0
		self.N.normalize()
		self.resetEquation()

	def distFromPlane(self, P):
		return self.A*P.x + self.B*P.y + self.C*P.z + self.D

	def __str__(self):
		return "Plane3D: %g*x + %g*y + %g*z + %g = 0"%(self.A, self.B, self.C, self.D)

if __name__ == '__main__':
	P = Plane3D(Point3D(1, 1, 1), Vector3D(1, 2, 3))
	print P
	angle = 30
	angle = angle*3.141/180.0
	(cosA, sinA) = [math.cos(angle), math.sin(angle)]
	A = Matrix4([cosA, -sinA, 0, 0,
				sinA, cosA, 0, 0,
				0, 0, 1, 0, 0, 0, 0, 1])
	verts = [Point3D(0, 0, 0), Point3D(0, 1, 0), Point3D(1, 1, 0), Point3D(2, 0.5, 0), Point3D(3.01, 0, 0), Point3D(1, -1, 0)]
	verts = [rotateAroundAxis(Point3D(1, 1, 0), Vector3D(1, 0, 0), angle, v) for v in verts]
	for v in verts:
		print v
	print are2DConvex(verts)
