from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from Shapes3D import *
import math

def getCameraMatrix(t, u, r, P):
	rotMat = Matrix4([r.x, u.x, -t.x, 0, r.y, u.y, -t.y, 0, r.z, u.z, -t.z, 0, 0, 0, 0, 1])
	rotMat = rotMat.Inverse()
	transMat = Matrix4([1, 0, 0, -P.x, 0, 1, 0, -P.y, 0, 0, 1, -P.z, 0, 0, 0, 1])
	#Translate first then rotate
	mat = rotMat*transMat
	return mat

#This function pushes a matrix onto the stack that puts everything
#in the frame of a camera which is centered at position "P",
#is pointing towards "t", and has vector "r" to the right
#t - towards vector
#u - up vector
#r - right vector
#P - Camera center
def gotoCameraFrame(t, u, r, P):
	rotMat = Matrix4([r.x, u.x, -t.x, 0, r.y, u.y, -t.y, 0, r.z, u.z, -t.z, 0, 0, 0, 0, 1])
	rotMat = rotMat.Inverse()
	transMat = Matrix4([1, 0, 0, -P.x, 0, 1, 0, -P.y, 0, 0, 1, -P.z, 0, 0, 0, 1])
	#Translate first then rotate
	mat = rotMat*transMat
	#OpenGL is column major and mine are row major so take transpose
	mat = mat.Transpose()
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	glMultMatrixd(mat.m)


class Key6DOFCamera(object):
	def __init__(self, eye, towards = Vector3D(0, 0, 1), up = Vector3D(0, 1, 0), yfov = 0.75):
		self.eye = eye
		self.towards = towards
		self.up = up
		self.yfov = yfov
		self.keys = {}
		self.LINEAR_RATE = 0.02
		self.ANGULAR_RATE = 0.01
		self.minTimestepLen = 1 #Minimum Timestep Length (msecs)
		self.zP = Point3D(0, 0, 0) #Zero Point (used for axis rotation)
		
	def gotoCameraFrame(self):
		gotoCameraFrame(self.towards, self.up, self.towards%self.up, self.eye)

	def Roll(self, theta):
		self.up = rotateAroundAxis(self.zP, self.towards, theta, self.up)
	
	def Pitch(self, theta):
		right = self.towards % self.up
		self.up = rotateAroundAxis(self.zP, right, theta, self.up)
		self.towards = rotateAroundAxis(self.zP, right, theta, self.towards)
	
	def Yaw(self, theta):
		self.towards = rotateAroundAxis(self.zP, self.up, theta, self.towards)

	def keyPressed(self, key1, key2 = 'null'):
		if (key1 in self.keys and self.keys[key1]) or (key2 in self.keys and self.keys[key2]):
			return True
		return False

	def handleUserInput(self, value):
		if self.keyPressed('w', 'W'):
			self.eye = self.eye + self.LINEAR_RATE*self.towards
		if self.keyPressed('s', 'S'):
			self.eye = self.eye - self.LINEAR_RATE*self.towards
		if self.keyPressed('a', 'A'):
			self.eye = self.eye - self.LINEAR_RATE*(self.towards%self.up)
		if self.keyPressed('d', 'D'):
			self.eye = self.eye + self.LINEAR_RATE*(self.towards%self.up)
		if self.keyPressed(' '):
			self.eye = self.eye + self.LINEAR_RATE*self.up
		if self.keyPressed('c', 'C'):
			self.eye = self.eye - self.LINEAR_RATE*self.up
		if self.keyPressed('e', 'E'):
			self.Roll(-self.ANGULAR_RATE)
		if self.keyPressed('q', 'Q'):
			self.Roll(self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_DOWN):
			self.Pitch(self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_UP):
			self.Pitch(-self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_LEFT):
			self.Yaw(self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_RIGHT):
			self.Yaw(-self.ANGULAR_RATE)
		glutTimerFunc(self.minTimestepLen, self.handleUserInput, 1)
		glutPostRedisplay()

	def __str__(self):
		right = self.towards % self.up
		str = "Eye: %s\n"%(self.eye)
		str = str+"Towards: (%g) %s\n"%(self.towards.Length(), self.towards)
		str = str + "Up: (%g) %s\n"%(self.up.Length(), self.up)
		str = str + "Right: (%g) %s"%(right.Length(), right)
		return str

class MousePolarCamera(object):
	#Coordinate system is defined as in OpenGL as a right
	#handed system with +z out of the screen, +x to the right,
	#and +y up
	#phi is CCW down from +y, theta is CCW away from +z
	def __init__(self, pixWidth, pixHeight, yfov = 0.75):
		self.pixWidth = pixWidth
		self.pixHeight = pixHeight
		self.yfov = yfov
		self.center = Point3D(0, 0, 0)
		self.R = 1
		self.theta = 0 
		self.phi = 0 
		self.updateVecsFromPolar()

	def centerOnMesh(self, mesh):
		bbox = mesh.getBBox()
		self.center = bbox.getCenter()
		self.R = bbox.getDiagLength()*3
		self.theta = -math.pi/2
		self.phi = math.pi/2
		self.updateVecsFromPolar()

	def updateVecsFromPolar(self):
		[sinT, cosT, sinP, cosP] = [math.sin(self.theta), math.cos(self.theta), math.sin(self.phi), math.cos(self.phi)]
		#Make the camera look inwards
		#i.e. towards is -dP(R, phi, theta)/dR, where P(R, phi, theta) is polar position
		self.towards = Vector3D(-sinP*cosT, -cosP, sinP*sinT)
		self.up = Vector3D(-cosP*cosT, sinP, cosP*sinT)
		self.eye = self.center - self.R*self.towards

	def gotoCameraFrame(self):
		gotoCameraFrame(self.towards, self.up, self.towards%self.up, self.eye)
	
	def orbitUpDown(self, dP):
		dP = 1.5*dP/float(self.pixHeight)
		self.phi = self.phi+dP
		self.updateVecsFromPolar()
	
	def orbitLeftRight(self, dT):
		dT = 1.5*dT/float(self.pixWidth)
		self.theta = self.theta-dT
		self.updateVecsFromPolar()
	
	def zoom(self, rate):
		rate = rate / float(self.pixHeight)
		self.R = self.R*pow(4, rate)
		self.updateVecsFromPolar()
	
	def translate(self, dx, dy):
		length = (self.center-self.eye).Length()*math.tan(self.yfov);
		dx = length*dx / float(self.pixWidth)
		dy = length*dy / float(self.pixHeight)
		r = self.towards % self.up
		self.center = self.center - dx*r - dy*self.up
		self.updateVecsFromPolar()

	def __str__(self):
		right = self.towards % self.up
		str = "R = %g, theta = %g, phi = %g\n"%(self.R, self.theta, self.phi)
		str = str + "Eye: %s\n"%(self.eye)
		str = str + "Towards: (%g) %s\n"%(self.towards.Length(), self.towards)
		str = str + "Up: (%g) %s\n"%(self.up.Length(), self.up)
		str = str + "Right: (%g) %s"%(right.Length(), right)
		str = str + "\nCenter: %s\n"%(self.center)
		str = str + "Eye.Dot(Towards) = %g\n"%(self.towards.Dot(self.eye)/self.eye.Length())
		return str


#This is different from a polar camera because it continually rotates the pole
#to align with where the user is looking (easier to use locally but harder
#to align to axes once you've stared moving)
class MouseSphericalCamera(object):
	def __init__(self, pixWidth, pixHeight, yfov = 0.75):
		self.pixWidth = pixWidth
		self.pixHeight = pixHeight
		self.yfov = yfov
		self.zP = Point3D(0, 0, 0) #Zero Point (used for axis rotation)
		self.center = Point3D(0, 0, 0)
		self.eye = Point3D(0, 0, 1)
		self.towards = Vector3D(0, 0, -1)
		self.up = Vector3D(0, 1, 0)

	def centerOnMesh(self, mesh):
		bbox = mesh.getBBox()
		self.center = bbox.getCenter()
		self.towards = Vector3D(0, 0, -1)
		self.up = Vector3D(0, 1, 0)
		self.eye = self.center - (bbox.getDiagLength()*3)*self.towards

	def gotoCameraFrame(self):
		gotoCameraFrame(self.towards, self.up, self.towards%self.up, self.eye)
	
	def orbitUpDown(self, dTheta):
		dTheta = 1.5*dTheta / float(self.pixHeight)
		right = self.towards % self.up
		self.eye = rotateAroundAxis(self.center, right, dTheta, self.eye)
		self.up = rotateAroundAxis(self.zP, right, dTheta, self.up)
		self.towards = rotateAroundAxis(self.zP, right, dTheta, self.towards)
	
	def orbitLeftRight(self, dTheta):
		dTheta = 1.5*dTheta / float(self.pixWidth)
		self.eye = rotateAroundAxis(self.center, self.up, -dTheta, self.eye)
		self.towards = rotateAroundAxis(self.zP, self.up, -dTheta, self.towards)

	def zoom(self, rate):
		rate = rate / float(self.pixHeight)
		R = self.eye - self.center
		R = pow(4, rate)*R
		self.eye = self.center + R
	
	def translate(self, dx, dy):
		length = (self.center-self.eye).Length()*math.tan(self.yfov);
		dx = length*dx / float(self.pixWidth)
		dy = length*dy / float(self.pixHeight)
		right = self.towards % self.up
		self.eye = self.eye - dx*right - dy*self.up
		self.center = self.center - dx*right - dy*self.up

	def __str__(self):
		right = self.towards % self.up
		str = "Eye: %s\n"%(self.eye)
		str = str + "Towards: (%g) %s\n"%(self.towards.Length(), self.towards)
		str = str + "Up: (%g) %s\n"%(self.up.Length(), self.up)
		str = str + "Right: (%g) %s"%(right.Length(), right)
		return str
