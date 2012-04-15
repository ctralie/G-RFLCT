#TODO: These cameras need some serious work; I screwed up the coordinate systems
#and made them left-handed by accident 
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from Shapes3D import *
import math

class Key6DOFCamera(object):
	def __init__(self, eye, towards = Vector3D(0, 0, 1), up = Vector3D(0, 1, 0), yfov = 0.75):
		self.eye = eye
		self.towards = towards
		self.up = up
		self.yfov = yfov
		self.keys = {}
		self.LINEAR_RATE = 0.15
		self.ANGULAR_RATE = 0.01
		self.minTimestepLen = 1 #Minimum Timestep Length (msecs)
		self.zP = Point3D(0, 0, 0) #Zero Point (used for axis rotation)
	
	#Return a column-major matrix that puts a point in the camera's coordinate frame
	def getColMajorMatrix(self):
		#TODO: Make this (slightly) more efficient by doing transposition up-front
		(t, u) = (self.towards, self.up)
		r = t%u
		t.normalize()
		u.normalize()
		r.normalize()
		eye = self.eye
		trans = Matrix4([1, 0, 0, -eye.x, 0, 1, 0, -eye.y, 0, 0, 1, -eye.z, 0, 0, 0, 1])
		rot = Matrix4([r.x, r.y, r.z, 0, u.x, u.y, u.z, 0, t.x, t.y, t.z, 0, 0, 0, 0, 1])
		cameraMatrix = rot * trans
		return cameraMatrix.Transpose()
		
	def gotoCameraFrame(self):
		t = self.towards
		u = self.up
		r = t % u
		mat = [r.x, u.x, t.x, 0, r.y, u.y, t.y, 0, r.z, u.z, t.z, 0, 0, 0, 0, 1]
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glMultMatrixd(mat)
		glTranslated(self.eye.x, self.eye.y, self.eye.z)

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
			self.eye = self.eye - self.LINEAR_RATE*self.towards
		if self.keyPressed('s', 'S'):
			self.eye = self.eye + self.LINEAR_RATE*self.towards
		if self.keyPressed('a', 'A'):
			self.eye = self.eye - self.LINEAR_RATE*(self.towards%self.up)
		if self.keyPressed('d', 'D'):
			self.eye = self.eye + self.LINEAR_RATE*(self.towards%self.up)
		if self.keyPressed(' '):
			self.eye = self.eye + self.LINEAR_RATE*self.up
		if self.keyPressed('c', 'C'):
			self.eye = self.eye - self.LINEAR_RATE*self.up
		if self.keyPressed('e', 'E'):
			self.Roll(self.ANGULAR_RATE)
		if self.keyPressed('q', 'Q'):
			self.Roll(-self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_DOWN):
			self.Pitch(-self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_UP):
			self.Pitch(self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_LEFT):
			self.Yaw(-self.ANGULAR_RATE)
		if self.keyPressed(GLUT_KEY_RIGHT):
			self.Yaw(self.ANGULAR_RATE)
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
	def __init__(self, pixWidth, pixHeight, yfov = 0.75):
		self.pixWidth = pixWidth
		self.pixHeight = pixHeight
		self.yfov = yfov
		self.zP = Point3D(0, 0, 0) #Zero Point (used for axis rotation)
		self.center = Point3D(0, 0, 0)
		self.R = 1
		self.theta = 0 #Angle down from (0, 0, 1) pole
		self.phi = 0 #Angle around (0, 0, 1) axis
		self.updateVecsFromPolar()

	def centerOnMesh(self, mesh):
		bbox = mesh.getBBox()
		self.center = bbox.getCenter()
		self.R = bbox.ZLen()*3
		self.theta = 0 #Polar vector
		self.phi = 0
		self.updateVecsFromPolar()
	
	def updateVecsFromPolar(self):
		[sinT, cosT, sinP, cosP] = [math.sin(self.theta), math.cos(self.theta), math.sin(self.phi), math.cos(self.phi)]
		self.towards = Vector3D(-sinT*cosP, -sinT*sinP, -cosT)
		self.up = Vector3D(sinP, -cosP, 0)
		if (self.towards.squaredMag() < EPS):
			self.towards = Vector3D(0, 0, 1)
			self.up = Vector3D(0, 1, 0)
		self.eye = self.center - self.R*self.towards

	def gotoCameraFrame(self):
		t = self.towards
		u = self.up
		r = t % u
		mat = [r.x, u.x, t.x, 0, r.y, u.y, t.y, 0, r.z, u.z, t.z, 0, 0, 0, 0, 1]
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glMultMatrixd(mat)
		glTranslated(self.eye.x, self.eye.y, self.eye.z)
	
	def orbitUpDown(self, dP):
		dP = 1.5*dP/float(self.pixHeight)
		self.phi = self.phi+dP
		self.updateVecsFromPolar()
	
	def orbitLeftRight(self, dT):
		dT = 1.5*dT/float(self.pixWidth)
		self.theta = self.theta+dT
		self.updateVecsFromPolar()
	
	def zoom(self, rate):
		rate = rate / float(self.pixHeight)
		self.R = self.R + rate*pow(4, rate/self.R)
		if self.R < 0:
			self.R = 0
		self.updateVecsFromPolar()
	
	def translate(self, dx, dy):
		length = (self.center-self.eye).Length()*math.tan(self.yfov);
		dx = length*dx / float(self.pixWidth)
		dy = length*dy / float(self.pixHeight)
		r = self.towards % self.up
		self.center = self.center + dx*r + dy*self.up
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
		self.eye = Point3D(1, 0, 0)
		self.towards = Vector3D(0, 0, -1)
		self.up = Vector3D(0, 1, 0)

	def centerOnMesh(self, mesh):
		bbox = mesh.getBBox()
		self.center = bbox.getCenter()
		self.towards = Vector3D(0, 0, -1)
		self.up = Vector3D(0, 1, 0)
		self.eye = self.center - (bbox.XLen()*3)*self.towards

	def gotoCameraFrame(self):
		t = self.towards
		u = self.up
		r = t%u
		P = self.eye
		rotMat = Matrix4([r.x, u.x, -t.x, 0, r.y, u.y, -t.y, 0, r.z, u.z, -t.z, 0, 0, 0, 0, 1])
		rotMat = rotMat.Inverse()
		transMat = Matrix4([1, 0, 0, -P.x, 0, 1, 0, -P.y, 0, 0, 1, -P.z, 0, 0, 0, 1])
		mat = rotMat*transMat
		mat = mat.Transpose()
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glMultMatrixd(mat.m)
	
	def orbitUpDown(self, dTheta):
		dTheta = 1.5*dTheta / float(self.pixHeight)
		right = self.towards % self.up
		self.eye = rotateAroundAxis(self.center, right, dTheta, self.eye)
		self.up = rotateAroundAxis(Point3D(0, 0, 0), right, dTheta, self.up)
		self.towards = rotateAroundAxis(Point3D(0, 0, 0), right, dTheta, self.towards)
	
	def orbitLeftRight(self, dTheta):
		dTheta = 1.5*dTheta / float(self.pixWidth)
		self.eye = rotateAroundAxis(self.center, self.up, -dTheta, self.eye)
		self.towards = rotateAroundAxis(Point3D(0, 0, 0), self.up, -dTheta, self.towards)
	
	def zoom(self, rate):
		rate = 6.0*rate / float(self.pixHeight)
		R = self.center - self.eye
		self.eye = self.eye + (rate*pow(2,rate/R.squaredMag()))*self.towards
		self.towards = self.center - self.eye
		self.towards.normalize()
	
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
