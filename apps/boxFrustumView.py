from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from PolyMesh import *
from Cameras3D import *
from math import *
from sys import argv

class Viewer(object):
	def __init__(self):
		#GLUT State variables
		self.GLUTwindow_height = 1400-640
		self.GLUTwindow_width = 1400
		self.GLUTmouse = [0, 0]
		self.GLUTButton = [0, 0, 0, 0, 0]
		self.GLUTModifiers = 0
		self.keys = {}
		
		#Set up box parameters here
		#focal parameters
		kinectFOV = 525.0
		xScale = 320.0/kinectFOV
		yScale = 240.0/kinectFOV
		[self.xScale, self.yScale] = [xScale, yScale]
		boxDim = 4	#Assume everything is in feet
		self.boxDim = boxDim
		self.frustumPoints = [Point3D(-xScale*boxDim, yScale*boxDim, boxDim), Point3D(xScale*boxDim, yScale*boxDim, boxDim), Point3D(xScale*boxDim, -yScale*boxDim, boxDim), Point3D(-xScale*boxDim, -yScale*boxDim, boxDim)]
		self.frustumTrans = 0.3 #Transparency of viewing frustum
		
		self.camera = MouseSphericalCamera(self.GLUTwindow_width, self.GLUTwindow_height)
		self.camera.center = Point3D(0, 0, -boxDim/2)
		self.camera.eye = Point3D(0, 0, 2*boxDim)
		self.camera.zoom(0)
		
		self.kinectPos = Point3D(0, 1, 0.1)
	
		self.initGL()

	def GLUTResize(self, w, h):
		glViewport(0, 0, w, h)
		self.GLUTwindow_width = w
		self.GLUTwindow_height = h
		self.camera.pixWidth = w
		self.camera.pixHeight = h
		glutPostRedisplay()

	def drawScene(self, drawFrustum = True):		
		glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
		glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
		
		glEnable(GL_LIGHTING)
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [100.0/256.0, 50.0/256.0, 0, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)

		#Do rendering here
		#STEP 1: Draw chamber box
		d = self.boxDim/2.0
		#No front face
		
		#Back Face
		glBegin(GL_POLYGON)
		glNormal3f(0, 0, 1)
		glVertex3f(-d, -d, d*2)
		glVertex3f(-d, d, d*2)
		glVertex3f(d, d, d*2)
		glVertex3f(d, -d, d*2)
		glEnd()
		
		#Left Face
		glBegin(GL_POLYGON)
		glNormal3f(-1, 0, 0)
		glVertex3f(-d, -d, 0)
		glVertex3f(-d, -d, 2*d)
		glVertex3f(-d, d, 2*d)
		glVertex3f(-d, d, 0)
		glEnd()
		
		#Right Face
		glBegin(GL_POLYGON)
		glNormal3f(1, 0, 0)
		glVertex3f(d, -d, 0)
		glVertex3f(d, -d, 2*d)
		glVertex3f(d, d, 2*d)
		glVertex3f(d, d, 0)
		glEnd()
		
		#Top Face
		glBegin(GL_POLYGON)
		glNormal3f(0, 1, 0)
		glVertex3f(-d, d, 0)
		glVertex3f(-d, d, 2*d)
		glVertex3f(d, d, 2*d)
		glVertex3f(d, d, 0)
		glEnd()
		
		#Bottom Face
		glBegin(GL_POLYGON)
		glNormal3f(0, -1, 0)
		glVertex3f(-d, -d, 0)
		glVertex3f(-d, -d, 2*d)
		glVertex3f(d, -d, 2*d)
		glVertex3f(d, -d, 0)
		glEnd()
		
		
		#STEP 2: Draw Spheres
		quadric = gluNewQuadric()
		glPushMatrix()
		glTranslatef(1, 2-17.0/12.0, 3)
		gluSphere(quadric, 0.25, 32, 32)
		glPopMatrix()
		
		if drawFrustum:
			#STEP 3: Draw Kinect frustum
			glDisable(GL_LIGHTING)
			glEnable(GL_BLEND)
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
			P0 = self.kinectPos
			Ps = [P + self.kinectPos for P in self.frustumPoints]
			glColor4f(0, 1, 0, self.frustumTrans)
			glBegin(GL_TRIANGLES)
			for i in range(0, 4):
				glVertex3f(P0.x, P0.y, P0.z)
				[P1, P2] = [Ps[i], Ps[(i+1)%4]]
				glVertex3f(P1.x, P1.y, P1.z)
				glVertex3f(P2.x, P2.y, P2.z)
			glEnd()
			glColor3f(0, 0, 1)
			glLineWidth(5)
			glBegin(GL_LINES)
			for P in Ps:
				glVertex3f(P0.x, P0.y, P0.z)
				glVertex3f(P.x, P.y, P.z)
			glEnd()

	def GLUTRedraw(self):
		glClearColor(0.0, 0.0, 0.0, 0.0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		#First draw the scene from the viewer's point of view on the left
		[w, h] = [self.GLUTwindow_width, self.GLUTwindow_height]
		glViewport(0, 0, w-640, h)
		glScissor(0, 0, w-640, h)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(180.0*self.camera.yfov/M_PI, float(w-640)/self.GLUTwindow_height, 0.01, 100.0)
		#Set up modelview matrix
		self.camera.gotoCameraFrame()
		self.drawScene()
		
		#Now drawn the scene from the camera's point of view on the right
		glViewport(w-640, 0, 640, 480)
		glScissor(w-640, 0, 640, 480)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		[near, far] = [0.1, 10.0]
		[cameraw, camerah] = [near*self.xScale, near*self.yScale]
		glFrustum(-cameraw, cameraw, -camerah, camerah, near, far)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		pos = self.kinectPos
		glMultMatrixd([-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
		gluLookAt(pos.x, pos.y, pos.z, pos.x, pos.y, self.boxDim, 0, 1, 0)
		self.drawScene(drawFrustum = False)
		
		glutSwapBuffers()

	def handleMouseStuff(self, x, y):
		y = self.GLUTwindow_height - y
		self.GLUTmouse[0] = x
		self.GLUTmouse[1] = y
		self.GLUTmodifiers = glutGetModifiers()
	
	def GLUTKeyboard(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.keys[key] = True
		glutPostRedisplay()
	
	def GLUTKeyboardUp(self, key, x, y):
		self.handleMouseStuff(x, y)
		if key in ['x', 'X']:
			self.frustumTrans = min(self.frustumTrans + 0.1, 1)
		elif key in ['z', 'Z']:
			self.frustumTrans = max(self.frustumTrans - 0.1, 0)
		elif key in ['a', 'A']:
			self.kinectPos = self.kinectPos + Point3D(-0.1, 0, 0)
		elif key in ['d', 'D']:
			self.kinectPos = self.kinectPos + Point3D(0.1, 0, 0)
		elif key in ['w', 'W']:
			self.kinectPos = self.kinectPos + Point3D(0, 0.1, 0)
		elif key in ['s', 'S']:
			self.kinectPos = self.kinectPos + Point3D(0, -0.1, 0)
		glutPostRedisplay()
	
	def GLUTSpecial(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.keys[key] = True
		glutPostRedisplay()
	
	def GLUTSpecialUp(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.keys[key] = False
		glutPostRedisplay()
		
	def GLUTMouse(self, button, state, x, y):
		buttonMap = {GLUT_LEFT_BUTTON:0, GLUT_MIDDLE_BUTTON:1, GLUT_RIGHT_BUTTON:2, 3:3, 4:4}
		if state == GLUT_DOWN:
			self.GLUTButton[buttonMap[button]] = 1
		else:
			self.GLUTButton[buttonMap[button]] = 0
		self.handleMouseStuff(x, y)
		glutPostRedisplay()

	def GLUTMotion(self, x, y):
		lastX = self.GLUTmouse[0]
		lastY = self.GLUTmouse[1]
		self.handleMouseStuff(x, y)
		dX = self.GLUTmouse[0] - lastX
		dY = self.GLUTmouse[1] - lastY
		if self.GLUTButton[2] == 1:
			self.camera.zoom(-dY)
		elif self.GLUTButton[1] == 1:
			self.camera.translate(dX, dY)
		else:
			self.camera.orbitLeftRight(-dX)
			self.camera.orbitUpDown(-dY)
		glutPostRedisplay()
	
	def initGL(self):
		glutInit('')
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
		glutInitWindowSize(self.GLUTwindow_width, self.GLUTwindow_height)
		glutInitWindowPosition(50, 50)
		glutCreateWindow('Frustum Viewer')
		glutReshapeFunc(self.GLUTResize)
		glutDisplayFunc(self.GLUTRedraw)
		glutKeyboardFunc(self.GLUTKeyboard)
		glutKeyboardUpFunc(self.GLUTKeyboardUp)
		glutSpecialFunc(self.GLUTSpecial)
		glutSpecialUpFunc(self.GLUTSpecialUp)
		glutMouseFunc(self.GLUTMouse)
		glutMotionFunc(self.GLUTMotion)
		
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glEnable(GL_LIGHT0)
		glLightfv(GL_LIGHT1, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
		glEnable(GL_LIGHT1)
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHTING)		
		
		glEnable(GL_DEPTH_TEST)
		
		glutMainLoop()

if __name__ == '__main__':
	viewer = Viewer()
