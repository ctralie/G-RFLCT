from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from PolyMesh import *
from Cameras3D import *
from sys import argv
import random

class Viewer(object):
	def __init__(self):
		#GLUT State variables
		self.GLUTwindow_height = 800
		self.GLUTwindow_width = 800
		self.GLUTmouse = [0, 0]
		self.GLUTButton = [0, 0, 0, 0, 0]
		self.GLUTModifiers = 0
		self.keys = {}
		self.drawEdges = 0
		self.drawVerts = 0
		self.drawNormals = 0
		
		#Camera state variables
		self.camera = MouseSphericalCamera(self.GLUTwindow_width, self.GLUTwindow_height)
		
		self.initGL()

	def GLUTResize(self, w, h):
		glViewport(0, 0, w, h)
		self.GLUTwindow_width = w
		self.GLUTwindow_height = h
		self.camera.pixWidth = w
		self.camera.pixHeight = h
		glutPostRedisplay()

	def GLUTRedraw(self):
		#Set up projection matrix
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(180.0*self.camera.yfov/M_PI, float(self.GLUTwindow_width)/self.GLUTwindow_height, 0.01, 100.0)
		
		#Set up modelview matrix
		self.camera.gotoCameraFrame()	
		glClearColor(0.0, 0.0, 0.0, 0.0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
		glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
		
		glEnable(GL_LIGHTING)
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [1, 0, 0, 1.0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [1, 0.2, 0.2, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)

		#Now draw tetrahedron
		#4 corner spheres
		quadric = gluNewQuadric()
		glPushMatrix()
		glTranslatef(-1, 1, 1)
		gluSphere(quadric, 0.5, 32, 32)
		glPopMatrix()
		glPushMatrix()
		glTranslatef(1, 1, -1)
		gluSphere(quadric, 0.5, 32, 32)
		glPopMatrix()
		glPushMatrix()
		glTranslatef(1, -1, 1)
		gluSphere(quadric, 0.5, 32, 32)
		glPopMatrix()
		glPushMatrix()
		glTranslatef(-1, -1, -1)
		gluSphere(quadric, 0.5, 32, 32)
		glPopMatrix()
		
		#Center sphere
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0, 0.5, 1, 1.0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [1, 0.2, 0.2, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
		gluSphere(quadric, 0.5, 32, 32)
		
		glColor3f(0.5, 0.5, 0.5)
		#Attachments
		glLineWidth(30)
		glBegin(GL_LINES)
		glVertex3f(-1, 1, 1)
		glVertex3f(0, 0, 0)
		glVertex3f(1, 1, -1)
		glVertex3f(0, 0, 0)
		glVertex3f(1, -1, 1)
		glVertex3f(0, 0, 0)
		glVertex3f(-1, -1, -1)
		glVertex3f(0, 0, 0)						
		glEnd()
		
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
		self.keys[key] = False
		if key in ['e', 'E']:
			self.drawEdges = 1 - self.drawEdges
		elif key in ['v', 'V']:
			self.drawVerts = 1 - self.drawVerts
		elif key in ['r', 'R']:
			self.mesh.splitFaces()
			print self.mesh
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
		glutCreateWindow('Viewer')
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
