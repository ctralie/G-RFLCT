from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from EMScene import *
from RayTraceImage import *
from sys import argv
import random

SHOWIMAGES = False
DRAWPATHS = False

class Viewer(object):
	def __init__(self, filename):
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
		self.scene = EMScene()
		self.scene.Read(filename)
		#scene = self.scene
		#node = scene.rootEMNode
		#boxMesh = getBoxMesh(5.1, 5.1, 2.5, Point3D(0, 0, 0))
		#EMMat = EMMaterial(0.9, 0)
		#sceneNode = EMNode(scene.rootEMNode, boxMesh, Matrix4(), EMMat, OpticalMaterial())
		#scene.rootEMNode.children.append(sceneNode)
		#scene.Source = Point3D(0, 0, -2.0)
		#scene.Receiver = Point3D(0, -1.0, 2.0)
		#scene.getMeshList()
		#scene.buildVirtualSourceTree(3)
		#scene.getPathsToReceiver()
		#response = scene.getSteadyStateSinusoid(915, 40, 10)
		#print "times = %s; signal = %s"%(response[0], response[1])
		
		self.camera = MouseSphericalCamera(self.GLUTwindow_width, self.GLUTwindow_height)
		random.seed()
		self.rayPoints = []
		self.rayNormals = []
		self.eyePoints = []
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

		self.scene.renderGL()
		
		glDisable(GL_LIGHTING)
		glColor3f(1, 0, 0)
		glBegin(GL_LINES)
		for P in self.rayPoints:
			glVertex3f(P.x, P.y, P.z)
		glEnd()
		
		glDisable(GL_LIGHTING)
		glColor3f(0, 0, 1)
		glBegin(GL_LINES)
		for P in self.rayNormals:
			glVertex3f(P.x, P.y, P.z)
		glEnd()		
		
		#self.eyePoints.append(self.camera.eye)
		#glDisable(GL_LIGHTING)
		#glPointSize(5)
		#glBegin(GL_POINTS)
		#for P in self.eyePoints:
		#	glVertex3f(P.x, P.y, P.z)
		#glEnd()
		
		
		if isinstance(self.scene.Source, Point3D):
			glDisable(GL_LIGHTING)
			glColor3f(1, 0, 0)
			P = self.scene.Source
			quadric = gluNewQuadric()
			glPushMatrix()
			glTranslatef(P.x, P.y, P.z)
			gluSphere(quadric, 0.1, 32, 32)
			glPopMatrix()
		
		if isinstance(self.scene.Receiver, Point3D):
			glDisable(GL_LIGHTING)
			glColor3f(0, 0, 1)
			P = self.scene.Receiver
			quadric = gluNewQuadric()
			glPushMatrix()
			glTranslatef(P.x, P.y, P.z)
			gluSphere(quadric, 0.1, 32, 32)
			glPopMatrix()
		
		if DRAWPATHS == True:
			glDisable(GL_LIGHTING)
			glColor3f(1, 0, 0)
			glBegin(GL_LINES)
			for path in self.scene.paths:
				for i in range(0, len(path)-1):
					P0 = path[i]
					P1 = path[(i+1)]
					glVertex3f(P0.x, P0.y, P0.z)
					glVertex3f(P1.x, P1.y, P1.z)
			glEnd()
		
		if SHOWIMAGES == True:
			if False:
				if isinstance(self.scene.Source, Point3D):
					glDisable(GL_LIGHTING)
					glColor3f(0, 1, 0)
					glBegin(GL_LINES)
					P0 = self.scene.Source
					for source in self.scene.vSources:
						P1 = source.pos
						glVertex3f(P0.x, P0.y, P0.z)
						glVertex3f(P1.x, P1.y, P1.z)				
					glEnd()
		
				if isinstance(self.scene.Receiver, Point3D):
					glDisable(GL_LIGHTING)
					glColor3f(0, 1, 1)
					glBegin(GL_LINES)
					P0 = self.scene.Receiver
					for source in self.scene.vSources:
						P1 = source.pos
						glVertex3f(P0.x, P0.y, P0.z)
						glVertex3f(P1.x, P1.y, P1.z)				
					glEnd()			
		
			glDisable(GL_LIGHTING)
			glColor3f(1, 0, 1)
			quadric = gluNewQuadric()
			for source in self.scene.vSources:
				P = source.pos
				glPushMatrix()
				glTranslatef(P.x, P.y, P.z)
				gluSphere(quadric, 0.1, 32, 32)
				glPopMatrix()		
		
			glDisable(GL_LIGHTING)
			glBegin(GL_LINES)
			for ray in self.scene.rays:
				P0 = ray.P0
				P1 = ray.P0 + 5*ray.V
				glColor3f(0, 0.1, 0)
				glVertex3f(P0.x, P0.y, P0.z)
				glColor3f(1, 1, 1)
				glVertex3f(P1.x, P1.y, P1.z)
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
			#Launch some rays for debugging
			res = 10
			self.rayPoints = []
			self.rayNormals = []
			towards = self.camera.towards
			up = self.camera.up
			right = towards % up
			for x in range(-res, res+1):
				for y in range(-res, res+1):
					direct = towards + float(x)/float(res)*right + float(y)/float(res)*up
					direct.normalize()
					ray = Ray3D(self.camera.eye, direct)
					intersection = self.scene.getRayIntersection(ray)
					if intersection != None:
						self.rayPoints.append(self.camera.eye)
						self.rayPoints.append(intersection[1])
						self.rayNormals.append(intersection[1])
						self.rayNormals.append(intersection[1]+0.1*intersection[2])
		elif key in ['t', 'T']:
			self.rayNormals = []
			(self.rayPoints, self.rayNormals) = RayTraceImage(self.scene, self.camera, 50, 50, "out.png")
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
			self.camera.orbitLeftRight(dX)
			self.camera.orbitUpDown(dY)
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
	if len(argv) < 2:
		print "Usage: sceneView <scene filepath>"
	else:
		viewer = Viewer(argv[1])
