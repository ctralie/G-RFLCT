from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from EMScene import *
from Beam3D import *
import Image
import matplotlib.pyplot as plt
from sys import argv

DRAW_BACKPROJECTED = True
KEYDELTA = 0.01

class Viewer(object):
	def updateCameraVars(self):
		self.yfov = self.camera.yfov
		self.xfov = self.yfov*float(self.pixWidth)/float(self.pixHeight)
		self.nearDist = 0.01
		self.farDist = 100.0
		self.xScale = math.tan(self.xfov/2.0)
		self.yScale = math.tan(self.yfov/2.0)

	def initBeamTreeAtPos(self):
		beamWidth = 50*math.pi/180
		orig = self.beamOrigin
		self.scene.Source = orig
		delt = math.tan(beamWidth/2.0)
		#def rotateAroundAxis(P0, axis, theta, V)
		beamVerts = [orig + Vector3D(-delt, -delt, -1), orig + Vector3D(delt, -delt, -1), orig + Vector3D(delt, delt, -1), orig + Vector3D(-delt, delt, -1)]
		#Rotate the beam 45 degrees
		for i in range(0, len(beamVerts)):
			beamVerts[i] = rotateAroundAxis(orig, Vector3D(0, 1, 0), -math.pi/4, beamVerts[i])
		self.beamTree = BeamTree(self.scene.Source, self.meshFaces, 1, [beamVerts])

	def recalculateBeams(self):
		self.beamsToDraw = []
		if self.drawRightBeams:
			path = [self.meshFaces[1], self.meshFaces[3]]
			toAdd = self.beamTree.getBeamsIntersectingFaces(path)
			self.beamsToDraw = self.beamsToDraw + toAdd
			if len(toAdd) > 0:
				beam = toAdd[0]
				face = beam.terminalFace
				P = face.getCentroid()
				self.path = self.beamTree.extractPathToOrigin(beam, P)		
		if self.drawBottomBeams:
			path = [self.meshFaces[4], self.meshFaces[3]]
			self.beamsToDraw = self.beamsToDraw + self.beamTree.getBeamsIntersectingFaces(path)
		if self.drawDirectBeams:
			path = [self.meshFaces[3]]
			newBeams = self.beamTree.getBeamsIntersectingFaces(path)
			self.beamsToDraw = self.beamsToDraw + newBeams

	def __init__(self, filename):
		#GLUT State variables
		self.pixHeight = 800
		self.pixWidth = 800
		self.GLUTmouse = [0, 0]
		self.GLUTButton = [0, 0, 0, 0, 0]
		self.GLUTModifiers = 0
		self.keys = {}
		
		#Camera and projection state variables
		self.camera = MouseSphericalCamera(self.pixWidth, self.pixHeight)
		self.updateCameraVars()
		
		self.scene = EMScene()
		self.scene.Read(filename, False)
		self.meshFaces = []
		for mesh in self.scene.meshes:
			self.meshFaces = self.meshFaces + mesh.faces
		
		self.sceneTransparent = True
		self.beamTrans = 0.3 #Beam transparency
		self.beamOrigin = Point3D(1.5, -1.5, -1.5) + Vector3D(-0.85, 0.95, 1.8288)
		self.initBeamTreeAtPos()
		self.drawRightBeams = False
		self.drawBottomBeams = False
		self.drawDirectBeams = False
		self.beamsToDraw = []
		self.path = []
		
		self.initGL()

	def GLUTResize(self, w, h):
		glViewport(0, 0, w, h)
		self.pixWidth = w
		self.pixHeight = h
		self.camera.pixWidth = w
		self.camera.pixHeight = h
		self.updateCameraVars()
		glutPostRedisplay()

	def GLUTRedraw(self):
		glClearColor(0, 0, 0, 1)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		#First draw the 3D beam scene on the left
		glViewport(0, 0, 800, 800)
		glScissor(0, 0, 800, 800)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(180.0*self.camera.yfov/math.pi, 1.0, 0.01, 100.0)
		
		self.camera.gotoCameraFrame()
		glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
		glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
	
		glEnable(GL_LIGHTING)
		
		if self.sceneTransparent:
			glDisable(GL_DEPTH_TEST)
		else:
			glEnable(GL_DEPTH_TEST)
		self.scene.renderGL()
		
		if self.scene.Source:
			glDisable(GL_LIGHTING)
			glColor3f(1, 0, 0)
			P = self.scene.Source
			quadric = gluNewQuadric()
			glPushMatrix()
			glTranslatef(P.x, P.y, P.z)
			gluSphere(quadric, 0.1, 32, 32)
			glPopMatrix()
		
		for beam in self.beamsToDraw:
			beam.drawBeam()	
		
		glDisable(GL_LIGHTING)
		glColor3f(1, 0, 0)
		glBegin(GL_LINES)
		for i in range(0, len(self.path)-1):
			P0 = self.path[i]
			P1 = self.path[(i+1)]
			glVertex3f(P0.x, P0.y, P0.z)
			glVertex3f(P1.x, P1.y, P1.z)
		glEnd()
		
		glutSwapBuffers()
	
	def handleMouseStuff(self, x, y):
		y = self.pixHeight - y
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
		if key == ' ':
			self.drawBeam = not self.drawBeam
		elif key in ['t', 'T']:
			self.sceneTransparent = not self.sceneTransparent
		#WASD Key control of beam origin
		elif key in ['a', 'A']:
			self.beamOrigin = self.beamOrigin + KEYDELTA*Vector3D(-1, 0, 0)
			self.initBeamTreeAtPos()
			self.recalculateBeams()
		elif key in ['d', 'D']:
			self.beamOrigin = self.beamOrigin + KEYDELTA*Vector3D(1, 0, 0)
			self.initBeamTreeAtPos()
			self.recalculateBeams()
		elif key in ['w', 'W']:
			self.beamOrigin = self.beamOrigin + KEYDELTA*Vector3D(0, 1, 0)
			self.initBeamTreeAtPos()
			self.recalculateBeams()
		elif key in ['s', 'S']:
			self.beamOrigin = self.beamOrigin + KEYDELTA*Vector3D(0, -1, 0)
			self.initBeamTreeAtPos()
			self.recalculateBeams()
		elif key in ['r', 'R']:
			#Draw beams that bounce from the right wall to the back wall
			self.drawRightBeams = not self.drawRightBeams
			self.recalculateBeams()
		elif key in ['b', 'B']:
			self.drawBottomBeams = not self.drawBottomBeams
			self.recalculateBeams()
		elif key in ['o', 'O']:
			self.drawDirectBeams = not self.drawDirectBeams
			self.recalculateBeams()
		elif key in ['p', 'P']:
			#Calculate and display interference pattern
			pattern = self.beamTree.getInterferencePatternOnFace(self.meshFaces[3], self.beamsToDraw, 400, 400, 10e9)
			plt.imshow(pattern)
			plt.show()
		#if key in ['e', 'E']:
		#	self.drawEdges = 1 - self.drawEdges
		#elif key in ['v', 'V']:
		#	self.drawVerts = 1 - self.drawVerts
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
		self.handleMouseStuff(x, y)
		if state == GLUT_DOWN:
			self.GLUTButton[buttonMap[button]] = 1
		else:
			self.GLUTButton[buttonMap[button]] = 0
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
		glutInitWindowSize(self.pixWidth, self.pixHeight)
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
		print "Usage: beamView <scene filepath>"
	else:
		viewer = Viewer(argv[1])
