from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from PolyMesh import *
from Cameras3D import *

class Viewer(object):
	def __init__(self):
		#GLUT State variables
		self.GLUTwindow_height = 800
		self.GLUTwindow_width = 800
		self.GLUTmouse = [0, 0]
		self.GLUTButton = [0, 0, 0]
		self.GLUTModifiers = 0
		self.keys = {}
		
		#Camera state variables
		#self.camera = Key6DOFCamera(Point3D(-0.230399, -0.0110676, -1.68559), Vector3D(0.135425, 0.0065054, 0.990766), Vector3D(-0.0479816, 0.998848, 0))
		self.camera = Key6DOFCamera(Point3D(-1, 0, 0), Vector3D(1, 0, 0), Vector3D(0, 1, 0))
		#Speeds of the UI controller
		self.dx = 0.1
		self.dTheta = 0.1
		
		self.mesh = PolyMesh()
		self.mesh.loadOffFile("meshes/homer.off")
		print "NVertices = %i, NEdges = %i, NFaces = %i"%(len(self.mesh.vertices), len(self.mesh.edges), len(self.mesh.faces))
		
		self.initGL()

	def GLUTResize(self, w, h):
		glViewport(0, 0, w, h)
		self.GLUTwindow_width = w
		self.GLUTwindow_height = h
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
		
		#print gluErrorString(glGetError())
		
		glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
		glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
		
		glEnable(GL_LIGHTING)
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)

		
		glBegin(GL_TRIANGLES)
		glVertex3f(-1, 0, 0)
		glVertex3f(1, 0, 0)
		glVertex3f(0.5, 1, 0)
		
		glVertex3f(0, -1, 0)
		glVertex3f(0, 1, 0)
		glVertex3f(0, 0, 1)
		glEnd()
		
		self.mesh.renderGL()
		
		glutSwapBuffers()
	
	def handleMouseStuff(self, x, y):
		y = self.GLUTwindow_height - y
		self.GLUTmouse[0] = x
		self.GLUTmouse[1] = y
		self.GLUTmodifiers = glutGetModifiers()
	
	def GLUTKeyboard(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.camera.keys[key] = True
		glutPostRedisplay()
	
	def GLUTKeyboardUp(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.camera.keys[key] = False
		glutPostRedisplay()
	
	def GLUTSpecial(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.camera.keys[key] = True
		glutPostRedisplay()
	
	def GLUTSpecialUp(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.camera.keys[key] = False
		glutPostRedisplay()
		
	def GLUTMouse(self, button, state, x, y):
		self.handleMouseStuff(x, y)
		glutPostRedisplay()

	def GLUTMotion(self, x, y):
		self.handleMouseStuff(x, y)
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
		glutTimerFunc(self.camera.minTimestepLen, self.camera.handleUserInput, 1)
		
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
