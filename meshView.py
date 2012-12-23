from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from PolyMesh import *
from Cameras3D import *
from sys import argv
import random

class Viewer(object):
	def __init__(self, infilename, outfilename = "out.obj"):
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
		self.drawCutPlane = 0
		
		#Camera state variables
		self.outfilename = outfilename
		self.mesh = PolyMesh()
		self.mesh.loadFile(infilename)
		self.meshBBox = self.mesh.getBBox()
		#self.mesh = getBoxMesh(1, 2, 1, Point3D(1, 100, 1), 0.55)
		#print self.mesh
		#self.mesh.truncate(0.2)
		#print self.mesh
		#self.camera = MouseSphericalCamera(self.GLUTwindow_width, self.GLUTwindow_height)
		self.camera = MousePolarCamera(self.GLUTwindow_width, self.GLUTwindow_height)
		self.camera.centerOnMesh(self.mesh)
		random.seed()
		self.cutPlane = None
		self.planeBorder = []
		
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
		farDist = (self.camera.eye - self.meshBBox.getCenter()).Length()*2
		#This is to make sure we can see on the inside
		farDist = max(farDist, self.meshBBox.getDiagLength())
		nearDist = farDist/50.0
		gluPerspective(180.0*self.camera.yfov/M_PI, float(self.GLUTwindow_width)/self.GLUTwindow_height, nearDist, farDist)
		
		#Set up modelview matrix
		self.camera.gotoCameraFrame()	
		glClearColor(0.0, 0.0, 0.0, 0.0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
		glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
		
		glEnable(GL_LIGHTING)
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)

		self.mesh.renderGL(self.drawEdges, self.drawVerts, self.drawNormals, self.planeBorder)
		#self.mesh.renderCCWEdgesDebug()
		
		if self.drawCutPlane:
			t = farDist*self.camera.towards
			r = farDist*(t % self.camera.up)
			dP0 = farDist / 10.0
			#dP0 = 1
			P0 = self.camera.eye - dP0*self.camera.up
			self.cutPlane = getRectMesh(P0 + t + r, P0 + t - r, P0 - t - r, P0 - t + r)
			glDisable(GL_LIGHTING)
			glColor3f(0, 1, 0)
			self.cutPlane.renderGL()
		
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
		if key in ['c', 'C']:
			if len(self.mesh.components) == 0:
				self.mesh.getConnectedComponents()
			counts = self.mesh.getConnectedComponentCounts()
			print "CONNECTED COMPONENTS: "
			total = 0
			for i in range(0, len(counts)):
				print "%i: %i"%(i, counts[i])
				total = total + counts[i]
			print "%i vertices in mesh, sum of components is %i"%(len(self.mesh.vertices), total)
			self.mesh.deleteAllButLargestConnectedComponent()
			print "Number of vertices now: %i"%len(self.mesh.vertices)
			print self.mesh
		elif key in ['e', 'E']:
			self.drawEdges = 1 - self.drawEdges
		elif key in ['h', 'H']:
			self.mesh.fillHoles()
		elif key in ['v', 'V']:
			self.drawVerts = 1 - self.drawVerts
		elif key in ['n', 'N']:
			self.drawNormals = 1 - self.drawNormals
		elif key in ['o', 'O']:
			#Save mesh
			self.mesh.saveFile(self.outfilename, True)
		elif key in ['p', 'P']:
			self.drawCutPlane = 1 - self.drawCutPlane
		elif key in ['r', 'R']:
			#Rotate mesh to align with viewpoint
			[t, u, r] = [self.camera.towards, self.camera.up, self.camera.towards%self.camera.up]
			rotMat = Matrix4([r.x, u.x, -t.x, 0, r.y, u.y, -t.y, 0, r.z, u.z, -t.z, 0, 0, 0, 0, 1])
			rotMat = rotMat.Inverse()
			self.mesh.Transform(rotMat)
			centroid = self.mesh.getCentroid()
			bbox = self.mesh.getBBox()
			minZ = min([v.pos.z for v in self.mesh.vertices])
			dV = -1*centroid
			dV.z = -minZ
			self.mesh.Translate(dV)
			#scale = 1.0/max([bbox.XLen(), bbox.YLen(), bbox.ZLen()])
			#self.mesh.Scale(scale, scale, scale)
			self.camera.centerOnMesh(self.mesh)
			print "minZ = %g"%(min([v.pos.z for v in self.mesh.vertices]))
			print "centroid = %s"%self.mesh.getCentroid()
		elif key in ['s', 'S']:
			self.mesh.sliceBelowPlane(self.cutPlane.faces[0].getPlane())
		elif key in ['t', 'T']:
			print "Triangulating mesh"
			print self.mesh
			self.mesh.minTrianglesRemesh()
			print self.mesh
		elif key in ['d', 'D']:
			#Print depth image
			width = self.GLUTwindow_width
			height = self.GLUTwindow_height
			pixels = glReadPixelsb(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT)
			fout = open('getDepth.m', 'w')
			fout.write("depth = [")
			for y in range(0, height):
				for x in range(0, width):
					fout.write("%s "%pixels[height-y-1][x])
				fout.write(";\n")
			fout.write("];")
		elif key in ['q', 'Q']:
			sys.exit()
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
			self.camera.zoom(-dY)#Want to zoom in as the mouse goes up
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
		print "Usage: meshView <mesh filepath> <out filename (optional)>"
	else:
		if len(argv) < 3:
			#The user has chosen an input name only
			viewer = Viewer(argv[1])
		else:
			#The user has chosen an output filename to save
			viewer = Viewer(argv[1], argv[2])
