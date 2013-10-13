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
from sys import argv

DRAW_BACKPROJECTED = True

def splitIntoRGBA(val):
	A = (0xff000000&val)>>24
	R = (0x00ff0000&val)>>16
	G = (0x0000ff00&val)>>8
	B = (0x000000ff&val)
	return [R, G, B, A]

def extractFromRGBA(R, G, B, A):
	A = 0
	return (((A<<24)&0xff000000) | ((R<<16)&0x00ff0000) | ((G<<8)&0x0000ff00) | (B&0x000000ff))

class Viewer(object):
	def updateCameraVars(self):
		self.yfov = self.camera.yfov
		self.xfov = self.yfov*float(self.pixWidth)/float(self.pixHeight)
		self.nearDist = 0.01
		self.farDist = 100.0
		self.xScale = math.tan(self.xfov/2.0)
		self.yScale = math.tan(self.yfov/2.0)
		self.pickingFace = False

	def __init__(self, filename):
		#GLUT State variables
		self.pixHeight = 800
		self.pixWidth = 1600
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
		
		self.selectedFace = None
		self.drawBeam = True
		self.drawChildren = True
		self.beamIndex = 0
		self.sceneTransparent = True
		self.toggleDrawSplits = False
		self.beamTrans = 0.3 #Beam transparency
		self.beamTree = BeamTree(self.scene.Source, self.meshFaces, 2)
		
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
		N = len(self.meshFaces)
		[Rclear, Gclear, Bclear, Aclear] = splitIntoRGBA(N)
		glClearColor(Rclear, Gclear, Bclear, Aclear)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		#First draw the 3D beam scene on the left
		glViewport(0, 0, 800, 800)
		glScissor(0, 0, 800, 800)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		gluPerspective(180.0*self.camera.yfov/math.pi, 1.0, 0.01, 100.0)
		
		self.camera.gotoCameraFrame()
		if self.pickingFace:
			self.pickFace()
		else:
			glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
			glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
		
			glEnable(GL_LIGHTING)
			
			if self.sceneTransparent:
				glDisable(GL_DEPTH_TEST)
			else:
				glEnable(GL_DEPTH_TEST)
			self.scene.renderGL()
			
			if self.selectedFace:
				glDisable(GL_LIGHTING)
				glDisable(GL_DEPTH_TEST)
				glColor3f(1, 0, 0)
				self.selectedFace.drawBorder()
				glEnable(GL_DEPTH_TEST)
			
			if self.scene.Source:
				glDisable(GL_LIGHTING)
				glColor3f(1, 0, 0)
				P = self.scene.Source
				quadric = gluNewQuadric()
				glPushMatrix()
				glTranslatef(P.x, P.y, P.z)
				gluSphere(quadric, 0.1, 32, 32)
				glPopMatrix()
			
			if self.drawBeam:
				beam = self.beamTree.root.children[self.beamIndex]
				beam.drawBeam()
				
				if self.drawChildren:
					for child in beam.children:
						if DRAW_BACKPROJECTED:
							child.drawBackProjected(self.meshFaces)
						child.drawBeam()
						for child2 in child.children:
							child2.drawBeam()
			
			#Next draw the 2D projection scene on the right
			dim = self.pixWidth - 800
			glViewport(800, 0, dim, dim)
			glScissor(800, 0, dim, dim)
			if len(beam.children) > 0:
				beam.children[0].drawProjectedMeshFaces(self.meshFaces, dim, self.toggleDrawSplits)
				
			glutSwapBuffers()
	
	def pickFace(self):
		glDisable(GL_LIGHTING)
		N = len(self.meshFaces)
		for i in range(0, len(self.meshFaces)):
			face = self.meshFaces[i]
			[R, G, B, A] = splitIntoRGBA(i)
			glColor4ub(R, G, B, A)
			face.drawFilled()
		glutSwapBuffers()
		[xdim, ydim] = [self.pixWidth, self.pixHeight]
		pixels = glReadPixelsb(0, 0, 800, 800, GL_RGBA, GL_UNSIGNED_BYTE)
		#im = Image.new("RGB", (800, 800))
		#pix = im.load()
		#for x in range(0, 800):
		#	for y in range(0, 800):
		#		pix[x, y] = (pixels[x][y][0], pixels[x][y][1], pixels[x][y][2])
		#im.save("out.png")
		#print len(pixels)
		#print len(pixels[0])
		[x, y] = [self.GLUTmouse[0], self.GLUTmouse[1]]
		pixel = pixels[y][x]
		faceIndex = extractFromRGBA(pixel[0], pixel[1], pixel[2], pixel[3])
		if faceIndex < N:
			self.selectedFace = self.meshFaces[faceIndex]
		else:
			print "ERROR: No face exists at that location (faceIndex %i)"%faceIndex
			[R, G, B, A] = splitIntoRGBA(faceIndex)
			print "(R, G, B, A) = (%i, %i, %i, %i)"%(R, G, B, A)
		self.pickingFace = False
		glutPostRedisplay()
	
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
		elif key in ['s', 'S']:
			self.toggleDrawSplits = not self.toggleDrawSplits
		elif key in ['t', 'T']:
			self.sceneTransparent = not self.sceneTransparent
		elif key in ['b', 'B']:
			self.beamIndex = (self.beamIndex + 1)%len(self.beamTree.root.children)
		elif key in ['r', 'R']:
			self.drawChildren = not self.drawChildren
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
			if button == GLUT_MIDDLE_BUTTON:
				self.pickingFace = True
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
