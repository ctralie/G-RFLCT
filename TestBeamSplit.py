from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from Primitives3D import *
from Beam3D import *
from sys import argv

class Viewer(object):
	def __init__(self):
		#GLUT State variables
		self.GLUTwindow_height = 800
		self.GLUTwindow_width = 800
		self.GLUTmouse = [0, 0]
		self.GLUTButton = [0, 0, 0, 0, 0]
		self.GLUTModifiers = 0
		self.keys = {}
		
		self.beamPoints = []
		self.polyPoints = []
		self.splitPolys = []
		self.selectingBeam = True
		
		
		self.initGL()

	def GLUTResize(self, w, h):
		glViewport(0, 0, w, h)
		self.GLUTwindow_width = w
		self.GLUTwindow_height = h
		glutPostRedisplay()

	def drawPolyOutline(self, poly):
		glPointSize(8)
		glBegin(GL_POINTS)
		for P in poly:
			glVertex2f(P.x, P.y)
		glEnd()
		glBegin(GL_LINES)
		for i in range(0, len(poly)):
			P1 = poly[i]
			P2 = poly[(i+1)%len(poly)]
			glVertex2f(P1.x, P1.y)
			glVertex2f(P2.x, P2.y)
		glEnd()		

	def GLUTRedraw(self):
		glClearColor(1.0, 1.0, 1.0, 1.0)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glDisable(GL_LIGHTING)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(0, self.GLUTwindow_width, self.GLUTwindow_height, 0, 0, 1)
		glDisable(GL_DEPTH_TEST)
		glMatrixMode (GL_MODELVIEW)
		glLoadIdentity()
		
		#Draw the beam points in red
		glColor3f(1, 0, 0)
		self.drawPolyOutline(self.beamPoints)

		#Draw the polygon points in blue
		glColor3f(0, 0, 1)
		self.drawPolyOutline(self.polyPoints)

		#Draw the individual split polygons
		for i in range(0, len(self.splitPolys)):
			poly = self.splitPolys[i]
			color = float(i)/float(len(self.splitPolys))
			glColor3f(color, color, color)
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
			glBegin(GL_POLYGON)
			for P in poly:
				glVertex2f(P.x, P.y)
			glEnd()
			glPointSize(4)
			glBegin(GL_POINTS)
			for P in poly:
				glVertex2f(P.x, P.y)
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
	
	def enforceCCW(self, poly):
		for i in range(0, len(poly)-2):
			P0 = poly[i]
			P1 = poly[i+1]
			P2 = poly[(i+2)%len(poly)]
			ccw = CCW2D(P0, P1, P2)
			if ccw == 1:
				print "REVERSING"
				poly.reverse()
				break
	
	def GLUTKeyboardUp(self, key, x, y):
		self.handleMouseStuff(x, y)
		self.keys[key] = False
		if key in ['s', 'S']:
			self.selectingBeam = not self.selectingBeam
		elif key in [' ']:
			self.enforceCCW(self.beamPoints)
			self.enforceCCW(self.polyPoints)
			beam = Beam3D(Point3D(0, 0, 0), self.beamPoints)
			beam.frustPoints = []
			for P in self.beamPoints:
				beam.frustPoints.append(P.Copy())
			self.polyPoints = beam.clipToFrustum(self.polyPoints)
			print "Clipped polygon has %i points"%len(self.polyPoints)
			self.splitPolys = splitBeam(beam, self.polyPoints)
			print "Got %i split polygons back: "%len(self.splitPolys),
			for poly in self.splitPolys:
				print "%i "%len(poly),
			print ""
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
		if state == GLUT_DOWN:
			poly = self.beamPoints
			if not self.selectingBeam:
				poly = self.polyPoints
			if button == GLUT_LEFT_BUTTON:
				poly.append(Point3D(x, y, 1))
			elif button == GLUT_RIGHT_BUTTON:
				if len(poly) > 0:
					poly.pop()
		self.handleMouseStuff(x, y)
		glutPostRedisplay()

	def GLUTMotion(self, x, y):
		lastX = self.GLUTmouse[0]
		lastY = self.GLUTmouse[1]
		self.handleMouseStuff(x, y)
		dX = self.GLUTmouse[0] - lastX
		dY = self.GLUTmouse[1] - lastY
	
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
		
		glutMainLoop()

if __name__ == '__main__':
	viewer = Viewer()
