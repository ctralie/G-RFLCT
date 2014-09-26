#Based off of http://wiki.wxpython.org/GLCanvas
#Lots of help from http://wiki.wxpython.org/Getting%20Started
from OpenGL.GL import *
import wx
from wx import glcanvas

from Primitives3D import *
from PolyMesh import *
from LaplacianMesh import *
from Cameras3D import *
from struct import *
from sys import exit, argv
import random
import numpy as np
import os
import math
import time
from time import sleep
from pylab import cm
import matplotlib.pyplot as plt

DEFAULT_SIZE = wx.Size(1200, 800)
DEFAULT_POS = wx.Point(10, 10)
PRINCIPAL_AXES_SCALEFACTOR = 1


def saveImageGL(mvcanvas, filename):
	view = glGetIntegerv(GL_VIEWPORT)
	img = wx.EmptyImage(view[2], view[3] )
	pixels = glReadPixels(0, 0, view[2], view[3], GL_RGB,
		             GL_UNSIGNED_BYTE)
	img.SetData( pixels )
	img = img.Mirror(False)
	img.SaveFile(filename, wx.BITMAP_TYPE_PNG)

def saveImage(canvas, filename):
	s = wx.ScreenDC()
	w, h = canvas.size.Get()
	b = wx.EmptyBitmap(w, h)
	m = wx.MemoryDCFromDC(s)
	m.SelectObject(b)
	m.Blit(0, 0, w, h, s, 70, 0)
	m.SelectObject(wx.NullBitmap)
	b.SaveFile(filename, wx.BITMAP_TYPE_PNG)

class MeshViewerCanvas(glcanvas.GLCanvas):
	def __init__(self, parent):
		attribs = (glcanvas.WX_GL_RGBA, glcanvas.WX_GL_DOUBLEBUFFER, glcanvas.WX_GL_DEPTH_SIZE, 24)
		glcanvas.GLCanvas.__init__(self, parent, -1, attribList = attribs)	
		self.context = glcanvas.GLContext(self)
		
		self.parent = parent
		#Camera state variables
		self.size = self.GetClientSize()
		self.camera = MousePolarCamera(self.size.x, self.size.y)
		#self.camera = MousePolarCamera(self.size.width, self.size.height)
		
		#Main state variables
		self.MousePos = [0, 0]
		self.initiallyResized = False

		self.bbox = BBox3D()
		self.unionbbox = BBox3D()
		random.seed()
		
		#State variables for saving screenshots
		self.filepath = None
		self.rotfilePrefix = "Rotation"
		self.rotstartAngle = -50
		self.rotendAngle = 50
		self.rotangleInc = 5
		self.rotAngle = 0
		self.zCenter = 0
		
		
		#Face mesh variables and manipulation variables
		self.mesh = None
		self.meshCentroid = None
		self.meshPrincipalAxes = None
		self.displayMeshFaces = True
		self.displayMeshEdges = False
		self.displayMeshVertices = False
		self.displayMeshNormals = False
		self.displayPrincipalAxes = False
		self.useLighting = False
		self.useTexture = True
		self.vertexColors = np.zeros(0)
		
		self.GLinitialized = False
		#GL-related events
		wx.EVT_ERASE_BACKGROUND(self, self.processEraseBackgroundEvent)
		wx.EVT_SIZE(self, self.processSizeEvent)
		wx.EVT_PAINT(self, self.processPaintEvent)
		#Mouse Events
		wx.EVT_LEFT_DOWN(self, self.MouseDown)
		wx.EVT_LEFT_UP(self, self.MouseUp)
		wx.EVT_RIGHT_DOWN(self, self.MouseDown)
		wx.EVT_RIGHT_UP(self, self.MouseUp)
		wx.EVT_MIDDLE_DOWN(self, self.MouseDown)
		wx.EVT_MIDDLE_UP(self, self.MouseUp)
		wx.EVT_MOTION(self, self.MouseMotion)		
		#self.initGL()
	
	def LoadMesh(self):
		self.mesh = LaplacianMesh()
		self.mesh.loadFile(self.filepath)
		print "Finished loading mesh\n %s"%self.mesh
		print self.mesh
		self.bbox = self.mesh.getBBox()
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
		self.zCenter = (self.bbox.zmax + self.bbox.zmin) / 2.0
		self.Refresh()
	
	def processEraseBackgroundEvent(self, event): pass #avoid flashing on MSW.

	def processSizeEvent(self, event):
		self.size = self.GetClientSize()
		self.SetCurrent(self.context)
		glViewport(0, 0, self.size.width, self.size.height)
		if not self.initiallyResized:
			#The canvas gets resized once on initialization so the camera needs
			#to be updated accordingly at that point
			self.camera = MousePolarCamera(self.size.width, self.size.height)
			self.camera.centerOnBBox(self.bbox, -math.pi/2, math.pi/2)
			self.initiallyResized = True

	def processPaintEvent(self, event):
		dc = wx.PaintDC(self)
		self.SetCurrent(self.context)
		if not self.GLinitialized:
			self.initGL()
			self.GLinitialized = True
		self.repaint()

	def repaint(self):
		if not self.mesh:
			self.LoadMesh()
			return
		#Set up projection matrix
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		farDist = (self.camera.eye - self.bbox.getCenter()).Length()*2
		#This is to make sure we can see on the inside
		farDist = max(farDist, self.unionbbox.getDiagLength()*2)
		nearDist = farDist/50.0
		gluPerspective(180.0*self.camera.yfov/M_PI, float(self.size.x)/self.size.y, nearDist, farDist)
		
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

		glTranslatef(0, 0, self.zCenter)
		glRotatef(self.rotAngle, 0, 1, 0)
		glTranslatef(0, 0, -self.zCenter)
		self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, self.useLighting, self.useTexture)
		saveImageGL(self, "%s%i.png"%(self.rotfilePrefix, self.rotAngle))
		self.rotAngle = self.rotAngle + self.rotangleInc
		if self.rotAngle > self.rotendAngle:
			exit(0)
		self.SwapBuffers()
		
		self.Refresh()
	
		self.SwapBuffers()
	
	def initGL(self):		
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glEnable(GL_LIGHT0)
		glLightfv(GL_LIGHT1, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
		glEnable(GL_LIGHT1)
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHTING)
		glEnable(GL_DEPTH_TEST)

	def handleMouseStuff(self, x, y):
		#Invert y from what the window manager says
		y = self.size.height - y
		self.MousePos = [x, y]

	def MouseDown(self, evt):
		state = wx.GetMouseState()
		x, y = evt.GetPosition()
		self.CaptureMouse()
		self.handleMouseStuff(x, y)
		self.Refresh()
	
	def MouseUp(self, evt):
		x, y = evt.GetPosition()
		self.handleMouseStuff(x, y)
		self.ReleaseMouse()
		self.Refresh()

	def MouseMotion(self, evt):
		state = wx.GetMouseState()
		x, y = evt.GetPosition()
		[lastX, lastY] = self.MousePos
		self.handleMouseStuff(x, y)
		dX = self.MousePos[0] - lastX
		dY = self.MousePos[1] - lastY
		self.Refresh()

class MeshViewerFrame(wx.Frame):
	
	def __init__(self, parent, id, title, filepath, outprefix, startAngle, angleInc, endAngle, pos=DEFAULT_POS, size=DEFAULT_SIZE, style=wx.DEFAULT_FRAME_STYLE, name = 'GLWindow'):
		style = style | wx.NO_FULL_REPAINT_ON_RESIZE
		super(MeshViewerFrame, self).__init__(parent, id, title, pos, size, style, name)
		#Initialize the menu
		self.CreateStatusBar()
		
		self.size = size
		self.pos = pos
		self.glcanvas = MeshViewerCanvas(self)
		self.glcanvas.filepath = filepath
		self.glcanvas.rotstartAngle = startAngle
		self.glcanvas.rotangleInc = angleInc
		self.glcanvas.rotendAngle = endAngle
		self.glcanvas.rotAngle = startAngle
		self.glcanvas.rotfilePrefix = outprefix
			
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.glcanvas, 2, wx.EXPAND)
		
		self.SetSizer(self.sizer)
		self.Layout()
		self.Show()

class MeshViewer(object):
	def __init__(self, fpath, oprefix, sAngle, aInc, eAngle):
		app = wx.App()
		frame = MeshViewerFrame(None, -1, 'MeshViewer', filepath = fpath, outprefix = oprefix, startAngle = sAngle, angleInc = aInc, endAngle = eAngle)
		frame.Show(True)
		app.MainLoop()
		app.Destroy()

if __name__ == '__main__':
	if len(argv) < 6:
		print "Usage: python meshRotate.py <mesh name> <fileoutprefix> <startAngle> <angleIncrement> <endAngle>"
		exit(0)
	filepath = argv[1]
	outprefix = argv[2]
	startAngle = float(argv[3])
	angleInc = float(argv[4])
	endAngle = float(argv[5])
	viewer = MeshViewer(filepath, outprefix, startAngle, angleInc, endAngle)
