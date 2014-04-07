#Based off of http://wiki.wxpython.org/GLCanvas
#Lots of help from http://wiki.wxpython.org/Getting%20Started
from OpenGL.GL import *
import wx
from wx import glcanvas

from Primitives3D import *
from PolyMesh import *
from LaplacianMesh import *
from Geodesics import *
from PointCloud import *
from Cameras3D import *
from ICP import *
from sys import exit, argv
import random
import numpy as np
import scipy.io as sio
from pylab import cm
import os
import subprocess
import math
import time
from sklearn import manifold

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
		#self.camera = MouseSphericalCamera(self.size.x, self.size.y)
		self.camera = MousePolarCamera(self.size.width, self.size.height)
		
		#Main state variables
		self.MousePos = [0, 0]
		self.initiallyResized = False

		self.bbox = BBox3D()
		self.unionbbox = BBox3D()
		random.seed()
		
		#Face mesh variables and manipulation variables
		self.mesh1 = None
		self.mesh1Dist = None
		self.mesh1DistLoaded = False
		self.mesh2 = None
		self.mesh2DistLoaded = False
		self.mesh2Dist = None
		#Holds the transformations of the best iteration in ICP
		self.transformations = []
		self.savingMovie = False
		self.movieIter = 0
		
		self.displayMeshFaces = True
		self.displayMeshEdges = False
		self.displayMeshVertices = False
		self.displayMeshNormals = False
		self.displayPrincipalAxes = False
		self.vertexColors = np.zeros(0)
		
		self.cutPlane = None
		self.displayCutPlane = False
		
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
	
	def initPointCloud(self, pointCloud):
		self.pointCloud = pointCloud
	
	def centerOnMesh1(self, evt):
		if not self.mesh1:
			return
		self.bbox = self.mesh1.getBBox()
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
		self.Refresh()
	
	def centerOnMesh2(self, evt):
		if not self.mesh2:
			return
		self.bbox = self.mesh2.getBBox()
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
		self.Refresh()
		
	def centerOnBoth(self, evt):
		if not self.mesh1 or not self.mesh2:
			return
		self.bbox = self.mesh1.getBBox()
		self.bbox.Union(self.mesh2.getBBox())
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
		self.Refresh()
		
	def MDSMesh1(self, evt):
		if not self.mesh1:
			print "ERROR: Mesh 1 not loaded yet"
			return
		if not self.mesh1DistLoaded:
			print "ERROR: Mesh 1 distance matrix not loaded"
			return
		mds = manifold.MDS(n_components=2, dissimilarity="precomputed", n_jobs=1)
		print "Doing MDS on mesh 1...."
		pos = mds.fit(self.mesh1Dist).embedding_
		print "Finished MDS on mesh 1"
		for i in range(pos.shape[0]):
			self.mesh1.vertices[i].pos = Point3D(pos[i, 0], pos[i, 1], pos[i, 2])
		self.mesh1.needsDisplayUpdate = True
		self.Refresh()
	
	def MDSMesh2(self, evt):
		if not self.mesh2:
			print "ERROR: Mesh 2 not loaded yet"
			return
		if not self.mesh2DistLoaded:
			print "ERROR: Mesh 2 distance matrix not loaded"
			return
		mds = manifold.MDS(n_components=2, dissimilarity="precomputed", n_jobs=1)
		print "Doing MDS on mesh 2..."
		pos = mds.fit(self.mesh2Dist).embedding_
		print "Finished MDS on mesh 2"
		for i in range(pos.shape[0]):
			self.mesh2.vertices[i].pos = Point3D(pos[i, 0], pos[i, 1], pos[i, 2])
		self.mesh2.needsDisplayUpdate = True
		self.Refresh()	
	
	def doGMDS(self, evt):
		if self.mesh1 and self.mesh2:
			print "TODO"
		else:
			print "ERROR: One or both meshes have not been loaded yet"
		self.Refresh()
	
	def displayMeshFacesCheckbox(self, evt):
		self.displayMeshFaces = evt.Checked()
		self.Refresh()

	def displayMeshEdgesCheckbox(self, evt):
		self.displayMeshEdges = evt.Checked()
		self.Refresh()
		
	def displayCutPlaneCheckbox(self, evt):
		self.displayCutPlane = evt.Checked()
		self.Refresh()

	def displayMeshVerticesCheckbox(self, evt):
		self.displayMeshVertices = evt.Checked()
		self.Refresh()

	def displayPrincipalAxesCheckbox(self, evt):
		self.displayPrincipalAxes = evt.Checked()
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
			self.camera.centerOnBBox(self.bbox, math.pi/2, math.pi/2)
			self.initiallyResized = True

	def processPaintEvent(self, event):
		dc = wx.PaintDC(self)
		self.SetCurrent(self.context)
		if not self.GLinitialized:
			self.initGL()
			self.GLinitialized = True
		self.repaint()

	def repaint(self):
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
		
		transform = np.eye(4)
		if len(self.transformations) > 0:
			transform = self.transformations[-1]
		if self.savingMovie:
			transform = self.transformations[self.movieIter]
			self.movieIter = self.movieIter + 1
		if self.mesh1:
			glPushMatrix()
			glMultMatrixd(transform.transpose().flatten())
			self.mesh1.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, None)
			glPopMatrix()
		if self.mesh2:
			self.mesh2.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, None)
		self.SwapBuffers()
		
		if self.savingMovie:
			saveImageGL(self, "ICP%i.png"%self.movieIter)
			if self.movieIter >= len(self.transformations):
				self.savingMovie = False
				#Convert the frames to an OGG video and delete the PNG images
				#system(sprintf('ffmpeg -f image2 -r 4 -i %s%s.tif -r 30 %s.mp4', moviename, '%d', moviename));
				#subprocess.call(["ffmpeg", "-f", "image2" "-r", "2", "-i", "ICP%d.png", "-r", "2", "ICP.ogg"])
				(stdin, stdout, stderr) = os.popen3("ffmpeg -f image2 -r 4 -i ICP%d.png -r 4 ICP.ogg")
				print stdout.readlines()
				print stderr.readlines()
				os.popen3("rm ICP*.png")
			self.Refresh()
	
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
		x, y = evt.GetPosition()
		[lastX, lastY] = self.MousePos
		self.handleMouseStuff(x, y)
		dX = self.MousePos[0] - lastX
		dY = self.MousePos[1] - lastY
		if evt.Dragging():
			if evt.MiddleIsDown():
				self.camera.translate(dX, dY)
			elif evt.RightIsDown():
				self.camera.zoom(-dY)#Want to zoom in as the mouse goes up
			elif evt.LeftIsDown():
				self.camera.orbitLeftRight(dX)
				self.camera.orbitUpDown(dY)
		self.Refresh()

class MeshViewerFrame(wx.Frame):
	(ID_LOADDATASET1, ID_LOADDATASET2, ID_SAVEDATASET, ID_SAVESCREENSHOT) = (1, 2, 3, 4)
	
	def __init__(self, parent, id, title, pos=DEFAULT_POS, size=DEFAULT_SIZE, style=wx.DEFAULT_FRAME_STYLE, name = 'GLWindow'):
		style = style | wx.NO_FULL_REPAINT_ON_RESIZE
		super(MeshViewerFrame, self).__init__(parent, id, title, pos, size, style, name)
		#Initialize the menu
		self.CreateStatusBar()
		
		self.size = size
		self.pos = pos
		print "MeshViewerFrameSize = %s, pos = %s"%(self.size, self.pos)
		
		filemenu = wx.Menu()
		menuOpenMesh1 = filemenu.Append(MeshViewerFrame.ID_LOADDATASET1, "&Load Mesh1","Load a polygon mesh")
		self.Bind(wx.EVT_MENU, self.OnLoadMesh1, menuOpenMesh1)
		menuOpenMesh2 = filemenu.Append(MeshViewerFrame.ID_LOADDATASET2, "&Load Mesh2","Load a polygon mesh")
		self.Bind(wx.EVT_MENU, self.OnLoadMesh2, menuOpenMesh2)
		menuSaveScreenshot = filemenu.Append(MeshViewerFrame.ID_SAVESCREENSHOT, "&Save Screenshot", "Save a screenshot of the GL Canvas")
		self.Bind(wx.EVT_MENU, self.OnSaveScreenshot, menuSaveScreenshot)
		menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")
		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
		
		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		self.glcanvas = MeshViewerCanvas(self)
		
		self.rightPanel = wx.BoxSizer(wx.VERTICAL)
		
		#Buttons to go to a default view
		viewPanel = wx.BoxSizer(wx.HORIZONTAL)
		center1Button = wx.Button(self, -1, "Mesh1")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.centerOnMesh1, center1Button)
		viewPanel.Add(center1Button, 0, wx.EXPAND)
		center2Button = wx.Button(self, -1, "Mesh2")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.centerOnMesh2, center2Button)
		viewPanel.Add(center2Button, 0, wx.EXPAND)
		bothButton = wx.Button(self, -1, "Both")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.centerOnBoth, bothButton)
		viewPanel.Add(bothButton, 0, wx.EXPAND)
		self.rightPanel.Add(wx.StaticText(self, label="Views"), 0, wx.EXPAND)
		self.rightPanel.Add(viewPanel, 0, wx.EXPAND)
		
		#Buttons for MDS
		MDSPanel = wx.BoxSizer(wx.HORIZONTAL)
		MDS1Button = wx.Button(self, -1, "MDS Mesh1")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.MDSMesh1, MDS1Button)
		MDSPanel.Add(MDS1Button, 0, wx.EXPAND)
		MDS2Button = wx.Button(self, -1, "MDS Mesh2")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.MDSMesh2, MDS2Button)
		MDSPanel.Add(MDS2Button, 0, wx.EXPAND)
		self.rightPanel.Add(wx.StaticText(self, label="MDS on Meshes"), 0, wx.EXPAND)
		self.rightPanel.Add(MDSPanel, 0, wx.EXPAND)				
		
		#Checkboxes for displaying data
		self.displayMeshFacesCheckbox = wx.CheckBox(self, label = "Display Mesh Faces")
		self.displayMeshFacesCheckbox.SetValue(True)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshFacesCheckbox, self.displayMeshFacesCheckbox)
		self.rightPanel.Add(self.displayMeshFacesCheckbox, 0, wx.EXPAND)
		self.displayMeshEdgesCheckbox = wx.CheckBox(self, label = "Display Mesh Edges")
		self.displayMeshEdgesCheckbox.SetValue(False)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshEdgesCheckbox, self.displayMeshEdgesCheckbox)
		self.rightPanel.Add(self.displayMeshEdgesCheckbox, 0, wx.EXPAND)
		self.displayMeshVerticesCheckbox = wx.CheckBox(self, label = "Display Mesh Points")
		self.displayMeshVerticesCheckbox.SetValue(False)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshVerticesCheckbox, self.displayMeshVerticesCheckbox)
		self.rightPanel.Add(self.displayMeshVerticesCheckbox)

		
		#Button for doing ICP
		GMDSButton = wx.Button(self, -1, "DO GMDS")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.doGMDS, GMDSButton)
		self.rightPanel.Add(GMDSButton, 0, wx.EXPAND)
		
		#Finally add the two main panels to the sizer		
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.glcanvas, 2, wx.EXPAND)
		self.sizer.Add(self.rightPanel, 0, wx.EXPAND)
		
		self.SetSizer(self.sizer)
		self.Layout()
		self.Show()
	
	def OnLoadMesh1(self, evt):
		print "Mesh 1 menu"
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "OBJ files (*.obj)|*.obj|OFF files (*.off)|*.off", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			print dirname
			self.glcanvas.mesh1 = LaplacianMesh()
			print "Loading mesh %s..."%filename
			self.glcanvas.mesh1.loadFile(filepath)
			print "Finished loading mesh 1\n %s"%self.glcanvas.mesh1
			print self.glcanvas.mesh1
			self.glcanvas.bbox = self.glcanvas.mesh1.getBBox()
			print "Mesh BBox: %s\n"%self.glcanvas.bbox
			self.glcanvas.camera.centerOnBBox(self.glcanvas.bbox, theta = -math.pi/2, phi = math.pi/2)
			#Now try to load in the distance matrix
			fileName, fileExtension = os.path.splitext(filepath)
			matfile = sio.loadmat("%s.mat"%fileName)
			if 'D' in matfile:
				self.glcanvas.mesh1Dist = matfile['D']
				self.glcanvas.mesh1DistLoaded = True
			else:
				print "ERROR: No distance matrix found for mesh 1"
			self.glcanvas.Refresh()
		dlg.Destroy()
		return

	def OnLoadMesh2(self, evt):
		print "Mesh 2 menu"
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "OBJ files (*.obj)|*.obj|OFF files (*.off)|*.off", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			print dirname
			self.glcanvas.mesh2 = LaplacianMesh()
			print "Loading mesh %s..."%filename
			self.glcanvas.mesh2.loadFile(filepath)
			print "Finished loading mesh 2\n %s"%self.glcanvas.mesh2
			print self.glcanvas.mesh2
			print "Mesh BBox: %s\n"%self.glcanvas.mesh2.getBBox()
			self.glcanvas.bbox.Union(self.glcanvas.mesh2.getBBox())
			fileName, fileExtension = os.path.splitext(filepath)
			matfile = sio.loadmat("%s.mat"%fileName)
			if 'D' in matfile:
				self.glcanvas.mesh2Dist = matfile['D']
				self.glcanvas.mesh2DistLoaded = True
			else:
				print "ERROR: No distance matrix found for mesh 2"
			self.glcanvas.Refresh()
		dlg.Destroy()
		return

	def OnSaveScreenshot(self, evt):
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			saveImageGL(self.glcanvas, filepath)
		dlg.Destroy()
		return

	def OnExit(self, evt):
		self.Close(True)
		return

class MeshViewer(object):
	def __init__(self, filename = None, ts = False, sp = "", ra = 0):
		app = wx.App()
		frame = MeshViewerFrame(None, -1, 'MeshViewer')
		frame.Show(True)
		app.MainLoop()
		app.Destroy()

if __name__ == '__main__':
	viewer = MeshViewer()
