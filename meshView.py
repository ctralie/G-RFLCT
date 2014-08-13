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
from struct import *
from PRST import PRST
from sys import exit, argv
import random
import numpy as np
import scipy.io as sio
from pylab import cm
import os
import math
import time
from time import sleep
from pylab import cm
import matplotlib.pyplot as plt

DEFAULT_SIZE = wx.Size(1200, 800)
DEFAULT_POS = wx.Point(10, 10)
PRINCIPAL_AXES_SCALEFACTOR = 1


#GUI States
(STATE_NORMAL, STATE_SAVEROTATIONS, STATE_SAVELIGHTING, STATE_CHOOSELAPLACEVERTICES) = (0, 1, 2, 3)
SUBSTATE_NONE = -1

#Laplacian substates
(CHOOSELAPLACE_WAITING, CHOOSELAPLACE_PICKVERTEX) = (0, 1)

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
		
		self.GUIState = STATE_NORMAL
		self.GUISubstate = SUBSTATE_NONE
		
		self.parent = parent
		#Camera state variables
		self.size = self.GetClientSize()
		self.camera = MouseSphericalCamera(self.size.x, self.size.y)
		#self.camera = MousePolarCamera(self.size.width, self.size.height)
		
		#Main state variables
		self.MousePos = [0, 0]
		self.initiallyResized = False

		self.bbox = BBox3D()
		self.unionbbox = BBox3D()
		random.seed()
		
		#State variables for saving screenshots
		self.rotAngle = 0
		self.lightPhi = 0
		self.lightTheta = 0
		self.lightIter = 0
		
		#State variables for laplacian mesh operations
		self.laplacianConstraints = {} #Elements will be key-value pairs (idx, Point3D(new position))
		self.laplaceCurrentIdx = -1
		self.eigvalues = np.zeros(0)
		self.eigvectors = np.zeros(0)
		
		#Face mesh variables and manipulation variables
		self.mesh = None
		self.meshCentroid = None
		self.meshPrincipalAxes = None
		self.displayMeshFaces = True
		self.displayMeshEdges = False
		self.displayMeshVertices = False
		self.displayMeshNormals = False
		self.displayPrincipalAxes = False
		self.useLighting = True
		self.vertexColors = np.zeros(0)

		self.PRSTPlaneMesh = None
		self.PRSTPointPairs = []
		
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
	
	def invalidateMesh(self):
		self.eigvalues = np.zeros(0)
		self.eigvectors = np.zeros(0)
		self.needsDisplayUpdate = True
		self.needsIndexDisplayUpdate = True
	
	def initPointCloud(self, pointCloud):
		self.pointCloud = pointCloud
	
	def viewFromFront(self, evt):
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
		self.Refresh()
	
	def viewFromTop(self, evt):
		self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = 0)
		self.Refresh()
	
	def viewFromSide(self, evt):
		self.camera.centerOnBBox(self.bbox, theta = -math.pi, phi = math.pi/2)
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

	def useLightingCheckbox(self, evt):
		self.useLighting = evt.Checked()
		self.needsDisplayUpdate = True
		self.Refresh()
	
	def CutWithPlane(self, evt):
		if self.cutPlane:
			self.mesh.sliceBelowPlane(self.cutPlane, False)
			self.mesh.starTriangulate() #TODO: This is a patch to deal with "non-planar faces" added
			self.Refresh()

	def doPRST(self, evt):
		if self.mesh:
			(plane, self.PRSTPointPairs) = PRST(self.mesh)
			R = self.mesh.getBBox().getDiagLength()
			P0 = plane.P0
			N = plane.N
			#Quick-n-dirty find directions orthogonal to normal
			right = Vector3D(-N.y, N.x, 0)
			if right.Length() == 0:
				right = Vector3D(0, 0, 1)
			right.normalize()
			up = right % N
			right = R*right
			up = R*up
			self.PRSTPlaneMesh = getRectMesh(P0 - right - up, P0 + right - up, P0 + right + up, P0 - right + up)
			self.Refresh()

	def FlipAcrossPlane(self, evt):
		if self.cutPlane:
			self.mesh.flipAcrossPlane(self.cutPlane)
			self.Refresh()

	def saveAutoRotatingScreenshots(self, evt):
		dlg = wx.TextEntryDialog(None,'File Prefix Name','Saving Rotation Sweep', 'Rotation')
		if dlg.ShowModal() == wx.ID_OK:
			self.rotatingScreenshotsPrefix = dlg.GetValue()
			self.GUIState = STATE_SAVEROTATIONS
			self.rotAngle = -80
			self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
			self.zCenter = (self.bbox.zmax + self.bbox.zmin) / 2.0
			self.Refresh()
		dlg.Destroy()

	def saveLightSweepScreenshots(self, evt):
		dlg = wx.TextEntryDialog(None,'File Prefix Name','Saving Lighting Sweep', 'Lighting')
		if dlg.ShowModal() == wx.ID_OK:
			self.lightingScreenshotsPrefix = dlg.GetValue()
			self.GUIState = STATE_SAVELIGHTING
			self.lightPhi = -math.pi/2
			self.lightTheta = 0
			self.lightIter = 0
			self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
			self.zCenter = (self.bbox.zmax + self.bbox.zmin) / 2.0
			self.Refresh()
		dlg.Destroy()

	def deleteConnectedComponents(self, evt):
		if self.mesh:
			self.mesh.deleteAllButLargestConnectedComponent()

	def splitFaces(self, evt):
		if self.mesh:
			self.mesh.splitFaces()
			self.invalidateMesh()
		self.Refresh()

	def FillHoles(self, evt):
		if self.mesh:
			self.mesh.fillHoles()
			self.invalidateMesh()
		self.Refresh()
	
	def Truncate(self, evt):
		if self.mesh:
			self.mesh.truncate(0.5)
			self.invalidateMesh()
		self.Refresh()
	
	def ComputeGeodesicDistances(self, evt):
		if not self.mesh:
			print "ERROR: Haven't loaded mesh yet"
			return
		D = getGeodesicDistancesFMM(self.mesh)
		D = D[0, :]
		minD = min(D)
		maxD = max(D)
		print "Finished computing geodesic distances"
		print "minD = %g, maxD = %g"%(minD, maxD)
		N = D.shape[0]
		cmConvert = cm.get_cmap('jet')
		self.vertexColors = np.zeros((N, 3))
		for i in range(0, N):
			self.vertexColors[i, :] = cmConvert((D[i] - minD)/(maxD - minD))[0:3]
		self.Refresh()
	
	def InterpolateRandomColors(self, evt):
		#constraints: [[(i1, w1), (i2, w2), ..., (in, wn)], ...], M constraints
		#deltaCoords: NxY numpy array, where Y is the dimension
		#g: MxY numpy array representing values of the constraints, where Y is the dimension
		#and M is the number of constraints
		N = len(self.mesh.vertices)
		constraints = [ [(0, 1)], [(N-1, 1)], [(int(math.floor(N/2)), 1) ] ]
		deltaCoords = np.zeros((N, 3))
		#Make the first vertex blue and the last vertex red
		g = np.array([ [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0] ])
		colors = self.mesh.solveFunctionWithConstraints(constraints, deltaCoords, g)
		for i in range(0, N):
			self.mesh.vertices[i].color = [a for a in colors[i]]
		self.mesh.needsDisplayUpdate = True
		self.mesh.needsIndexDisplayUpdate = True
		self.Refresh()
	
	def doLaplacianMeshSelectVertices(self, evt):
		if self.mesh:
			self.GUIState = STATE_CHOOSELAPLACEVERTICES
			self.GUISubstate = CHOOSELAPLACE_WAITING
			self.Refresh()
	
	def clearLaplacianMeshSelection(self, evt):
		self.laplacianConstraints.clear()
		self.Refresh()
	
	def doLaplacianSolveWithConstraints(self, evt):
		if self.mesh:
			self.mesh.solveVertexPositionsWithConstraints(self.laplacianConstraints)
		self.invalidateMesh()
		self.Refresh()

	def doLaplacianMembraneWithConstraints(self, evt):
		if self.mesh:
			self.mesh.createMembraneSurface(self.laplacianConstraints)
		self.invalidateMesh()
		self.Refresh()

	def getHeatKernelSignature(self, evt):
		dlg = wx.TextEntryDialog(None,'Time Parameter','Choose Time Parameter', '1')
		if dlg.ShowModal() == wx.ID_OK:
			t = float(dlg.GetValue())
			k = min(100, len(self.mesh.vertices)-1)
			cmConvert = cm.get_cmap("jet")
			(x, self.eigvalues, self.eigvectors) = self.mesh.getHKS(k, t, self.eigvalues, self.eigvectors)
			plt.subplot(1, 2, 1)
			plt.plot(self.eigvalues)
			plt.title('Laplace Beltrami Eigenvalues')
			plt.subplot(1, 2, 2)
			plt.plot(np.exp(-t*self.eigvalues))
			plt.title('Exponential Factors')
			plt.show()
			x = x/np.max(x)
			for i in range(len(self.mesh.vertices)):
				self.mesh.vertices[i].color = cmConvert(x[i])
			self.mesh.needsDisplayUpdate = True
			self.Refresh()
		dlg.Destroy()

	def getHeatFlow(self, evt):
		dlg = wx.TextEntryDialog(None,'Time Parameter','Choose Time Parameter', '1')
		if dlg.ShowModal() == wx.ID_OK:
			indices = [a for a in self.laplacianConstraints]
			t = float(dlg.GetValue())
			k = min(100, len(self.mesh.vertices)-1)
			cmConvert = cm.get_cmap("jet")
			(x, self.eigvalues, self.eigvectors) = self.mesh.getHeatFlowFromPoints(indices, k, t, self.eigvalues, self.eigvectors)
			x = x/np.max(x)
			for i in range(len(self.mesh.vertices)):
				self.mesh.vertices[i].color = cmConvert(x[i])
			self.mesh.needsDisplayUpdate = True
			self.Refresh()
		dlg.Destroy()
	
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

	def drawMeshStandard(self):
		if self.mesh:
			glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
			glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
			glEnable(GL_LIGHTING)
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0]);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
			glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
			self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, self.useLighting, None)
			#self.mesh.renderGLTriBuffers()
			if self.displayPrincipalAxes:
				(Axis1, Axis2, Axis3, maxProj, minProj, axes)	= self.meshPrincipalAxes
				C = self.meshCentroid
				glDisable(GL_LIGHTING)
				glLineWidth(5)
				
				glColor3f(1, 0, 0)
				glBegin(GL_LINES)
				glVertex3f(C.x, C.y, C.z)
				A1 = C + maxProj[0]*PRINCIPAL_AXES_SCALEFACTOR*Axis1
				glVertex3f(A1.x, A1.y, A1.z)
				glEnd()
				
				glColor3f(0, 0, 1)
				glBegin(GL_LINES)
				glVertex3f(C.x, C.y, C.z)
				A2 = C + maxProj[1]*PRINCIPAL_AXES_SCALEFACTOR*Axis2
				glVertex3f(A2.x, A2.y, A2.z)
				glEnd()				
				
				glColor3f(0, 1, 0)
				glBegin(GL_LINES)
				glVertex3f(C.x, C.y, C.z)
				A3 = C + maxProj[2]*PRINCIPAL_AXES_SCALEFACTOR*Axis3
				glVertex3f(A3.x, A3.y, A3.z)
				glEnd()				
				
				glEnable(GL_LIGHTING)

		if self.PRSTPlaneMesh:
			self.PRSTPlaneMesh.renderGL()
			glDisable(GL_LIGHTING)
			glColor3f(1, 0, 0)
			glPointSize(10)
			glBegin(GL_POINTS)
			for pair in self.PRSTPointPairs:
				P = pair[0]
				glVertex3f(P.x, P.y, P.z)
			glColor3f(0, 1, 0)
			for pair in self.PRSTPointPairs:
				P = pair[1]
				glVertex3f(P.x, P.y, P.z)
			glEnd()
			glColor3f(0, 0, 1)
			glLineWidth(5)
			glBegin(GL_LINES)
			for pair in self.PRSTPointPairs:
				glVertex3f(pair[0].x, pair[0].y, pair[0].z)
				glVertex3f(pair[1].x, pair[1].y, pair[1].z)
			glEnd()
		
		if self.displayCutPlane:
			#TODO: Many invalid states are possible if laplacian constraints aren't
			#cleared before other things happen
			self.laplacianConstraints.clear()
			farDist = (self.camera.eye - self.bbox.getCenter()).Length()*2
			#This is to make sure we can see on the inside
			farDist = max(farDist, self.unionbbox.getDiagLength()*2)
			t = farDist*self.camera.towards
			r = t % self.camera.up
			u = farDist*self.camera.up
			dP0 = farDist / 10.0
			P0 = self.camera.eye - (dP0/farDist/10.0)*r
			cutPlaneMesh = getRectMesh(P0 + t + u, P0 + t - u, P0 - t - u, P0 - t + u)
			glDisable(GL_LIGHTING)
			glColor3f(0, 1, 0)
			cutPlaneMesh.renderGL(lightingOn = False)
			self.cutPlane = Plane3D(P0, r)
			print "%s, %s"%(P0, r)

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
		
		if self.GUIState == STATE_NORMAL:
			self.drawMeshStandard()
		elif self.GUIState == STATE_SAVELIGHTING:
			lightPos = np.array([3.0, 4.0, 5.0])
			Phi = self.lightPhi
			RotMatrix = np.array([[1, 0, 0], [0, math.cos(Phi), -math.sin(Phi)], [0, math.sin(Phi), math.cos(Phi)] ])
			lightPos = RotMatrix.dot(lightPos)
			glLightfv(GL_LIGHT0, GL_POSITION, [lightPos[0], lightPos[1], lightPos[2], 0.0]);
			glEnable(GL_LIGHTING)
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0]);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
			glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
			self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, self.useLighting, None)
			saveImageGL(self, "%s%i.png"%(self.lightingScreenshotsPrefix, self.lightIter))
			self.lightPhi = self.lightPhi + 0.05
			self.lightIter = self.lightIter + 1
			if self.lightPhi > math.pi/2:
				self.GUIState = STATE_NORMAL
			self.SwapBuffers()
			self.Refresh()
		elif self.GUIState == STATE_SAVEROTATIONS:
			glLightfv(GL_LIGHT0, GL_POSITION, [3.0, 4.0, 5.0, 0.0]);
			glLightfv(GL_LIGHT1, GL_POSITION,  [-3.0, -2.0, -3.0, 0.0]);
			glEnable(GL_LIGHTING)
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0]);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
			glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
			self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, self.useLighting, None)
			glTranslatef(0, 0, self.zCenter)
			glRotatef(self.rotAngle, 0, 1, 0)
			glTranslatef(0, 0, -self.zCenter)
			self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshNormals, self.displayMeshFaces, self.useLighting, None)
			saveImageGL(self, "%s%i.png"%(self.rotatingScreenshotsPrefix, self.rotAngle))
			self.rotAngle = self.rotAngle + 5
			if self.rotAngle > 80:
				self.GUIState = STATE_NORMAL
			self.SwapBuffers()
			self.Refresh()
		elif self.GUIState == STATE_CHOOSELAPLACEVERTICES:
			if self.GUISubstate == CHOOSELAPLACE_WAITING:
				self.drawMeshStandard()
				glDisable(GL_LIGHTING)
				glPointSize(10)
				glBegin(GL_POINTS)
				for idx in self.laplacianConstraints:
					if idx == self.laplaceCurrentIdx:
						glColor3f(1, 0, 0)
					else:
						glColor3f(0, 0, 1)
					P = self.mesh.vertices[idx].pos
					glVertex3f(P.x, P.y, P.z)
					P = self.laplacianConstraints[idx]
					glVertex3f(P.x, P.y, P.z)
				glEnd()
				glColor3f(1, 1, 0)
				glBegin(GL_LINES)
				for idx in self.laplacianConstraints:
					P1 = self.mesh.vertices[idx].pos
					P2 = self.laplacianConstraints[idx]
					glVertex3f(P1.x, P1.y, P1.z)
					glVertex3f(P2.x, P2.y, P2.z)
				glEnd()
				
			elif self.GUISubstate == CHOOSELAPLACE_PICKVERTEX:
				glClearColor(0.0, 0.0, 0.0, 0.0)
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
				self.mesh.renderGLIndices()
				pixel = glReadPixels(self.MousePos[0], self.MousePos[1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE)
				[R, G, B, A] = [int(pixel.encode("hex")[i*2:(i+1)*2], 16) for i in range(4)]
				idx = extractFromRGBA(R, G, B, 0) - 1
				if idx >= 0 and idx < len(self.mesh.vertices):
					if idx in self.laplacianConstraints:
						#De-select if it's already selected
						self.laplaceCurrentIdx = -1
						self.laplacianConstraints.pop(idx, None)
					else:
						self.laplacianConstraints[idx] = self.mesh.vertices[idx].pos.Copy()
						self.laplaceCurrentIdx = idx
				self.GUISubstate = CHOOSELAPLACE_WAITING
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
		if self.GUIState == STATE_CHOOSELAPLACEVERTICES:
			if state.ShiftDown():
				#Pick vertex for laplacian mesh constraints
				self.GUISubstate = CHOOSELAPLACE_PICKVERTEX
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
		if evt.Dragging():
			idx = self.laplaceCurrentIdx
			if self.GUIState == STATE_CHOOSELAPLACEVERTICES and state.ControlDown() and self.laplaceCurrentIdx in self.laplacianConstraints:
				#Move up laplacian mesh constraint based on where the user drags
				#the mouse
				t = self.camera.towards
				u = self.camera.up
				r = t % u
				P0 = self.laplacianConstraints[idx]
				#Construct a plane going through the point which is parallel to the
				#viewing plane
				plane = Plane3D(P0, t)
				#Construct a ray through the pixel where the user is clicking
				tanFOV = math.tan(self.camera.yfov/2)
				scaleX = tanFOV*(self.MousePos[0] - self.size.x/2)/(self.size.x/2)
				scaleY = tanFOV*(self.MousePos[1] - self.size.y/2)/(self.size.y/2)
				V = t + scaleX*r + scaleY*u
				ray = Ray3D(self.camera.eye, V)
				self.laplacianConstraints[idx] = ray.intersectPlane(plane)[1]
			else:
				#Translate/rotate shape
				if evt.MiddleIsDown():
					self.camera.translate(dX, dY)
				elif evt.RightIsDown():
					self.camera.zoom(-dY)#Want to zoom in as the mouse goes up
				elif evt.LeftIsDown():
					self.camera.orbitLeftRight(dX)
					self.camera.orbitUpDown(dY)
		self.Refresh()

class MeshViewerFrame(wx.Frame):
	(ID_LOADDATASET, ID_SAVEDATASET, ID_SAVEDATASETMETERS, ID_SAVESCREENSHOT, ID_CONNECTEDCOMPONENTS, ID_SPLITFACES, ID_TRUNCATE, ID_FILLHOLES, ID_GEODESICDISTANCES, ID_PRST, ID_INTERPOLATECOLORS, ID_SAVEROTATINGSCREENSOTS, ID_SAVELIGHTINGSCREENSHOTS, ID_SELECTLAPLACEVERTICES, ID_CLEARLAPLACEVERTICES, ID_SOLVEWITHCONSTRAINTS, ID_MEMBRANEWITHCONSTRAINTS, ID_GETHKS, ID_GETHEATFLOW) = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
	
	def __init__(self, parent, id, title, pos=DEFAULT_POS, size=DEFAULT_SIZE, style=wx.DEFAULT_FRAME_STYLE, name = 'GLWindow'):
		style = style | wx.NO_FULL_REPAINT_ON_RESIZE
		super(MeshViewerFrame, self).__init__(parent, id, title, pos, size, style, name)
		#Initialize the menu
		self.CreateStatusBar()
		
		self.size = size
		self.pos = pos
		print "MeshViewerFrameSize = %s, pos = %s"%(self.size, self.pos)
		self.glcanvas = MeshViewerCanvas(self)
		
		#####File menu
		filemenu = wx.Menu()
		menuOpenMesh = filemenu.Append(MeshViewerFrame.ID_LOADDATASET, "&Load Mesh","Load a polygon mesh")
		self.Bind(wx.EVT_MENU, self.OnLoadMesh, menuOpenMesh)
		menuSaveMesh = filemenu.Append(MeshViewerFrame.ID_SAVEDATASET, "&Save Mesh", "Save the edited polygon mesh")
		self.Bind(wx.EVT_MENU, self.OnSaveMesh, menuSaveMesh)
		menuSaveMeshMeters = filemenu.Append(MeshViewerFrame.ID_SAVEDATASETMETERS, "&Save Mesh in Meters", "Save the edited polygon mesh; convert from millimeters to meters")
		self.Bind(wx.EVT_MENU, self.OnSaveMeshMeters, menuSaveMeshMeters)
		menuSaveScreenshot = filemenu.Append(MeshViewerFrame.ID_SAVESCREENSHOT, "&Save Screenshot", "Save a screenshot of the GL Canvas")
		self.Bind(wx.EVT_MENU, self.OnSaveScreenshot, menuSaveScreenshot)
		
		menuSaveRotating = filemenu.Append(MeshViewerFrame.ID_SAVEROTATINGSCREENSOTS, "&Save Rotating Screenshots", "Save Auto Rotating Screenshots")
		self.Bind(wx.EVT_MENU, self.glcanvas.saveAutoRotatingScreenshots, menuSaveRotating)
		menuSaveLighting = filemenu.Append(MeshViewerFrame.ID_SAVELIGHTINGSCREENSHOTS, "&Save Lighting Sweep Screenshots", "Save Lighting Sweep Screenshots")
		self.Bind(wx.EVT_MENU, self.glcanvas.saveLightSweepScreenshots, menuSaveLighting)	
		
		menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")
		self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
		
		#####Operations menu
		operationsMenu = wx.Menu()
		#Menu option for deleting all but largest connected component
		menuConnectedComponents = operationsMenu.Append(MeshViewerFrame.ID_CONNECTEDCOMPONENTS, "&Delete all but largest connected component", "Delete all but largest connected component")
		self.Bind(wx.EVT_MENU, self.glcanvas.deleteConnectedComponents, menuConnectedComponents)

		#Menu option for splitting faces
		menuSplitFaces = operationsMenu.Append(MeshViewerFrame.ID_SPLITFACES, "&Split Faces", "Split Faces")
		self.Bind(wx.EVT_MENU, self.glcanvas.splitFaces, menuSplitFaces)
		
		#Menu option for truncating mesh
		menuTruncate = operationsMenu.Append(MeshViewerFrame.ID_TRUNCATE, "&Truncate", "Truncate")
		self.Bind(wx.EVT_MENU, self.glcanvas.Truncate, menuTruncate)
		
		#Menu option for filling holes
		menuFillHoles = operationsMenu.Append(MeshViewerFrame.ID_FILLHOLES, "&Fill Holes", "Fill Holes")
		self.Bind(wx.EVT_MENU, self.glcanvas.FillHoles, menuFillHoles)
		
		#Menu option for computing geodesic distance
		menuGeodesicDistances = operationsMenu.Append(MeshViewerFrame.ID_GEODESICDISTANCES, "&Compute Geodesic Distances", "Compute Geodesic Distances")
		self.Bind(wx.EVT_MENU, self.glcanvas.ComputeGeodesicDistances, menuGeodesicDistances)
		
		#Menu option for interpolating colors
		menuInterpolateColors = operationsMenu.Append(MeshViewerFrame.ID_INTERPOLATECOLORS, "&Interpolate Random Colors", "Interpolate Random Colors")
		self.Bind(wx.EVT_MENU, self.glcanvas.InterpolateRandomColors, menuInterpolateColors)		
		
		#Max planar symmetry
		menuPRST = operationsMenu.Append(MeshViewerFrame.ID_PRST, "&Max Planar Symmetry", "Max Planar Symmetry")
		self.Bind(wx.EVT_MENU, self.glcanvas.doPRST, menuPRST)
		
		#####Laplacian Mesh Menu
		laplacianMenu = wx.Menu()
		menuSelectLaplaceVertices = laplacianMenu.Append(MeshViewerFrame.ID_SELECTLAPLACEVERTICES, "&Select Laplace Vertices", "Select Laplace Vertices")
		self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianMeshSelectVertices, menuSelectLaplaceVertices)
		
		menuClearLaplaceVertices = laplacianMenu.Append(MeshViewerFrame.ID_CLEARLAPLACEVERTICES, "&Clear vertex selection", "Clear Vertex Selection")
		self.Bind(wx.EVT_MENU, self.glcanvas.clearLaplacianMeshSelection, menuClearLaplaceVertices)
		
		menuSolveWithConstraints = laplacianMenu.Append(MeshViewerFrame.ID_SOLVEWITHCONSTRAINTS, "&Solve with Constraints", "Solve with Constraints")
		self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianSolveWithConstraints, menuSolveWithConstraints)

		menuMembraneWithConstraints = laplacianMenu.Append(MeshViewerFrame.ID_MEMBRANEWITHCONSTRAINTS, "&Membrane Surface with Constraints", "Membrane Surface with Constraints")
		self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianMembraneWithConstraints, menuMembraneWithConstraints)
		
		#Get heat kernel signature
		menuHKS = laplacianMenu.Append(MeshViewerFrame.ID_GETHKS, "&Get Heat Kernel Signature", "Get Heat Kernel Signature")
		self.Bind(wx.EVT_MENU, self.glcanvas.getHeatKernelSignature, menuHKS)
		
		menuHeatFlow = laplacianMenu.Append(MeshViewerFrame.ID_GETHEATFLOW, "&Get Heat Flow from Vertices", "Get Heat Flow from Selected Vertices")
		self.Bind(wx.EVT_MENU, self.glcanvas.getHeatFlow, menuHeatFlow)
		
		
		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
		menuBar.Append(operationsMenu,"&Operations")
		menuBar.Append(laplacianMenu, "&MeshLaplacian")
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		
		self.rightPanel = wx.BoxSizer(wx.VERTICAL)
		
		#Buttons to go to a default view
		viewPanel = wx.BoxSizer(wx.HORIZONTAL)
		topViewButton = wx.Button(self, -1, "Top")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromTop, topViewButton)
		viewPanel.Add(topViewButton, 0, wx.EXPAND)
		sideViewButton = wx.Button(self, -1, "Side")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromSide, sideViewButton)
		viewPanel.Add(sideViewButton, 0, wx.EXPAND)
		frontViewButton = wx.Button(self, -1, "Front")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromFront, frontViewButton)
		viewPanel.Add(frontViewButton, 0, wx.EXPAND)
		self.rightPanel.Add(wx.StaticText(self, label="Views"), 0, wx.EXPAND)
		self.rightPanel.Add(viewPanel, 0, wx.EXPAND)
		
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
		self.rightPanel.Add(self.displayMeshVerticesCheckbox, 0, wx.EXPAND)
		self.displayPrincipalAxesCheckbox = wx.CheckBox(self, label = "Display Principal Axes")
		self.displayPrincipalAxesCheckbox.SetValue(False)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayPrincipalAxesCheckbox, self.displayPrincipalAxesCheckbox)
		self.rightPanel.Add(self.displayPrincipalAxesCheckbox, 0, wx.EXPAND)
		self.useLightingCheckbox = wx.CheckBox(self, label = "Use Lighting")
		self.useLightingCheckbox.SetValue(True)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.useLightingCheckbox, self.useLightingCheckbox)
		self.rightPanel.Add(self.useLightingCheckbox, 0, wx.EXPAND)

		#Checkboxes and buttons for manipulating the cut plane
		self.rightPanel.Add(wx.StaticText(self, label="Cutting Plane"), 0, wx.EXPAND)
		self.displayCutPlaneCheckbox = wx.CheckBox(self, label = "Display Cut Plane")
		self.displayCutPlaneCheckbox.SetValue(False)
		self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayCutPlaneCheckbox, self.displayCutPlaneCheckbox)
		self.rightPanel.Add(self.displayCutPlaneCheckbox, 0, wx.EXPAND)
		CutWithPlaneButton = wx.Button(self, -1, "Cut With Plane")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.CutWithPlane, CutWithPlaneButton)
		self.rightPanel.Add(CutWithPlaneButton)
		FlipAcrossPlaneButton = wx.Button(self, -1, "Flip Across Plane")
		self.Bind(wx.EVT_BUTTON, self.glcanvas.FlipAcrossPlane, FlipAcrossPlaneButton)
		self.rightPanel.Add(FlipAcrossPlaneButton)

		#Finally add the two main panels to the sizer		
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		#cubecanvas = CubeCanvas(self)
		#self.sizer.Add(cubecanvas, 2, wx.EXPAND)
		self.sizer.Add(self.glcanvas, 2, wx.EXPAND)
		self.sizer.Add(self.rightPanel, 0, wx.EXPAND)
		
		self.SetSizer(self.sizer)
		self.Layout()
		#self.SetAutoLayout(1)
		#self.sizer.Fit(self)
		self.Show()
	
	def OnLoadMesh(self, evt):
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "OBJ files (*.obj)|*.obj|OFF files (*.off)|*.off", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			print dirname
			self.glcanvas.mesh = LaplacianMesh()
			print "Loading mesh %s..."%filename
			self.glcanvas.mesh.loadFile(filepath)
			self.glcanvas.meshCentroid = self.glcanvas.mesh.getCentroid()
			self.glcanvas.meshPrincipalAxes = self.glcanvas.mesh.getPrincipalAxes()
			print "Finished loading mesh\n %s"%self.glcanvas.mesh
			print self.glcanvas.mesh
			self.glcanvas.bbox = self.glcanvas.mesh.getBBox()
			print "Mesh BBox: %s\n"%self.glcanvas.bbox
			self.glcanvas.camera.centerOnBBox(self.glcanvas.bbox, theta = -math.pi/2, phi = math.pi/2)
			self.glcanvas.Refresh()
		dlg.Destroy()
		return

	def OnSaveMesh(self, evt):
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			self.glcanvas.mesh.saveFile(filepath, True)
			self.glcanvas.Refresh()
		dlg.Destroy()
		return

	def OnSaveMeshMeters(self, evt):
		for V in self.glcanvas.mesh.vertices:
			V.pos = 0.001*V.pos
		dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.SAVE)
		if dlg.ShowModal() == wx.ID_OK:
			filename = dlg.GetFilename()
			dirname = dlg.GetDirectory()
			filepath = os.path.join(dirname, filename)
			self.glcanvas.mesh.saveFile(filepath, True)
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
