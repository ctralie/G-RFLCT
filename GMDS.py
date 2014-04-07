from PolyMesh import *
from Primitives3D import *
from Utilities2D import *
from Shapes3D import *
from OpenGL.GL import *
from Graphics3D import *
import sys
import re
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr
import scipy.spatial as spatial
#from sklearn import manifold
from ICP import *
from Utilities2D import *

GMDS_MAXITER = 1000

#Denoted as Dy(ti, tj) in my notes
def getDyij(MeshY, DY, ti, tj):
	itri = [v.id for v in MeshY.faces[ti].getVertices()]
	jtri = [v.id for v in MeshY.faces[tj].getVertices()]
	Dyij = np.zeros((3, 3))
	for j in range(0, 3):
		for i in range(0, 3):
			Dyij[j, i] = DY[itri[i], jtri[j]]
	return Dyij	

#Get the stress due to vertex i
def getStress(DX, DY, MeshY, ts, us, i):
	stress = 0.0
	ti = ts[i]
	ui = us[i, :]
	for i in range(DX.shape[0]):
		if i == j:
			continue
		tj = ts[j]
		uj = us[j, :]
		Dyij = getDYij(MeshY, DY, ti, tj)
		stress = stress + (DX[i, j] - uj.dot(Dyij.dot(ui)))**2
	return stress

#Get the gradient matrices corresponding to the vertex i
def getAandb(DX, DY, MeshY, ts, us, i):
	A = np.zeros((3, 3))
	b = np.zeros((1, 3))
	ti = ts[i]
	ui = us[i, :]
	for i in range(DX.shape[0]):
		if i == j:
			continue
		tj = ts[j]
		uj = us[j, :]
		Dyij = getDYij(MeshY, DY, ti, tj)
		b = b - (2*DX[i, j])*uj.dot(Dyij)
		APart = (uj.transpose()).dot(Dyij)
		A = A + np.outer(APart, APart)
	return (A, b)

#VX: NX x 3 numpy array of vertex positions in the X point set
#DX: NX x NX matrix of geodesic distances in the point set
#MeshY: PolyMesh object representing the target mesh
#DY: NY x NY matrix of geodesic distances along the mesh
#NOTE: This method assumes MeshY is a triangular mesh
def GMDSPointsToMesh(VX, DX, MeshY, DY):
	####Step 1: Set up the mesh and point cloud variables
	NX = VX.shape[0] #Number of points in the point set
	NY = DY.shape[0] #Number of points in the mesh
	MY = len(MeshY.faces) #Number of triangles in the mesh
	#Copy over points to a point cloud object to be used with
	#the ICP class
	PX = PointCloud()
	for i in range(NX):
		PX.points.append(Point3D(VX[i, 0], VX[i, 1], VX[i, 2]))
		PX.colors.append([1, 0, 0])
	#Copy over the mesh points to a numpy array
	VY = np.zeros((NY, 0))
	for i in range(NY):
		thisP = mesh.vertices[i].pos
		VY[i, :] = np.array([thisP.x, thisP.y, thisP.z])

	####Step 2: Rigidly align X and Y using ICP
	AllTransformations, minIndex, error = ICP_PointsToMesh(PX, MeshY, False, False, True)
	
	####Step 3: Project onto the axis of least variation and use barycentric coordinates
	#in triangles as the intialization for GMDS
	ts = np.zeros((NX, 1))
	us = np.zeros((NX, 3))
	us[:, 0] = 1
	(MAxis1, MAxis2, MAxis3, maxProj, minProj, MAxes) = M.getPrincipalAxes()
	#Transform the mesh and the the coordinate system of the principal axes of the mesh
	#where the Z coordinate can be ignored (since it's the axis of least variation)
	T = np.eye(4)
	T[0:3, 0:3] = MAxes.transpose()
	VX2D = transformPoints(T, VX)
	VX2D = VX2D[:, 0:2]
	VY2D = transformPoints(T, VY)
	VY2D = VY2D[:, 0:2]
	Y2DKDTree = spatial.KDTree(VY2D)
	dists, idx = Y2DKDTree.query(VX2D)
	for i in range(len(idx)):
		thisVertex = mesh.vertices[idx[i]]
		thisP = Point3D(VX2D[i, 0], VX2D[i, 1], VX2D[i, 2])
		for f in thisVertex.getAttachedFaces():
			YPoints = [Point3D(VY2D[v.ID, 0], VY2D[v.ID, 1], VY2D[v.ID, 2] for v in f.getVertices()]
			if pointInsideConvexPolygon2D(YPoints, thisP):
				#Use the first triangle that contains this point as the initial
				#guess, and use its barycentric coordinates in 2D
				[A, B, C] = YPoints
				ts[i] = f.ID
				us[i] = getBarycentricCoords(A, B, C, thisP)
	
	####Step 4: Run the GMDS algorithm to refine t and u
	thisiter = 0
	while thisiter < GMDS_MAXITER:
		#Find the vertex which gives rise to the maximum gradient
		#when everything else is fixed
		maxIndex = 0
		maxStressGradMag = 0
		AMax = np.eye(3)
		bMax = np.zeros((1, 3))
		for i in range(NX):
			(A, b) = getAandb(DX, DY, MeshY, ts, us, i)
			stressGradMag = A.dot(us[i, :]) + b.dot(us[i, :])
			if stressGradMag > maxStressGradMag:
				maxStressGradMag = stressGradMag
				maxIndex = i
				AMax = A
				bMax = b
		#Solve the constrained optimization problem to update
		#ui and ti with all j neq i vertices fixed
		
		thisiter = thisiter + 1
	
	#Return the coordinates
	return (ts, us)
