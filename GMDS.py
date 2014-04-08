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
import numpy.linalg as linalg
import scipy.spatial as spatial
import matplotlib.pyplot as plt
#from sklearn import manifold
from ICP import *
from Utilities2D import *

GMDS_MAXITER = 1000

#Denoted as Dy(ti, tj) in my notes
def getDyij(MeshY, DY, ti, tj):
	itri = [v.id for v in MeshY.faces[ti].getVertices()]
	jtri = [v.id for v in MeshY.faces[tj].getVertices()]
	Dyij = np.zeros((3, 3))
	for j in range(3):
		for i in range(3):
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
#due to all of the fixed vertices in "js"
def getAandb(DX, DY, MeshY, ts, us, i, js):
	A = np.zeros((3, 3))
	b = np.zeros((1, 3))
	ti = ts[i]
	ui = us[i, :]
	for j in js:
		if i == j:
			continue
		tj = ts[j]
		uj = us[j, :]
		Dyij = getDYij(MeshY, DY, ti, tj)
		b = b - (2*DX[i, j])*uj.dot(Dyij)
		APart = (uj.transpose()).dot(Dyij)
		A = A + np.outer(APart, APart)
	return (A, b)

def getInitialGuess2DProjection(VX, MeshY):
	NX = VX.shape[0]
	NY = len(MeshY.vertices)
	#Copy over the mesh points to a numpy array
	VY = np.zeros((NY, 3))
	for i in range(NY):
		thisP = MeshY.vertices[i].pos
		VY[i, :] = np.array([thisP.x, thisP.y, thisP.z])	
	ts = np.zeros((NX, 1))
	us = np.zeros((NX, 3))
	us[:, 0] = 1
	(MAxis1, MAxis2, MAxis3, maxProj, minProj, MAxes) = MeshY.getPrincipalAxes()
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
		thisVertex = MeshY.vertices[idx[i]]
		thisP = Point3D(VX2D[i, 0], VX2D[i, 1], 0)
		for f in thisVertex.getAttachedFaces():
			YPoints = [Point3D(VY2D[v.ID, 0], VY2D[v.ID, 1], 0) for v in f.getVertices()]
			if pointInsideConvexPolygon2D(YPoints[0:3], thisP, 0):
				#Use the first projected triangle that contains this point
				#as the initial guess, and use its barycentric coordinates in 2D
				[A, B, C] = YPoints
				ts[i] = f.ID
				us[i] = getBarycentricCoords(A, B, C, thisP)
	print "There are %i/%i nonzero triangle indices"%((ts != 0).sum(), len(ts))
	return ts, us

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

	####Step 2: Rigidly align X and Y using ICP
	AllTransformations, minIndex, error = ICP_PointsToMesh(PX, MeshY, False, False, True)
	
	####Step 3: Project onto the axis of least variation and use barycentric coordinates
	#in triangles as the intialization for GMDS
	ts, us = getInitialGuess2DProjection(VX, MeshY)
	
	####Step 4: Run the GMDS algorithm to refine t and u
	thisiter = 0
	#Precompute all of the stress matrices, and incrementally update
	#them as they change (more efficient than recomputing each time)
	As = [None]*NX
	bs = [None]*NX
	print "Computing initial gradients..."
	for i in range(NX):
		(As[i], bs[i]) = getAandb(DX, DY, MeshY, ts, us, i, range(0, NX))
	print "Finished computing initial gradients"
	#Do the optimization iterations
	while thisiter < GMDS_MAXITER:
		#Find the vertex which gives rise to the maximum magnitude gradient
		#when everything else is fixed
		maxIndex = 0
		maxStressGradMag = 0
		for i in range(NX):
			stressGrad = 2*As[i].dot(us[i, :]) + bs[i]
			stressGradMag = linalg.norm(stressGrad)
			if stressGradMag > maxStressGradMag:
				maxStressGradMag = stressGradMag
				maxIndex = i
		#Solve the constrained optimization problem to update
		#ui and ti with all j neq i vertices fixed
		triFace = MeshY.faces[ts[maxIndex]]
		S = np.zeros((4, 4))
		S[0:3, 0:3] = As[maxIndex]
		S[0, 3] = -1
		S[3, 0:3] = 1
		blambda = np.ones((4, 1))
		blambda[0:3] = -bs[maxIndex]
		ulambda = solve(S, blambda)
		unew = ulambda[0:3]
		tnew = ts[maxIndex]
		if (unew >= 0).sum() < 3:
			#The constrained global optimum is outside of the triangle
			#First search over all edges
			minStress = np.inf
			validEdgeFound = False
			minEdge = [0, 0]
			for i in range(0, 3):
				edge = np.array([i, (i+1)%3, 3])
				Sedge = S[edge, :][:, edge]
				blambdaedge = blambda[edge]
				ulambda = solve(Sedge, blambdaedge)
				thisu = np.zeros((1, 4))
				thisu[edge] = ulambda
				thisu = thisu[0:3]
				if (thisu >= 0).sum() == 3:
					#A valid max has been found within the edge
					validEdgeFound = True
					thisStress = thisu.dot(As[maxIndex].dot(thisu)) + bs[maxIndex].dot(thisu)
					if thisStress < minStress:
						#TODO: Check to make sure it makes sense to take min stress
						minStress = thisStress
						minEdge = [i, (i+1)%3]
						unew = thisu
			if validEdgeFound:
				#Figure out which triangle to translate to (if any)
				faceVertices = triFace.getVertices()
				v1 = faceVertices[minEdge[0]]
				bary1 = unew[minEdge[0]]
				v2 = faceVertices[minEdge[1]]
				bary2 = unew[minEdge[1]]
				edge = MeshY.getEdge(v1, v2)
				nextFace = edge.faceAcross(triFace)
				tnext = nextFace.ID
				#Translate the barycentric coordinates to the new face
				#where they are possibly in a different order
				nextFaceVertices = nextFace.getVertices()
				unext = np.zeros((1, 3))
				for i in range(0, 3):
					if nextFaceVertices[i] == v1:
						unext[i] = bary1
					elif nextFaceVertices[i] == v2:
						unext[i] = bary2
				#Check to make sure the new gradient would point inside of the triangle
				#TODO: Finish this gradient consistency check
				#Anew, bnew = getAandb(DX, DY, MeshY, ts, us, i, js)
				unew = unext
				tnew = tnext
			else:
				#Search over all of the vertices if no valid edge has been found
				minStress = np.inf
				minVertex = 0
				for i in range(0, 3):
					thisu = np.zeros((0, 3))
					thisu[i] = 1
					thisStress = thisu.dot(As[maxIndex].dot(thisu)) + bs[maxIndex].dot(thisu)
					if thisStress < minStress:
						minStress = thisStress
						minVertex = i
				minFace = 0
		#Update the stress matrices of the other coordinates to reflect
		#the change in this variable
		usNew = us.copy()
		usNew[maxIndex, :] = unew
		tsNew = ts.copy()
		tsNew[maxIndex] = tnew
		for i in range(NX):
			if i == maxIndex:
				continue
			AOld, bOld = getAandb(DX, DY, MeshY, ts, us, i, [maxIndex])
			ANew, bNew = getAandb(DX, DY, MeshY, tsNew, usNew, i, [maxIndex])
			As[i] = As[i] - AOld + ANew
			bs[i] = bs[i] - bOld + bNew
		ts[maxIndex] = tnew
		us[maxIndex, :] = unew
		thisiter = thisiter + 1
	
	#Return the coordinates
	return (ts, us)
