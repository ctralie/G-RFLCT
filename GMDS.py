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
import matplotlib.cm as cm
#from sklearn import manifold
from ICP import *
from Utilities2D import *

GMDS_MAXITER = 100

#Denoted as Dy(ti, tj) in my notes
def getDyij(MeshY, DY, ti, tj):
	itri = [v.ID for v in MeshY.faces[ti].getVertices()]
	jtri = [v.ID for v in MeshY.faces[tj].getVertices()]
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
	ti = int(ts[i].flatten()[0])
	ui = us[i, :]
	for j in js:
		if i == j:
			continue
		tj = int(ts[j].flatten()[0])
		uj = us[j, :]
		uj = uj.flatten() #Make uj a column vector
		Dyij = getDyij(MeshY, DY, ti, tj)
		b = b - (2*DX[i, j])*uj.transpose().dot(Dyij)
		APart = Dyij.transpose().dot(uj)
		A = A + np.outer(APart, APart)
	return (A, b)

#Find the approximate nearest neighbors on the mesh by searching for the closest
#points and looking for the closest point in the attached faces to that point
def getInitialGuessClosestPoints(VX, MeshY):
	NX = VX.shape[0]
	NY = len(MeshY.vertices)
	#Copy over the mesh points to a numpy array
	VY = np.zeros((NY, 3))
	for i in range(NY):
		thisP = MeshY.vertices[i].pos
		VY[i, :] = np.array([thisP.x, thisP.y, thisP.z])
	MKDTree = spatial.KDTree(VY)
	dists, idx = MKDTree.query(VX)
	ts = np.zeros((NX, 1))
	us = np.zeros((NX, 3))
	for i in range(len(idx)):
		thisPoint = Point3D(VX[i, 0], VX[i, 1], VX[i, 2])
		faces = MeshY.vertices[idx[i]].getAttachedFaces()
		minDist = np.inf
		closestPoint = None
		closestFace = None
		for f in faces:
			thisClosestP = f.getClosestPoint(thisPoint)
			distSqr = (thisClosestP - thisPoint).squaredMag()
			if distSqr < minDist:
				minDist = distSqr
				closestPoint = thisClosestP
				closestFace = f
		ts[i] = closestFace.ID
		#Find barycentric coordinates of the closest point in the face
		Vs = [v.pos for v in closestFace.getVertices()]
		Axis1 = Vs[1] - Vs[0]
		Axis2 = Vs[2] - Vs[0]
		Axis2 = Axis1.projPerp(Axis2)
		Axis1.normalize()
		Axis2.normalize()
		A = Point3D(0, 0, 0)
		B = Point3D(Axis1.Dot(Vs[1] - Vs[0]), Axis2.Dot(Vs[1] - Vs[0]), 0)
		C = Point3D(Axis1.Dot(Vs[2] - Vs[0]), Axis2.Dot(Vs[2] - Vs[0]), 0)
		X = Point3D(Axis1.Dot(closestPoint - Vs[0]), Axis2.Dot(closestPoint - Vs[0]), 0)
		us[i, :] = getBarycentricCoords(A, B, C, X, False) #Due to numerical precision skip the validity check;
		if us[i, :].sum() > 3:
			for V in Vs:
				print V
			print us[i, :]
		#points close on the edge are certainly within the triangle
	return ts, us	

def getInitialGuess2DProjection(VX, MeshY):
	NX = VX.shape[0]
	NY = len(MeshY.vertices)
	#Copy over the mesh points to a numpy array
	VY = np.zeros((NY, 3))
	for i in range(NY):
		thisP = MeshY.vertices[i].pos
		VY[i, :] = np.array([thisP.x, thisP.y, thisP.z])
	ts = 0*np.ones((NX, 1))
	us = np.zeros((NX, 3))
	us[:, 0] = 1
	(MAxis1, MAxis2, MAxis3, maxProj, minProj, MAxes) = MeshY.getPrincipalAxes()
	MCentroid = MeshY.getCentroid()
	#Transform the mesh and the the coordinate system of the principal axes of the mesh
	#where the Z coordinate can be ignored (since it's the axis of least variation)
	Trans = np.eye(4)
	Trans[0, 3] = -MCentroid.x
	Trans[1, 3] = -MCentroid.y
	Trans[2, 3] = -MCentroid.z
	Rot = np.eye(4)
	Rot[0:3, 0:3] = MAxes.transpose()
	T = Rot.dot(Trans)
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
	print "There are %i/%i valid projected triangle indices"%((ts != -1).sum(), len(ts))
	#For the points that did not fall within a triangle, find the closest points on the mesh
	idxOutside = (ts == -1).flatten()
	ts[idxOutside], us[idxOutside, :] = getInitialGuessClosestPoints(VX[idxOutside, :], MeshY)
	return ts, us

#Determine whether a vector incident on an an edge points
#inside of the corresponding triangle
def vectorInsideTriangleEdge(vec, edge, face):
	[v1, v2] = [edge.v1, edge.v2]
	v3 = None
	for v in face.getVertices():
		if not (v == v1 or v == v2):
			v3 = v
			break
	[v1, v2, v3] = [v1.pos, v2.pos, v3.pos]
	dV = v3 - v1
	dEdge = v2 - v1
	dV = dEdge.projPerp(dV) #This now represents the edge normal
	#pointing inside of the triangle
	if vec.Dot(dV) > 0:
		return True
	return False

def vectorInsideTriangleVertex(vec, v1, face):
	[v2, v3] = [None, None]
	verts = face.getVertices()
	starti = 0
	for i in range(len(verts)):
		if verts[i] == v1:
			starti = i
			break
	v2 = verts[(starti+1)%len(verts)]
	v3 = verts[(starti+2)%len(verts)]
	[v1, v2, v3] = [v1.pos, v2.pos, v3.pos]
	vec1 = v2 - v1
	vec2 = v3 - v1
	#Check to make sure vec is in between vec1 and vec2
	n = vec1 % vec2
	#Project vec onto the plane spanned by vec1 and vec2
	vec = n.projPerp(vec)
	#Get 2D coordinates for each of vec, vec1, and vec2
	Axis1 = vec1
	Axis2 = Axis1.projPerp(vec2)
	Axis1.normalize()
	Axis2.normalize()
	u = Axis1.Dot(vec)
	v = Axis2.Dot(vec)
	u2 = Axis1.Dot(vec2)
	v2 = Axis2.Dot(vec2)
	det1 = np.zeros((3, 3))
	det1[0, :] = np.array([1, 0, 0])
	det1[1, :] = np.array([1, 1, 0])
	det1[2, :] = np.array([u, v, 0])
	det1 = linalg.det(det1)
	det2 = np.zeros((3, 3))
	det2[0, :] = np.array([1, 0, 0])
	det2[1, :] = np.array([u, v, 0])
	det2[2, :] = np.array([u2, v2, 0])
	det2 = linalg.det(det2)
	#Check to see if the sign of the determinants is the same
	if np.sign(det1) == np.sign(det2):
		return True
	return False

#For debugging, plot a heatmap of the stress function evaluated
#at different points in the triangle
def plotRasterizedStress(A, b, triBefore, triNew, ubefore, unew):
	#Spread some barycentric coordinates out
	#and get colors proportional to the stress at
	#those coordinates
	N = 100
	Extent = 1.2
	ux, uy = np.mgrid[-Extent:Extent:np.complex(0, N), -Extent:Extent:np.complex(0, N)]
	uz = 1 - (ux + uy)
	u = np.zeros((3, ux.flatten().shape[0]))
	u[0, :] = ux.flatten()
	u[1, :] = uy.flatten()
	u[2, :] = uz.flatten()
	Stresses = A.dot(u)
	Stresses = (u*Stresses).sum(0)
	Stresses = Stresses + b.dot(u)
	Stresses = Stresses.reshape(N, N)
	#Stresses = np.log(Stresses - Stresses.min() + 0.1)
	plt.imshow(Stresses, extent = [-Extent, Extent, -Extent, Extent], origin = 'lower')
	plt.title('Stresses in triangle')
	plt.hold(True)
	plt.plot(unew[1], unew[0], 'rx')
	plt.plot([0, 1, 0, 0], [1, 0, 0, 1], 'r')
	plt.colorbar()
	plt.show()

def plotRasterizedStressEdges(A, b, tri, ubefore, unew):
	N = 100
	Extent = 1.2
	ux, uy = np.mgrid[-Extent:Extent:np.complex(0, N), -Extent:Extent:np.complex(0, N)]
	uz = 1 - (ux + uy)
	u = np.zeros((3, ux.flatten().shape[0]))
	u[0, :] = ux.flatten()
	u[1, :] = uy.flatten()
	u[2, :] = uz.flatten()
	Stresses = A.dot(u)
	Stresses = (u*Stresses).sum(0)
	Stresses = Stresses + b.dot(u)
	Stresses = Stresses.reshape(N, N)
	#Stresses = np.log(Stresses - Stresses.min() + 0.1)
	plt.subplot(2, 2, 0)
	plt.imshow(Stresses, extent = [-Extent, Extent, -Extent, Extent], origin = 'lower')
	plt.title('Stresses in triangle')
	plt.hold(True)
	plt.plot(unew[1], unew[0], 'rx')
	plt.plot(ubefore[1], ubefore[0], 'gx')
	plt.plot([1, 0], [0, 1], 'r')
	plt.plot([0, 1], [0, 0], 'g')
	plt.plot([0, 0], [1, 0], 'b')
	Stress = unew.transpose().dot(A.dot(unew)) + b.dot(unew)
	StressBefore = ubefore.transpose().dot(A.dot(ubefore)) + b.dot(ubefore)
	StressDiff = Stress - StressBefore
	plt.title("Stress: %g, StressDiff: %g"%(Stress, StressDiff))
	print unew
	plt.colorbar()
	#The edge (0, 1, 0) - (1, 0, 0): Red
	#The edge (0, 0, 1) - (0, 1, 0): Green
	#The edge (0, 0, 1) - (1, 0, 0): Blue
	ua = np.arange(0, 1, 0.01)
	ub = 1 - ua
	colors = ['r', 'g', 'b']
	for i in range(3):
		u = np.zeros((3, len(ua)))
		u[i, :] = ua
		u[(i+1)%3, :] = ub
		Stresses = A.dot(u)
		Stresses = (u*Stresses).sum(0)
		Stresses = Stresses + b.dot(u)
		Stresses = Stresses.flatten()
		plt.subplot(2, 2, i+1)
		plt.plot(ua, Stresses, colors[i])
		plt.title('Min: %g'%np.min(Stresses))
	plt.show()	

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
		PX.points.append(Point3D(VX[i, 0], VX[i, 1], -VX[i, 2]))
		PX.colors.append([1, 0, 0])

	####Step 2: Rigidly align X and Y using ICP
	print "Doing ICP..."
	AllTransformations, minIndex, error = ICP_PointsToMesh(PX, MeshY, False, False, True)
	VXTrans = np.zeros(VX.shape)
	for i in range(0, len(PX.points)):
		VXTrans[i, :] = np.array([PX.points[i].x, PX.points[i].y, PX.points[i].z])
	print "Finished ICP"
	
	####Step 3: Project onto the axis of least variation and use barycentric coordinates
	#in triangles as the intialization for GMDS
	print "Getting initial guess..."
	ts, us = getInitialGuessClosestPoints(VXTrans, MeshY)
	#ts, us = getInitialGuess2DProjection(VXTrans, MeshY)
	print "Finished Getting Initial guess"
	
	####Step 4: Run the GMDS algorithm to refine t and u
	thisiter = 0
	#Precompute all of the stress matrices, and incrementally update
	#them as they change (more efficient than recomputing each time)
	As = [None]*NX
	bs = [None]*NX
	print "Computing initial gradients..."
	for i in range(NX):
		(As[i], bs[i]) = getAandb(DX, DY, MeshY, ts, us, i, range(0, NX))
		print ".",
		if i%30 == 0:
			print ""
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
		print "thisiter = %i, maxIndex = %i, maxStressGradMag = %g"%(thisiter, maxIndex, maxStressGradMag)
		
		#Solve the constrained optimization problem to update
		#ui and ti with all j neq i vertices fixed
		S = np.zeros((4, 4))
		S[0:3, 0:3] = 2*As[maxIndex]
		S[0:3, 3] = -1
		S[3, 0:3] = 1
		blambda = np.ones(4)
		blambda[0:3] = -bs[maxIndex].flatten()
		ulambda = linalg.solve(S, blambda)
		unew = ulambda[0:3]
		thisStress = unew.dot(As[maxIndex].dot(unew)) + bs[maxIndex].dot(unew)
		print "thisStress = %g"%thisStress
		tnew = int(ts[maxIndex].flatten()[0])
		triFace = MeshY.faces[tnew]
		
		
		if (unew >= 0).sum() < 3:
			print "CONSTRAINTS ARE ACTIVE"
			#The constrained global optimum is outside of the triangle
			#Search over all edges an vertices to find the valid solution 
			#that contributes to the minimum amount of stress
			minStress = np.inf
			validEdgeFound = False
			minEdge = [0, 0]
			for i in range(0, 3):
				edge = np.array([i, (i+1)%3, 3])
				Sedge = S[edge, :][:, edge]
				blambdaedge = blambda[edge]
				ulambda = linalg.solve(Sedge, blambdaedge)
				thisu = np.zeros(4)
				thisu[edge] = ulambda
				thisu = thisu[0:3]
				if (thisu >= 0).sum() == 3:
					#A valid max has been found within the edge
					validEdgeFound = True
					thisStress = thisu.dot(As[maxIndex].dot(thisu)) + bs[maxIndex].dot(thisu)
					if thisStress < minStress:
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
				if nextFace: #If this isn't at a boundary
					tnext = nextFace.ID
					#Translate the barycentric coordinates to the new face
					#where they are possibly in a different order
					nextFaceVertices = nextFace.getVertices()
					unext = np.zeros(3)
					for i in range(0, 3):
						if nextFaceVertices[i] == v1:
							unext[i] = bary1
						elif nextFaceVertices[i] == v2:
							unext[i] = bary2
					#Check to make sure the negative of the new gradient 
					#would point inside of the triangle
					usnew = us.copy()
					tsnew = ts.copy()
					us[maxIndex, :] = unext
					tsnew[maxIndex, :] = tnext
					Anew, bnew = getAandb(DX, DY, MeshY, tsnew, usnew, maxIndex, range(0, NX))
					grad = (2*Anew.dot(unext) + bnew).flatten()
					if vectorInsideTriangleEdge(Vector3D(-grad[0], -grad[1], -grad[2]), edge, nextFace):
						unew = unext
						tnew = tnext
						print "Moving across triangle edge"
				#Otherwise, keep the point on the edge in the current triangle
				
			#Also search over all vertices
			minVertexIdx = -1
			for i in range(0, 3):
				thisu = np.zeros(3)
				thisu[i] = 1
				thisStress = thisu.dot(As[maxIndex].dot(thisu)) + bs[maxIndex].dot(thisu)
				if thisStress < minStress:
					minStress = thisStress
					minVertex = i
					unew = thisu
			if minVertexIdx >= 0:
				#If a vertex has been found with stress less than the min stress edge
				minVertex = triFace.getVertices()[minVertexIdx]
				print "minVertex = %i"%minVertexIdx
				#Search for a valid face to translate to, where the negative gradient
				#points inside of the face
				for f in minVertex.getAttachedFaces():
					if f.ID == tnew:
						continue #Don't check the current face
					newi = 0
					verts = f.getVertices()
					#Figure out what the barycentric coordinate would be
					#in the new face
					for i in range(0, 3):
						if verts[i] == minVertex:
							newi = i
					usnew = us.copy()
					tsnew = ts.copy()
					usnew[maxIndex, :] = 0
					usnew[maxIndex, newi] = 1
					tsnew[maxIndex] = f.ID
					Anew, bnew = getAandb(DX, DY, MeshY, tsnew, usnew, maxIndex, range(0, NX))
					grad = (2*Anew.dot(unew) + bnew).flatten()
					if vectorInsideTriangleVertex(Vector3D(-grad[0], -grad[1], -grad[2]), minVertex, f):
						print "Moving across vertex"
						unew = usnew[maxIndex, :]
						tnew = f.ID
						break
		
		#Update the stress matrices of the other coordinates to reflect
		#the change in this variable
		usNew = us.copy()
		usNew[maxIndex, :] = unew
		tsNew = ts.copy()
		tsNew[maxIndex] = tnew
		#plotRasterizedStress(As[maxIndex], bs[maxIndex], triFace, MeshY.faces[tnew], us[maxIndex, :].flatten(), unew)
		for i in range(NX):
			if i == maxIndex:
				continue
			#Only update the parts of the sum corresponding to the maxIndex
			AOld, bOld = getAandb(DX, DY, MeshY, ts, us, i, [maxIndex])
			ANew, bNew = getAandb(DX, DY, MeshY, tsNew, usNew, i, [maxIndex])
			As[i] = As[i] - AOld + ANew
			bs[i] = bs[i] - bOld + bNew
		#plotRasterizedStressEdges(As[maxIndex], bs[maxIndex], MeshY.faces[tnew], us[maxIndex, :], unew)
		print "tbefore = %i"%ts[maxIndex].flatten()
		print "tnew = %i"%tnew
		print "ubefore = %s"%us[maxIndex, :].flatten()
		print "unew = %s\n\n"%unew
		
		ts[maxIndex] = tnew
		us[maxIndex, :] = unew
		thisiter = thisiter + 1
	
	#Return the coordinates
	return (ts, us)
