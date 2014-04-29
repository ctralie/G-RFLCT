from PolyMesh import *
from Primitives3D import *
from Shapes3D import *
from OpenGL.GL import *
from Graphics3D import *
import sys
import re
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr, cg
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

#Note: This class assumes mesh is triangular
#(use e.g. PolyMesh.starTriangulate() if this is not the case)

def getCotangent(v1, v2, v3):
	P1 = v1.pos
	P2 = v2.pos
	P3 = v3.pos
	dV1 = P1 - P3
	dV2 = P2 - P3
	cosAngle = dV1.Dot(dV2)/(dV1.Length()*dV2.Length())
	if cosAngle**2 >= 0.99:
		return 10
	sinAngle = math.sqrt(1 - cosAngle**2)
	if sinAngle == 0:
		return 10
	return cosAngle/sinAngle

class LaplacianMesh(PolyMesh):
	#Return the sparse NxN upper part of the matrix representing the
	#Laplacian constraints, in the sparse coordinate format
	def getLaplacianSparseMatrixCoords(self, overwriteRows = None, useCotangentWeights = True):
		I = []
		J = []
		V = []
		overwriteIdx = 0
		for v1 in self.vertices:
			i = v1.ID
			if overwriteRows:
				if overwriteIdx < len(overwriteRows):
					if i == overwriteRows[overwriteIdx]:
						#If this is a row that's supposed to be overwritten
						#then overwrite it with the relevant value
						(j, v) = overwriteRows[overwriteIdx]
						I.append(i)
						J.append(j)
						V.append(v)
						overwriteIdx = overwriteIdx+1
						continue
			#Otherwise, put on laplacian constraints
			totalWeight = 0.0
			for e in v1.edges:
				v2 = e.vertexAcross(v1)
				j = v2.ID
				#Use cotangent weights
				cotAlpha = 0
				cotBeta = 0
				if e.f1:
					Vs = e.f1.getVertices()
					Vs.remove(v1)
					Vs.remove(v2)
					if len(Vs) < 1:
						print "LaplacianMesh Warning: Face with only 2 vertices encountered"
					v3 = Vs[0]
					if useCotangentWeights:
						cotAlpha = getCotangent(v1, v2, v3)
				if e.f2:
					Vs = e.f2.getVertices()
					Vs.remove(v1)
					Vs.remove(v2)
					if len(Vs) < 1:
						print "LaplacianMesh Warning: Face with only 2 vertices encountered"
					v3 = Vs[0]
					if useCotangentWeights:
						cotBeta = getCotangent(v1, v2, v3)
				wij = 1
				if useCotangentWeights:
					wij = 0.5*(cotAlpha + cotBeta)
				totalWeight = totalWeight + wij
				I.append(i)
				J.append(j)
				V.append(-wij)			
			I.append(i)
			J.append(i)
			V.append(totalWeight)
		return (I, J, V)
				

	#Inputs:
	#constraints: [[(i1, w1), (i2, w2), ..., (in, wn)], ...], M constraints
	#deltaCoords: NxY numpy array, where Y is the dimension
	#g: MxY numpy array representing values of the constraints, where Y is the dimension
	#and M is the number of constraints
	#overWriteRows: [(i1, w1), (i2, w2), ..., (in, wn)].  Incides and values of rows in
	#the upper square matrix to overwrite.  It is assumed that the indices are in increasing order
	#useCotangentWeights: Use the cotangent weights (umbrella weights used if false)
	def solveFunctionWithConstraints(self, constraints, deltaCoords, g, overwriteRows = None, useCotangentWeights = True):
		#TODO: Implement overwriteRows
		(I, J, V) = self.getLaplacianSparseMatrixCoords(overwriteRows, useCotangentWeights)
		NVerts = len(self.vertices)
		NConstraints = g.shape[0]
		Y = g.shape[1]
		for i in range(NConstraints):
			constraint = constraints[i]
			for elem in constraint:
				(index, val) = elem
				I.append(i + NVerts)
				J.append(index)
				V.append(val)
		N = NVerts + NConstraints
		M = NVerts
		b = np.append(deltaCoords, g, axis = 0)
		I = np.array(I)
		J = np.array(J)
		V = np.array(V)
		#print "V.shape = %s, N = %i, M = %i"%(V.shape, N, M)
		A = sparse.coo_matrix((V, (I, J)), shape=(N,M)).tocsr()
		ret = np.zeros((NVerts, Y))
		for i in range(Y):
			print "Solving column %i..."%i
			thisColumn = lsqr(A, b[:, i])[0]
			ret[:, i] = thisColumn
		return ret
	
	#Make a call to solveFunctionWithConstraints and update positions
	#constraints are in the form [(index, Position), ...]
	#Update the vertex positions with the results
	def solveVertexPositionsWithConstraints(self, constraints):
		#TODO: Some wasted computation computing laplacian matrix twice
		[I, J, V] = self.getLaplacianSparseMatrixCoords()
		N = len(self.vertices)
		A = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()
		x = np.zeros((N, 3))
		for i in range(N):
			P = self.vertices[i].pos
			x[i] = [P.x, P.y, P.z]
		deltaCoords = A*x
		constraintsToPass = []
		g = np.zeros((len(constraints), 3))
		for i in range(len(constraints)):
			(index, P) = constraints[i]
			constraintsToPass.append([(index, 1)])
			g[i] = [P.x, P.y, P.z]
		newPos = self.solveFunctionWithConstraints(constraintsToPass, deltaCoords, g)
		for i in range(N):
			self.vertices[i].pos = Point3D(newPos[i, 0], newPos[i, 1], newPos[i, 2])

	#anchoredVertices: array of the form [(index, Position)]
	def createMembraneSurface(self, anchoredVertices):
		print "TODO"



#An application of Laplacian meshes for blending color smoothly across a mesh
#TargetMesh: The target mesh
#CX: an array of colors for the transplanted point set
#tx: Indices of the triangles in "TargetMesh" where each of the points in X falls
#ux: Barycentric coordinates of each point
def transplantColorsLaplacianUsingBarycentric(TargetMesh, CX, tx, ux):
	NX = CX.shape[0]
	#Create a new mesh that will hold the result of the color transplant
	NewMesh = LaplacianMesh()
	#Step 1: Copy the mesh
	for v in TargetMesh.vertices:
		NewMesh.addVertex(v.pos)
	for f in TargetMesh.faces:
		verts = [NewMesh.vertices[v.ID] for v in f.getVertices()]
		NewMesh.addFace(verts)
					
	#Step 2: Setup the laplacian constraints and solve the system
	print "Solving Laplacian Mesh system for colors..."
	constraints = []
	coloredVertices = []
	for i in range(tx.shape[0]):
		triIndex = int(tx[i].flatten()[0])
		vIDs = [v.ID for v in NewMesh.faces[triIndex].getVertices()]
		constraints.append([(vIDs[0], ux[i, 0]), (vIDs[1], ux[i, 1]), (vIDs[2], ux[i, 2])])
	g = CX.copy()
	deltaCoords = np.zeros((len(NewMesh.vertices), 3))
	CY = NewMesh.solveFunctionWithConstraints(constraints, deltaCoords, g)
	#Make sure colors don't go above 1
	#CY[CY > 1] = 1
	CY = CY - CY.min()
	CY = CY/CY.max()
	#Copy colors over
	for i in range(0, CY.shape[0]):
		C = CY[i, :]
		NewMesh.vertices[i].color = [C[0], C[1], C[2]] 
	print "Finished Laplacian Mesh system..."
	print CY
	return NewMesh

#An application of Laplacian meshes for blending color smoothly across a mesh
#TargetMesh: The target mesh
#CX: an array of colors for the transplanted point set
#tx: Indices of the triangles in "TargetMesh" where each of the points in X falls
#ux: Barycentric coordinates of each point
#TODO: FINISH THIS
def transplantColorsLaplacianUsingBarycentricSubdivision(TargetMesh, CX, tx, ux):
	NX = CX.shape[0]
	#Create a new mesh that will hold the result of the color transplant
	NewMesh = LaplacianMesh()
	#Step 1: Copy over the target mesh
	for v in TargetMesh.vertices:
		NewMesh.addVertex(v.pos)
	for f in TargetMesh.faces:
		verts = [NewMesh.vertices[v.ID] for v in f.getVertices()]
		NewMesh.addFace(verts)	
	
	#Step 2: Find all of the vertices that are contained within the same
	#triangle, and subdivide that triangle
	#Add the vertices of the transplanted mesh that fall inside
	#different triangles
	tris = []
	for i in range(len(TargetMesh.faces)):
		tris.append([])
	for i in range(len(tx)):
		triIndex = int(tx[i].flatten()[0])
		triVertices = [v.pos for v in TargetMesh.faces[triIndex].getVertices()]
		newPos = Point3D(0, 0, 0)
		for k in range(len(triVertices)):
			newPos = newPos + ux[i, k]*triVertices[k]
		newVertex = NewMesh.addVertex(newPos, [CX[i, 0], CX[i, 1], CX[i, 2]])
		#Store along the color with this new vertex
		tris[triIndex].append(newVertex)	
	
	#Step 3: Setup the laplacian constraints and solve the system
	print "Solving Laplacian Mesh system for colors..."
	constraints = []
	coloredVertices = []
	for i in range(tx.shape[0]):
		triIndex = int(tx[i].flatten()[0])
		vIDs = [v.ID for v in NewMesh.faces[triIndex].getVertices()]
		constraints.append([(vIDs[0], ux[i, 0]), (vIDs[1], ux[i, 1]), (vIDs[2], ux[i, 2])])
	g = CX.copy()
	deltaCoords = np.zeros((len(NewMesh.vertices), 3))
	CY = NewMesh.solveFunctionWithConstraints(constraints, deltaCoords, g)
	#Make sure colors don't go above 1
	#CY[CY > 1] = 1
	CY = CY - CY.min()
	CY = CY/CY.max()
	#Copy colors over
	for i in range(0, CY.shape[0]):
		C = CY[i, :]
		NewMesh.vertices[i].color = [C[0], C[1], C[2]] 
	print "Finished Laplacian Mesh system..."
	print CY
	return NewMesh

#An application of Laplacian meshes for blending color smoothly across a mesh
#TargetMesh: The target mesh
#CX: an array of colors for the transplanted point set
#tx: Indices of the triangles in "TargetMesh" where each of the points in X falls
#ux: Barycentric coordinates of each point
def transplantColorsLaplacianUsingDelaunay(TargetMesh, CX, tx, ux):
	NX = CX.shape[0]
	#Create a new mesh that will hold the result of the color transplant
	NewMesh = LaplacianMesh()
	#Step 1: Find all of the vertices that are contained within the same
	#triangle, and sub-triangulate that triangle
	
	#Add the original vertices of the mesh
	for v in TargetMesh.vertices:
		NewMesh.addVertex(v.pos)
	#Add the vertices of the transplanted mesh that fall inside
	#different triangles
	tris = []
	for i in range(len(TargetMesh.faces)):
		tris.append([])
	for i in range(len(tx)):
		triIndex = int(tx[i].flatten()[0])
		triVertices = [v.pos for v in TargetMesh.faces[triIndex].getVertices()]
		newPos = Point3D(0, 0, 0)
		for k in range(len(triVertices)):
			newPos = newPos + ux[i, k]*triVertices[k]
		newVertex = NewMesh.addVertex(newPos, [CX[i, 0], CX[i, 1], CX[i, 2]])
		#Store along the color with this new vertex
		tris[triIndex].append(newVertex)
	for i in range(0, len(tris)):
		#Get the vertices in the new mesh corresponding to this triangle
		vertices = [NewMesh.vertices[v.ID] for v in TargetMesh.faces[i].getVertices()]
		if len(tris[i]) > 0:
			#Sub-triangulate the triangle where points fell
			Vs = vertices + tris[i]
			Ps = [v.pos for v in Vs]
			#Transform the Ps into an orthogonal planar coordinate system
			#to set it up for the Delaunay triangulation
			Axis1 = Ps[1] - Ps[0]
			Axis2 = Ps[2] - Ps[0]
			Axis2 = Axis1.projPerp(Axis2)
			Axis1.normalize()
			Axis2.normalize()
			PsMat = np.zeros((len(Ps), 2))
			for k in range(len(Ps)):
				dP = Ps[k] - Ps[0]
				x = Axis1.proj(dP).Length()
				y = Axis2.proj(dP).Length()
				PsMat[k, 0] = x
				PsMat[k, 1] = y
			#Do the Delaunay triangulation and add the faces
			tri = Delaunay(PsMat)
			for k in range(tri.vertices.shape[0]):
				newFace = [Vs[j] for j in tri.vertices[k, :]]
				NewMesh.addFace(newFace)
		else:
			#For every face that wasn't added, be sure to put the
			#original triangles back
			NewMesh.addFace(vertices)
					
	#Step 2: Setup the laplacian constraints and solve the system
	print "Solving Laplacian Mesh system for colors..."
	constraints = []
	coloredVertices = []
	for v in NewMesh.vertices:
		if v.color:
			coloredVertices.append(v)
			constraints.append([(v.ID, 1)])
	g = np.zeros((len(coloredVertices), 3))
	for i in range(0, len(coloredVertices)):
		g[i, :] = np.array(coloredVertices[i].color)
	deltaCoords = np.zeros((len(NewMesh.vertices), 3))
	CY = NewMesh.solveFunctionWithConstraints(constraints, deltaCoords, g)
	#Make sure colors don't go above 1
	CY[CY > 1] = 1
	CY = CY/CY.max()
	#Copy colors over
	for i in range(0, CY.shape[0]):
		C = CY[i, :]
		NewMesh.vertices[i].color = [C[0], C[1], C[2]] 
	print "Finished Laplacian Mesh system..."
	print CY
	return NewMesh

if __name__ == '__main__':
	mesh1 = getBoxMesh()
	mesh1.starTriangulate()
	mesh = LaplacianMesh()
	mesh.vertices = mesh1.vertices
	mesh.edges = mesh1.edges
	mesh.faces = mesh1.faces
	constraints = [[(1, 1)], [(3, 1)]]
	g = np.array([[1], [1]])
	x = mesh.solveFunctionWithConstraints(constraints, g)
	print x
