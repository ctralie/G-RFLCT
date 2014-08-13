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
from scipy.sparse.linalg import lsqr, cg, eigsh
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
	sinAngle = math.sqrt(1 - cosAngle**2)
	return cosAngle/sinAngle

def LaplaceBeltramiWeightFunc(v1, v2, v3, v4):
	w = 0.0
	if v3:
		w = w + getCotangent(v1, v2, v3)
	if v4:
		w = w + getCotangent(v1, v2, v4)
	#TODO: Fix area scaling
	#return (3.0/(2.0*v1.oneRingArea))*w
	return w

def UmbrellaWeightFunc(v1, v2, v3, v4):
	#Very simple function that just returns 1 for the umbrella weights
	return 1

class LaplacianMesh(PolyMesh):
	#Return the sparse NxN upper part of the matrix representing the
	#Laplacian constraints, in the sparse coordinate format
	def getLaplacianSparseMatrixCoords(self, overwriteRows = None, weightFunction = LaplaceBeltramiWeightFunc):
		I = []
		J = []
		V = []
		overwriteIdx = 0
		for v1 in self.vertices:
			#Precompute 1-ring areas
			v1.oneRingArea = v1.getOneRingArea()
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
				v3 = None
				v4 = None
				if e.f1:
					Vs = e.f1.getVertices()
					Vs.remove(v1)
					Vs.remove(v2)
					if len(Vs) < 1:
						print "LaplacianMesh Warning: Face with only 2 vertices encountered"
					v3 = Vs[0]
				if e.f2:
					Vs = e.f2.getVertices()
					Vs.remove(v1)
					Vs.remove(v2)
					if len(Vs) < 1:
						print "LaplacianMesh Warning: Face with only 2 vertices encountered"
					v4 = Vs[0]
				wij = weightFunction(v1, v2, v3, v4)
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
	def solveFunctionWithConstraints(self, constraints, deltaCoords, g, overwriteRows = None):
		#TODO: Implement overwriteRows
		(I, J, V) = self.getLaplacianSparseMatrixCoords(overwriteRows, LaplaceBeltramiWeightFunc)
		NVerts = len(self.vertices)
		NConstraints = g.shape[0]
		Y = deltaCoords.shape[1]
		for i in range(NConstraints):
			constraint = constraints[i]
			for elem in constraint:
				(index, val) = elem
				I.append(i + NVerts)
				J.append(index)
				V.append(val)
		N = NVerts + NConstraints
		M = NVerts
		b = deltaCoords
		if g.shape[0] > 0:
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
	#constraints are a dictionary of the form {index: Position}
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
		i = 0
		for index in constraints:
			P = constraints[index]
			constraintsToPass.append([(index, 1)])
			g[i] = [P.x, P.y, P.z]
			i = i+1
		newPos = self.solveFunctionWithConstraints(constraintsToPass, deltaCoords, g)
		for i in range(N):
			self.vertices[i].pos = Point3D(newPos[i, 0], newPos[i, 1], newPos[i, 2])
		self.needsDisplayUpdate = True
		self.needsIndexDisplayUpdate = True

	#Make a "soap bubble" surface
	#anchoredVertices: dictionary of the form {index: Position}
	def createMembraneSurface(self, constraints):
		[I, J, V] = self.getLaplacianSparseMatrixCoords()
		N = len(self.vertices)
		deltaCoords = np.zeros((N, 3))
		constraintsToPass = []
		#g = np.zeros((len(constraints), 3))
		g = np.zeros((0, 0))
		i = 0
		for index in constraints:
			P = constraints[index]
			constraintsToPass.append([(index, 1)])
			deltaCoords[i] = [P.x, P.y, P.z]
			i = i+1
		newPos = self.solveFunctionWithConstraints([], deltaCoords, g, constraintsToPass)
		for i in range(N):
			self.vertices[i].pos = Point3D(newPos[i, 0], newPos[i, 1], newPos[i, 2])
		self.needsDisplayUpdate = True
		self.needsIndexDisplayUpdate = True
	
	#Return the first k eigenvectors of the Laplace-Beltrami operator
	def getHKSEigenvectors(self, k):
		N = len(self.vertices)
		[I, J, V] = self.getLaplacianSparseMatrixCoords()
		A = sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()
		#Make symmetric so that the eigenvalues are real
		A = 0.5*(A + A.transpose())
		print "Computing %i eigenvalues for a mesh with %i vertices..."%(k, N)
		(eigvalues, eigvectors) = eigsh(A, k, which='SM')
		print "Finished computing eigenvalues"
		eigvalues[0] = 0
		return (eigvalues, eigvectors)
	
	#Return the Heat Kernel Signature of the shape
	#k: The number of eigenvalues to take to approximate the HKS
	#t: The timescale of the approximation
	#eigalues/eigvectors: If they have been precomputed don't recompute them
	def getHKS(self, k, t, eigvalues = None, eigvectors = None):
		N = len(self.vertices)
		if eigvalues.shape[0] == 0 or eigvectors.shape[0] == 0:
			(eigvalues, eigvectors) = self.getHKSEigenvectors(k)
		expScale = np.exp(-eigvalues*t)
		V = np.array(np.exp(-eigvalues*t))*eigvectors
		V = V*V
		V = np.sqrt(V.sum(1))
		return (V, eigvalues, eigvectors)
	
	#Return the heat flow from the vertices at "indices" after a specified
	#amount of time
	def getHeatFlowFromPoints(self, indices, k, t, eigvalues = None, eigvectors = None):
		N = len(self.vertices)
		initConditions = np.zeros((N, 1))
		initConditions[indices] = 1 #Start heat flow off at the selected vertices
		if eigvalues.shape[0] == 0 or eigvectors.shape[0] == 0:
			(eigvalues, eigvectors) = self.getHKSEigenvectors(k)
		expScale = np.exp(-eigvalues*t)
		coeffs = expScale*( ((eigvectors.transpose()).dot(initConditions)).flatten() )
		V = coeffs*eigvectors
		return (V.sum(1), eigvalues, eigvectors)

if __name__ == '__main__':
	mesh1 = getBoxMesh()
	mesh1.starTriangulate()
	mesh = LaplacianMesh()
	mesh.vertices = mesh1.vertices
	mesh.edges = mesh1.edges
	mesh.faces = mesh1.faces
	N = len(mesh.vertices)
	print N
	(eigenvalues, eigenvectors) = mesh.getHKSEigenvectors(6)
	print eigenvalues
