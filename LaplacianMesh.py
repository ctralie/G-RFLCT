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
from scipy.sparse.linalg import lsqr

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

class LaplacianMesh(PolyMesh):
	#Return the sparse NxN upper part of the matrix representing the
	#Laplacian constraints, in the sparse coordinate format
	def getLaplacianSparseMatrixCoords(self):
		I = []
		J = []
		V = []
		for v1 in self.vertices:
			i = v1.ID
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
						print "Warning: Face with only 2 vertices encountered"
					v3 = Vs[0]
					cotAlpha = getCotangent(v1, v2, v3)
				if e.f2:
					Vs = e.f2.getVertices()
					Vs.remove(v1)
					Vs.remove(v2)
					if len(Vs) < 1:
						print "Warning: Face with only 2 vertices encountered"
					v3 = Vs[0]
					cotBeta = getCotangent(v1, v2, v3)
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
	def solveFunctionWithConstraints(self, constraints, deltaCoords, g):
		(I, J, V) = self.getLaplacianSparseMatrixCoords()
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
			constraintsToPass.append([index, 1])
			g[i] = [P.x, P.y, P.z]
		newPos = self.solveFunctionWithConstraints(constraintsToPass, deltaCoords, g)
		for i in range(N):
			self.vertices[i].pos = Point3D(newPos[i, 0], newPos[i, 1], newPos[i, 2])

	#anchoredVertices: array of the form [(index, Position)]
	def createMembraneSurface(self, anchoredVertices):
		print "TODO"
			

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
