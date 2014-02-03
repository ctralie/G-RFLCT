import math
from Primitives3D import *
from PolyMesh import *
import numpy as np
import Queue

FMM_GREEN_TYPE = 0
FMM_BLACK_TYPE = -1
FMM_RED_TYPE = 1


#Return the vertex which has been changed
#TODO: Deal with obtuse angles 
def updateFMM(v1, v2, face):
	verts = face.getVertices()
	v3 = None
	for v in verts:
		if not (v == v1 or v == v2):
			v3 = v
			break
	if not v3:
		print "Warning: No v3 found"
		return None
	if v3.FMMType == FMM_BLACK_TYPE:
		return None #Can't update points that have already been reached
	
	#Put v1, v2, and v3 on a plane, and set v3 to be 0
	[V1, V2, V3] = [v.pos for v in [v1, v2, v3]]
	V1 = V1 - V3
	V2 = V2 - V3
	V3 = V3 - V3
	
	V2Perp = V1.proj(V2)
	V2Par = V2 - V2Perp
	V = np.zeros((2, 2))
	V[0, 0] = V1.Length()
	V[1, 0] = 0
	V[0, 1] = V2Par.Length()
	V[1, 1] = V2Perp.Length()
	Q = np.linalg.inv( V.transpose().dot(V) )
	col1 = np.ones((2, 1))
	row1 = np.ones((1, 2))
	d = np.ones((2, 1))
	d[0] = v1.FMMDist
	d[1] = v2.FMMDist
	#if v1.FMMDist < np.inf and v2.FMMDist < np.inf:
	#	print "d[0] = %g, d[1] = %g"%(v1.FMMDist, v2.FMMDist)
	#Quadratic formula
	a = row1.dot( Q.dot(col1) )
	b = -2*row1.dot( Q.dot(d) )
	c = d.transpose().dot( (Q.dot(d)) ) - 1
	print "Q = %s\n\nd = %s\n\nb = %s  \n\na = %s  \n\nc = %s"%(Q, d, b, a, c)
	if b**2 < 4*a*c:
		return None #(No solution)
	p = b + math.sqrt(b**2 - 4*a*c)
	p = p/(2*a)
	p = p*np.ones((2, 1))
	p = p.flatten()[0]
	print "%g "%p
	n = np.linalg.inv(V.transpose())*(d - p)
	d3 = v3.FMMDist
	ConsistencyCheck = Q.dot( V.transpose().dot(n) ).flatten()
	if ConsistencyCheck[0] < 0 and ConsistencyCheck[1] < 0:
		d3 = min(d3, p)
	else:
		x1Len = V[0, 0]
		x2Len = math.sqrt(V[0, 1]**2 + V[1, 1]**2)
		d3 = min(d3, d3 + x1Len, d3 + x2Len)
	v3.FMMDist = d3
	v3.FMMType = FMM_RED_TYPE#This vertex is red now
	return v3

#Given the poly mesh "mesh" with N vertices returns an NxN 2D numpy array
#representing the pairwise distances between vertices of the mesh.
#Performs this computation using the Fast Marching Method (FMM)
#Returns a matrix "D" which represents the pairwise geodesic distances
#between points
def getGeodesicDistancesFMM(mesh):
	V = mesh.vertices
	N = len(V)
	D = np.Inf*np.ones((N, N)) #Initialize all distances to Infinity
	for i in range(0, 1):  #Fill in the rows of D one at a time (propagate from each vertex i)
		for k in range(0, N):
			V[k].FMMType = FMM_GREEN_TYPE
			V[k].FMMDist = np.Inf
		#Make the first point a red point at distance 0
		V[i].FMMDist = 0
		V[i].FMMType = FMM_RED_TYPE
		
		Q = Queue.PriorityQueue() #Set up priority queue
		#Add V[i] and the neighbors of V[i]
		Q.put( (V[i].FMMDist, V[i]) )
		for VN in V[i].getVertexNeighbors():
			length = (VN.pos - V[i].pos).Length()
			VN.FMMDist = length
			Q.put( (length , VN ))
		
		while not Q.empty():
			(d1, v1) = Q.get()
			if not (v1.FMMType == FMM_RED_TYPE):
				continue
			vNeighbors = v1.getVertexNeighbors()
			for v2 in vNeighbors:
				#TODO: Make faster by looking at attached faces instead of vertices?
				e = getEdgeInCommon(v1, v2)
				if not e:
					break
				if e.f1:
					v3 = updateFMM(v1, v2, e.f1)
					#Do lazy priority queue adding of new distance at vertex
					if v3:
						#print "%g "%v3.FMMDist,
						Q.put( (v3.FMMDist, v3) )
				if e.f2:
					v3 = updateFMM(v1, v2, e.f2)
					if v3:
						#print "%g "%v3.FMMDist,
						Q.put( (v3.FMMDist, v3) )
			v1.FMMType = FMM_BLACK_TYPE
		#Now copy over the distances to the matrix
		for k in range(0, N):
			D[i, k] = V[i].FMMDist
	return D
