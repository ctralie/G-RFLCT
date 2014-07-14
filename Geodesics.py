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
	
	print "======================="
	#Put v1, v2, and v3 on a plane, and set v3 to be 0
	[V1, V2, V3] = [v.pos for v in [v1, v2, v3]]
	print "V1 = %s (dist = %g), V2 = %s (dist = %g), V3 = %s (dist = %g)"%(V1, v1.FMMDist, V2, v2.FMMDist, V3, v3.FMMDist)
	V1 = V1 - V3
	V2 = V2 - V3
	V3 = V3 - V3
	print "V1 = %s (dist = %g), V2 = %s (dist = %g), V3 = %s (dist = %g)"%(V1, v1.FMMDist, V2, v2.FMMDist, V3, v3.FMMDist)
	
	#Get the coordinates of V1 and V2 on the plane spanned
	#by V1, V2, and V3, with origin V3
	Axis1 = V1
	Axis2 = V1.projPerp(V2)
	Axis1.normalize()
	Axis2.normalize()
	x1 = Axis1.Dot(V1)
	y1 = Axis2.Dot(V1)
	x2 = Axis1.Dot(V2)
	y2 = Axis2.Dot(V2)
	V = np.zeros((2, 2))
	V[0, 0] = x1
	V[1, 0] = y1
	V[0, 1] = x2
	V[1, 1] = y2
	print "V = %s\n\n"%V
	Q = np.linalg.inv( V.transpose().dot(V) )
	col1 = np.ones((2, 1))
	row1 = np.ones((1, 2))
	d = np.ones((2, 1))
	d[0] = v1.FMMDist
	d[1] = v2.FMMDist
	#Quadratic formula
	a = row1.dot( Q.dot(col1) )
	a = a.flatten()[0]
	b = -2*row1.dot( Q.dot(d) )
	b = b.flatten()[0]
	c = d.transpose().dot( (Q.dot(d)) ) - 1
	c = c.flatten()[0]
	print "Q = %s\n\nd = %s\n\nb = %s  \n\na = %s  \n\nc = %s"%(Q, d, b, a, c)
	if b**2 < 4*a*c:
		return None #(No solution)
	p = -b + math.sqrt(b**2 - 4*a*c)
	p = p/(2*a)
	p = p*np.ones((2, 1))
	p = p.flatten()[0]
	print "p = %g "%p
	n = np.linalg.inv(V.transpose())*(d - p)
	d3 = v3.FMMDist
	ConsistencyCheck = Q.dot( V.transpose().dot(n) ).flatten()
	#Check for monotonicity
	if ConsistencyCheck[0] < 0 and ConsistencyCheck[1] < 0:
		d3 = min(d3, p)
	else:
		print "NOT MONOTONE"
		d3 = min(d3, v1.FMMDist + V1.Length(), v2.FMMDist + V2.Length())
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
		#Make the first point a black point at distance 0
		V[i].FMMDist = 0
		V[i].FMMType = FMM_BLACK_TYPE
		
		Q = Queue.PriorityQueue() #Set up priority queue
		#Add the neighbors of V[i]
		for VN in V[i].getVertexNeighbors():
			length = (VN.pos - V[i].pos).Length()
			VN.FMMDist = length
			VN.FMMType = FMM_RED_TYPE
			Q.put( (length , VN ))
		
		while not Q.empty():
			(d1, v1) = Q.get()
			print "\n\n\nv1: %s, dist = %g"%(v1.pos, d1)
			if v1.FMMType == FMM_BLACK_TYPE:
				continue
			if v1.FMMDist == np.Inf:
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
						print "Updating %s to %g "%(v3.pos, v3.FMMDist)
						Q.put( (v3.FMMDist, v3) )
				if e.f2:
					v3 = updateFMM(v1, v2, e.f2)
					if v3:
						print "Updating %s to %g "%(v3.pos, v3.FMMDist)
						Q.put( (v3.FMMDist, v3) )
			v1.FMMType = FMM_BLACK_TYPE
		#Now copy over the distances to the matrix
		for k in range(0, N):
			D[i, k] = V[i].FMMDist
	return D

if __name__ == '__main__':
	mesh = PolyMesh()
	v1 = mesh.addVertex(Point3D(0, 0, 0))
	v2 = mesh.addVertex(Point3D(1, 0, 0))
	v3 = mesh.addVertex(Point3D(1, 1, 0))
	v4 = mesh.addVertex(Point3D(0, 1, 0))
	v5 = mesh.addVertex(Point3D(1, 2, 0))
	v6 = mesh.addVertex(Point3D(0, 2, 0))
	mesh.addFace([v1, v2, v3])
	mesh.addFace([v1, v3, v4])
	mesh.addFace([v3, v5, v6])
	mesh.addFace([v4, v3, v6])
	
	mesh = getOctahedronMesh()
	
	getGeodesicDistancesFMM(mesh)
	mesh.saveOffFile("out.off")
