#Planar Reflective Symmetry Transform based on
#Podolak, Joshua, et al. "A planar-reflective symmetry transform for 3D shapes." 
#ACM Transactions on Graphics (TOG). Vol. 25. No. 3. ACM, 2006.

import math
from Primitives3D import *
from PolyMesh import *
import numpy as np
import numpy.linalg as linalg

def normalToSpherical(N):
	(X, Y, Z) = (N[0], N[1], N[2])
	R = math.sqrt(X**2 + Y**2 + Z**2)
	theta = math.acos(Z/R)
	phi = math.atan2(Y, X)
	#Follow what Podolak2006 suggests for avoiding a singularity at the origin
	#(also folding in redundancy)
	if theta > math.pi/2:
		#Flip the normal to the top of the sphere if it points in the bottom
		theta = theta - math.pi/2
		phi = (phi + math.pi)%(2*math.pi)
		R = -R
	return (R, theta, phi)

def normalFromSpherical(theta, phi):
	[sinTheta, cosTheta] = [math.sin(theta), math.cos(theta)]
	[sinPhi, cosPhi] = [math.sin(phi), math.cos(phi)]
	return [sinTheta*cosPhi, sinTheta*sinPhi, cosTheta]

def getBinIndex(minVal, maxVal, NIndices, x):
	dx = x - minVal
	index = int(math.floor(NIndices*dx/(maxVal-minVal)))
	if index == NIndices:
		index = NIndices - 1
	return index

def getBinAverage(minVal, maxVal, NIndices, idx):
	binRange = maxVal - minVal
	x = minVal + binRange*idx/NIndices
	return x + binRange/(2*NIndices) #Put it in the center of the bin

#TODO: Sample the points evenly on the mesh
#Nx, Ny, Nz are the bin sizes
#Returns: Tuple of (plane, pairs), where pairs is an array of point pairs 
#that voted on the plane
def PRST(mesh, NR = 64, NTheta = 64, NPhi = 64, MaxIter = 100000):
	#Center the mesh so that the binning will be more uniform
	centroid = mesh.getCentroid()
	#Copy over the points and work with them
	NPoints = len(mesh.vertices)
	Points = np.zeros((NPoints, 3))
	for i in range(NPoints):
		pos = mesh.vertices[i].pos
		Points[i, 0:3] = [pos.x - centroid.x, pos.y - centroid.y, pos.z - centroid.z]
	#Setup the grid and index into it later based on the extent
	#of R, theta, and Phi
	RExtent = mesh.getBBox().getDiagLength()
	grid = np.zeros((NR, NTheta, NPhi))
	pairs = {}
	
	#Now do the Monte Carlo algorithm
	for i in range(MaxIter):
		#Randomly sample two points from the mesh
		idx1 = np.random.randint(0, NPoints-1, 1)
		idx2 = np.random.randint(0, NPoints-2, 1)
		if idx2 >= idx1:
			idx2 = idx2 + 1
		P1 = Points[idx1, :].flatten()
		P2 = Points[idx2, :].flatten()
		N = P2 - P1 #Plane Normal (un-normalized)
		d = linalg.norm(N) #distance between two points, d, as written in paper
		(NormR, NormTheta, NormPhi) = normalToSpherical(N)
		N = N/d
		P = P1 + 0.5*(P2 - P1) #Point on plane
		R = P.dot(N) #Signed distance of the plane from the origin
		
		#Now figure out bin indices of the plane
		#R ranges from [-RExtent, RExtent]
		#theta ranges from [0, pi/2]
		#phi ranges from [0, 2pi]
		RIdx = getBinIndex(-RExtent, RExtent, NR, R)
		ThetaIdx = getBinIndex(0, math.pi/2, NTheta, NormTheta)
		PhiIdx = getBinIndex(0, 2*math.pi, NPhi, NormPhi)
		#TODO: Fix weighting
		#weight = 1.0/( d*d*math.sin(getBinAverage(0, math.pi/2, NTheta, ThetaIdx)) )
		weight = 1.0/(d*d)
		#weight = 1.0
		grid[RIdx][ThetaIdx][PhiIdx] = grid[RIdx][ThetaIdx][PhiIdx] + weight
		triple = (RIdx, ThetaIdx, PhiIdx)
		if not triple in pairs:
			pairs[triple] = []
		pairs[triple].append([idx1, idx2])
	
	#Find the bin of the max plane
	maxPair = (0, 0, 0)
	maxVal = 0
	for thispair in pairs:
		(RIdx, ThetaIdx, PhiIdx) = thispair
		if grid[RIdx][ThetaIdx][PhiIdx] > maxVal:
			maxPair = thispair
			maxVal = grid[RIdx][ThetaIdx][PhiIdx]
	(RIdx, ThetaIdx, PhiIdx) = maxPair
	R = getBinAverage(-RExtent, RExtent, NR, RIdx)
	Theta = getBinAverage(0, math.pi/2, NTheta, ThetaIdx)
	Phi = getBinAverage(0, 2*math.pi, NPhi, PhiIdx)
	
	#Convert the normal coordinates into cartesian and put the plane
	#back in the frame of the mesh
	N = np.array(normalFromSpherical(Theta, Phi))
	P0 = R*N + [centroid.x, centroid.y, centroid.z]
	N = Vector3D(N[0], N[1], N[2])
	P0 = Point3D(P0[0], P0[1], P0[2])
	print "(RIdx, ThetaIdx, PhiIdx) = (%i, %i, %i)"%(RIdx, ThetaIdx, PhiIdx)
	print "(NPairs, maxVal) = (%i, %g)"%(len(pairs[(RIdx, ThetaIdx, PhiIdx)]), maxVal)

		
	maxPairs = pairs[(RIdx, ThetaIdx, PhiIdx)]
	for i in range(len(maxPairs)):
		for j in [0, 1]:
			maxPairs[i][j] = mesh.vertices[maxPairs[i][j]].pos
	return (Plane3D(P0, N), maxPairs)
