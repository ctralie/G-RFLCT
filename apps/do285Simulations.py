from Primitives3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from EMScene import *
import numpy 
import pylab
import pickle

def plotImpulseResponse(response):
	imptimes = [1e9*(response[i])[1] for i in range(0, len(response))]
	impamps = [(response[i])[0] for i in range(0, len(response))]
	(times, sinusoid) = scene.getSteadyStateSinusoid(fMhz, 40, 10)
	pylab.stem(imptimes, impamps)
	pylab.xlabel('Time (nanoseconds)')
	pylab.ylabel('Attenuation')
	pylab.show()

#W, L, H are dimensions of box
#rx and ry are the receiver x and y coordinates
#R is the reflection coefficient of the box
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

def decipherKeyString(arg):
	return arg.split("_")

def getGroundTruthForBox(W, L, H, R, groundTruth):
	fMhz = 915
	scene = EMScene()
	boxMesh = getBoxMesh(W, L, H, Point3D(0, 0, 0))
	EMMat = EMMaterial(R, 0)
	sceneNode = EMNode(scene.rootEMNode, boxMesh, Matrix4(), EMMat, OpticalMaterial())
	scene.rootEMNode.children.append(sceneNode)
	scene.getMeshList()
	#The roles of the source and receiver are actually interchanged here
	#since the whole thing is symmetric and precomputing virtual receiver images
	#is faster because the receiver stays fixed
	scene.Source = Point3D(0.3, -1.2, 2.0) #This is the receiver (which doesn't move)
	#Precomupte receiver images
	scene.buildVirtualSourceTree(4)
	for rx in [-1.0, 0.0, 1.0]:
		for ry in [-1.0, 0.0, 1.0]:
			keyStr = getKeyString(W, L, H, rx, ry, R)
			scene.Receiver = Point3D(rx, ry, -2.0)
			scene.getPathsToReceiver()
			impulseResponse = scene.getImpulseResponse(fMhz)
			sinusoidResponse = scene.getSteadyStateSinusoid(fMhz, 40, 10)
			#Store a tuple with times, signal, impulse response
			groundTruth[keyStr] = (sinusoidResponse[0], sinusoidResponse[1], impulseResponse)
			print "Finished simulation %i"%len(groundTruth)

if __name__ == '__main__':
	#Use a dictionary to store the ground truth signals
	groundTruth = {}
	boxUncertainty = [i/100.0 for i in range(-5, 6)]
	boxUncertainty = [0]
	[W, L, H] = [5.0, 5.0, 2.5]
	R = 0.9
	setting = 1
	for dW in boxUncertainty:
		for dL in boxUncertainty:
			for dH in boxUncertainty:
				getGroundTruthForBox(W+dW, L+dL, H+dH, R, groundTruth)
				print "Finished composite setting %i of %i (dW = %g, dL = %g, dH = %g)"%(setting, pow(len(boxUncertainty), 3), dW, dL, dH)
				setting = setting+1
	pickle.dump(groundTruth, open("groundTruth.dat", "wb"))
	print "Finished uncertainty around box size (saving results)"
	print "Now doing reflection uncertainty"
	RUncertainty = [i/20.0 for i in range(1, 21)]
	for R in RUncertainty:
		getGroundTruthForBox(W, L, H, R, groundTruth)
	print "Finished uncertainty around reflection"
	pickle.dump(groundTruth, open("groundTruth.dat", "wb"))
