from Primitives3D import *
from Graphics3D import *
from PolyMesh import *
from Cameras3D import *
from EMScene import *
import numpy 
import pylab

def plotImpulseResponse(response):
	imptimes = [1e9*(response[i])[1] for i in range(0, len(response))]
	impamps = [(response[i])[0] for i in range(0, len(response))]
	(times, sinusoid) = scene.getSteadyStateSinusoid(fMhz, 40, 10)
	pylab.stem(imptimes, impamps)
	pylab.xlabel('Time (nanoseconds)')
	pylab.ylabel('Attenuation')
	pylab.show()

if __name__ == '__main__':
	fMhz = 915
	scene = EMScene()
	boxMesh = getBoxMesh(5.1, 5.1, 2.5, Point3D(0, 0, 0))
	EMMat = EMMaterial(0.9, 0)
	sceneNode = EMNode(scene.rootEMNode, boxMesh, Matrix4(), EMMat, OpticalMaterial())
	scene.rootEMNode.children.append(sceneNode)
	scene.Source = Point3D(0, 0, -2.0)
	scene.Receiver = Point3D(0, -1.0, 2.0)
	scene.getMeshList()
	scene.buildVirtualSourceTree(4)
	scene.getPathsToReceiver()
	response = scene.getImpulseResponse(fMhz)
	plotImpulseResponse(response)
