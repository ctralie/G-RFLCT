import numpy 
import pylab
import pickle
from sys import argv

def plotImpulseResponse(response, W, L, H, rx, ry, R):
	imptimes = [1e9*(response[i])[1] for i in range(0, len(response))]
	impamps = [(response[i])[0] for i in range(0, len(response))]
	pylab.stem(imptimes, impamps)
	pylab.xlabel('Time (nanoseconds)')
	pylab.ylabel('Attenuation')
	pylab.title('Impulse Response for Box [%gx%gx%g] Source at (%g, %g, -2) Reflection Coefficient %g'%(W, L, H, rx, ry, R))
	pylab.show()

#W, L, H are dimensions of box
#rx and ry are the receiver x and y coordinates
#R is the reflection coefficient of the box
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
	if len(argv) < 7:
		print "Usage: plotImpulseResponse W L H rx ry R"
	else:
		[W, L, H, rx, ry, R] = [float(argv[i]) for i in range(1, 7)]
		print "Loading data..."
		groundTruth = pickle.load(open("groundTruth915.dat", "rb"))
		print "Finished loading data..."
		data = groundTruth[getKeyString(W, L, H, rx, ry, R)]
		plotImpulseResponse(data[2], W, L, H, rx, ry, R)
