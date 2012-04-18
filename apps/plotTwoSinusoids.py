import numpy 
import pylab
import pickle
from sys import argv

def plotSinusoids(times, data1, data2, str1, str2):
	print str1
	print str2
	times = [1e9*times[i] for i in range(0, len(times))]
	p1 = pylab.plot(times, data1)
	p2 = pylab.plot(times, data2)
	pylab.xlabel('Time (nanoseconds)')
	pylab.ylabel('Received Signal')
	pylab.title('Comparison of Received Sinusoids')
	pylab.legend([p1, p2], [str1, str2])
	pylab.show()

#W, L, H are dimensions of box
#rx and ry are the receiver x and y coordinates
#R is the reflection coefficient of the box
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
	if len(argv) < 13:
		print "Usage: plotImpulseResponse W1 L1 H1 rx1 ry1 R1 W2 L2 rx2 ry2 R2"
	else:
		print "Loading data..."
		groundTruth = pickle.load(open("groundTruthAllBefore.dat", "rb"))
		print "Finished loading data..."
		[W1, L1, H1, rx1, ry1, R1] = [float(argv[i]) for i in range(1, 7)]
		str1 = getKeyString(W1, L1, H1, rx1, ry1, R1)
		data1 = groundTruth[str1]
		[W2, L2, H2, rx2, ry2, R2] = [float(argv[i]) for i in range(7, 13)]
		str2 = getKeyString(W2, L2, H2, rx2, ry2, R2)
		data2 = groundTruth[str2]
		plotSinusoids(data1[0], data1[1], data2[1], "Box:[%gx%gx%g] Pos:(%g, %g, -2), R:%g"%(W1, L1, H1, rx1, ry1, R1), "Box:[%gx%gx%g] Pos:(%g, %g, -2), R:%g"%(W2, L2, H2, rx2, ry2, R2))
