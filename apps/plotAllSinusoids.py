import numpy 
import pylab
import pickle
from sys import argv

def plotSinusoids(times, data):
	times = [1e9*times[i] for i in range(0, len(times))]
	plots = []
	legend = ["-1 -1", "-1 0", "-1 1", "0 -1", "0 0", "0 1", "1 -1", "1 0", "1 1"]
	for dataPoint in data:
		plots.append(pylab.plot(times, dataPoint))
	pylab.xlabel('Time (nanoseconds)')
	pylab.ylabel('Received Signal')
	pylab.title('Comparison of Received Sinusoids')
	pylab.legend(plots, legend)
	pylab.show()

#W, L, H are dimensions of box
#rx and ry are the receiver x and y coordinates
#R is the reflection coefficient of the box
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
	print "Loading data..."
	groundTruth = pickle.load(open("groundTruth.dat", "rb"))
	print "Finished loading data..."
	data = []
	times = []
	for rx in [-1, 0, 1]:
		for ry in [-1, 0, 1]:
			keyStr = getKeyString(5, 5, 2.5, rx, ry, 0.9)
			thisdata = groundTruth[keyStr]
			times = thisdata[0]
			data.append(thisdata[1])
	plotSinusoids(times, data)
