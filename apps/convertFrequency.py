import pickle
import math

def getSteadyStateSinusoid(fMhz, impulses, SampsPerCycle = 40, NCycles = 10):
	T = 1/(fMhz*1e6)
	dt = T/SampsPerCycle
	signal = [0]*SampsPerCycle*NCycles
	times = [dt*i for i in range(0, len(signal))]
	for impulse in impulses:
		A = impulse[0]
		phi = impulse[1]
		pathLen = impulse[2]
		#Attenuation also needs to be updated
		A = A*pow(915.0/fMhz, 2)
		signal = [signal[i] + A*math.cos(2*math.pi*fMhz*1e6*(dt*i - phi)) for i in range(0, len(signal))]
	return (times, signal)

if __name__ == '__main__':
	print "Loading data..."
	groundTruth = pickle.load(open("groundTruth915.dat", "rb"))
	print "Finished Loading data"
	fMhz = 15*1000 #15 Ghz
	item = 0
	for key, data in groundTruth.items():
		impulses = data[2]
		(times, signal) = getSteadyStateSinusoid(fMhz, impulses)
		groundTruth[key] = (times, signal, data[2])
		print "Finished conversion %i of %i"%(item, len(groundTruth))
		item = item + 1
	pickle.dump(groundTruth, open('groundTruth.dat', 'wb'))
