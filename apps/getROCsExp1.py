import numpy 
import pylab
import pickle
import math
from sys import argv

def lambdaSKE(r, s, Es, sigmaNSqr):
	arg = -Es/(2*sigmaNSqr)
	arg = arg + sum([r[i]*s[i] for i in range(0, len(s))])/sigmaNSqr
	return math.exp(arg)

def lambdaComposite(signals, energies, sigmaNSqr, r):
	return sum([lambdaSKE(r, signals[i], energies[i], sigmaNSqr)/float(len(signals)) for i in range(0, len(signals))])
	
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
	[W, L, H, R] = [5.0, 5.0, 2.5, 0.9]
	print "Loading data..."
	groundTruth = pickle.load(open("groundTruth915.dat", "rb"))
	print "Finished Loading data"
	#First calculate all of the signal energies
	maxEnergy = 0.0
	dSquared = 160.0*18.7/31.97 #Average dSquared should be about 18.7 for 915 and 15
	signalEnergies = {}
	#groundTruth[keyStr] = (sinusoidResponse[0], sinusoidResponse[1], impulseResponse)
	allEnergy = 0.0
	for key, data in groundTruth.items():
		signal = data[1]
		Energy = sum(val**2 for val in signal)
		allEnergy = allEnergy + Energy
		signalEnergies[key] = Energy
		if Energy > maxEnergy:
			maxEnergy = Energy
	allEnergy = allEnergy / len(groundTruth)
	sigmaNSqr = maxEnergy/dSquared
	print "Average dSquared = %g"%(allEnergy/sigmaNSqr)
	sigmaN = math.sqrt(sigmaNSqr)
	print "sigmaN = %g\n"%sigmaN
	
	#Get all of the possible composite signals
	signals = []
	energies = []
	for rx in [-1, 0, 1]:
		for ry in [-1, 0, 1]:
			keyStr = getKeyString(W, L, H, rx, ry, R)
			s = (groundTruth[keyStr])[1]
			signals.append(s)
			energies.append(signalEnergies[keyStr])
	
	#Now calculate ROC curve
	print "sigmaNSqr = %g\n"%sigmaNSqr
	print "maxEnergy = %g\n"%maxEnergy
	maxLambda = 2*math.exp((maxEnergy**2-maxEnergy/2)/sigmaNSqr)
	print "maxLambda = %g"%maxLambda
	NROCPoints = 100
	NExperiments = 100
	#detectionThresh = -2.5*Es/nVar:(4.0*Es/nVar)/100.0:1.5*Es/nVar;
	lnLambdas = [0.5*i*maxEnergy/sigmaN for i in range(-100, 200)]
	lambdas = [math.exp(i) for i in lnLambdas]
	Pf = [0]*len(lambdas)
	Pd = [0]*len(lambdas)
	for i in range(0, len(lambdas)):
		cutoff = lambdas[i]
		detected = 0
		false = 0
		num = 0
		#First run trials for H1 and cycle through the possible signals
		for s in signals:
			for trial in range(0, NExperiments):
				noise = pylab.randn(len(s), 1)
				r = [s[k] + sigmaN*noise[k][0] for k in range(0, len(s))]
				val = lambdaComposite(signals, energies, sigmaNSqr, r)
				if val > cutoff:
					detected = detected + 1
				num = num + 1
		Pd[i] = float(detected) / float(num)
		num = 0
		#Now run trials for H0
		for trial in range(0, NExperiments*len(signals)):
			noise = pylab.randn(len(signals[0]), 1)
			r = [sigmaN*noise[k][0] for k in range(0, len(noise))]
			val = lambdaComposite(signals, energies, sigmaNSqr, r)
			if val > cutoff:
				false = false + 1
			num = num + 1
		Pf[i] = float(false) / float(num)
		print "Finished lambda = %g (%i of %i): Pf = %g, Pd = %g"%(cutoff, i+1, len(lambdas), Pf[i], Pd[i])
		pickle.dump((Pf, Pd), open('ROCsExp1highsnr915.txt', 'wb'))
	pylab.plot(Pf, Pd, 'b.')
	pylab.show()
