import numpy 
import pylab
import pickle
import math
from sys import argv
import random

random.seed()
sigma = 3.0 #Standard deviation 3cm
gausspdf = [math.exp(-t*t/(2*sigma*sigma)) for t in range(-5, 6)]
PDFNorm = sum(gausspdf)
gausspdf = [i/PDFNorm for i in gausspdf]
print gausspdf
gausscdf = [0]*len(gausspdf)
gausscdf[0] = gausspdf[0]
[W, L, H, R] = [5.0, 5.0, 2.5, 0.9]
boxUncertainty = [i/100.0 for i in range(-5, 6)]
posUncertainty = [-1, 0, 1]
dBoxProbs = {}
for i in range(0, len(boxUncertainty)):
	dBoxProbs[boxUncertainty[i]] = gausspdf[i]
print dBoxProbs

for i in range(1, len(gausspdf)):
	gausscdf[i] = gausspdf[i] + gausscdf[i-1]

def sampleBoxUncertainty():
	cdfval = random.random()
	for i in range(0, len(gausscdf)):
		if cdfval <= gausscdf[i]:
			return boxUncertainty[i]

def lambdaSKE(r, s, Es, sigmaNSqr):
	arg = -Es/(2*sigmaNSqr)
	arg = arg + sum([r[i]*s[i] for i in range(0, len(s))])/sigmaNSqr
	return math.exp(arg)

def lambdaComposite(r, sigmaNSqr, groundTruth, signalEnergies):
	ret = 0.0
	for dW in boxUncertainty:
		for dL in boxUncertainty:
			for dH in boxUncertainty:
				for rx in [-1, 0, 1]:
					for ry in [-1, 0, 1]:
						prob = dBoxProbs[dW]*dBoxProbs[dL]*dBoxProbs[dH]/9.0
						keyStr = getKeyString(W+dW, L+dL, H+dH, rx, ry, R)
						s = (groundTruth[keyStr])[1]
						Es = signalEnergies[keyStr]
						ret = ret + prob*lambdaSKE(r, s, Es, sigmaNSqr)
	return ret
	#return sum([lambdaSKE(r, signals[i], energies[i], sigmaNSqr)/float(len(signals)) for i in range(0, len(signals))])
	
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
#The commended code below was used to verify that Gaussian sampling
#was working correctly
#	for i in range(0, 1000):
#		print "%g, "%sampleBoxUncertainty(),
	print "Loading data..."
	groundTruth = pickle.load(open("groundTruth.dat", "rb"))
	print "Finished Loading data"
	#First calculate all of the signal energies
	maxEnergy = 0.0
	dSquared = 16.0
	signalEnergies = {}
	for key, data in groundTruth.items():
		signal = data[1]
		Energy = sum(val**2 for val in signal)
		signalEnergies[key] = Energy
		if Energy > maxEnergy:
			maxEnergy = Energy
	sigmaNSqr = maxEnergy/dSquared
	sigmaN = math.sqrt(sigmaNSqr)
	print "sigmaN = %g\n"%sigmaN

	#Now calculate ROC curve
	print "sigmaNSqr = %g\n"%sigmaNSqr
	print "maxEnergy = %g\n"%maxEnergy
	NExperiments = 100
	#detectionThresh = -2.5*Es/nVar:(4.0*Es/nVar)/100.0:1.5*Es/nVar;
	lnLambdas = [i*50*maxEnergy/sigmaN for i in range(-500, 500)]
	lambdas = [math.exp(i) for i in lnLambdas]
	print "Min Lambda = %g, Max lambda = %g"%(lambdas[0], lambdas[-1])
	
	simH1s = NExperiments*[0]
	simH0s = NExperiments*[0]
	#def lambdaComposite(r, sigmaNSqr, groundTruth, signalEnergies):
	for i in range(0, NExperiments):
		#Randomly choose one of the signals
		rx = posUncertainty[int(math.floor(random.random()*3))%3]
		ry = posUncertainty[int(math.floor(random.random()*3))%3]
		thisW = W + sampleBoxUncertainty()
		thisL = L + sampleBoxUncertainty()
		thisH = H + sampleBoxUncertainty()
		print "Box: [%g x %g x %g], Receiver: (%g, %g, -2)"%(thisW, thisH, thisL, rx, ry)
		keyStr = getKeyString(thisW, thisL, thisH, rx, ry, R)
		s = (groundTruth[keyStr])[1]
		noise = pylab.randn(len(s), 1)
		r = [s[k] + sigmaN*noise[k][0] for k in range(0, len(s))]
		simH1s[i] = lambdaComposite(r, sigmaNSqr, groundTruth, signalEnergies)
		print "Sim Lambda: %g"%simH1s[i]
		noise = pylab.randn(len(s), 1)
		r = [sigmaN*noise[k][0] for k in range(0, len(s))]
		simH0s[i] = lambdaComposite(r, sigmaNSqr, groundTruth, signalEnergies)
		print "Noise Lambda: %g"%simH0s[i]		
		print "Finished experiment %i of %i"%(i+1, NExperiments)
	
	#Now calculate ROC curves
	Pf = [0]*len(lambdas)
	Pd = [0]*len(lambdas)
	for i in range(0, len(lambdas)):
		cutoff = lambdas[i]
		detected = 0
		false = 0
		for k in range(0, NExperiments):
			if simH1s[k] > cutoff:
				detected = detected + 1
			if simH0s[k] > cutoff:
				false = false + 1
		Pf[i] = float(false) / float(NExperiments)
		Pd[i] = float(detected) / float(NExperiments)
	pickle.dump((Pf, Pd), open('ROCsExp2.txt', 'wb'))
	pylab.plot(Pf, Pd, 'b.')
	pylab.show()
