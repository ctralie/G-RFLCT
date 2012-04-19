import numpy 
import pylab
import pickle
import math
from sys import argv
import random

def lambdaSKE(r, s, sigmaNSqr):
	arg = -s*s/(2*sigmaNSqr)
	arg = arg + r*s/sigmaNSqr
	return math.exp(arg)
	
def getKeyString(W, L, H, rx, ry, R):
	return "%g_%g_%g_%g_%g_%g"%(W, L, H, rx, ry, R)

if __name__ == '__main__':
	[W, L, H, R] = [5.0, 5.0, 2.5, 0.9]
	print "Loading data..."
	groundTruth = pickle.load(open("groundTruth.dat", "rb"))
	print "Finished Loading data"
	#First calculate all of the signal energies
	maxEnergy = 0.0
	dSquared = 16.0
	signalEnergies = {}
	#groundTruth[keyStr] = (sinusoidResponse[0], sinusoidResponse[1], impulseResponse)
	for key, data in groundTruth.items():
		signal = data[1]
		Energy = sum(val**2 for val in signal)
		signalEnergies[key] = Energy
		if Energy > maxEnergy:
			maxEnergy = Energy
	sigmaNSqr = maxEnergy/dSquared
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
	
	for choiceSignal in range(0, 9):
		print "postx%i = ["%choiceSignal,
		s = signals[choiceSignal]
		noise = pylab.randn(len(s), 1)
		r = [s[k] + sigmaN*noise[k][0] for k in range(0, len(s))]
		prior = [1.0/9.0]*9
		posts = [] #Posterior distributions at different timesteps
		lastPost = prior
		lastLambdaConds = [1]*9 #lambdas conditioned 
		lambdas = []
		lastLambda = 1
		for k in range(0, len(r)):
			rVal = r[k]
			newLambdas = [lambdaSKE(rVal, signals[i][k], sigmaNSqr) for i in range(0, 9)]
			lambdaKCond = sum([newLambdas[i]*lastPost[i] for i in range(0, 9)])
			lastLambda = lastLambda*lambdaKCond
			lambdas.append(lastLambda)
			#def lambdaSKE(r, s, sigmaNSqr)
			lastPost = [newLambdas[i]*lastPost[i]/lambdaKCond for i in range(0, 9)]
			posts.append(lastPost)
			for i in range(0, len(lastPost)):
				print "%g"%lastPost[i],
				if i < len(lastPost) - 1:
					print ",",
				else:
					print ";",
		print "];"
