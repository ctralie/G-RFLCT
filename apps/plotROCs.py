import numpy 
import pylab
import pickle

if __name__ == '__main__':
	[Pf1, Pd1] = pickle.load(open("ROCsExp1.txt", "rb"))
	[Pf2, Pd2] = pickle.load(open("ROCsExp2.txt", "rb"))
	print "Pf1 = %s; Pd1 = %s; Pf2 = %s; Pd2 = %s;"%(Pf1, Pd1, Pf2, Pd2)
	ROC1 = pylab.plot(Pf1, Pd1)
	ROC2 = pylab.plot(Pf2, Pd2)
	pylab.legend([ROC1, ROC2], ['BoxDims Known Exactly', 'Uncertainty in BoxDims'])
	pylab.show()
