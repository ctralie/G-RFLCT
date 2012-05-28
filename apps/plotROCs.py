import numpy 
import pylab
import pickle

if __name__ == '__main__':
	suffix = "highsnr15"
	[Pf1, Pd1] = pickle.load(open("ROCsExp1%s.txt"%suffix, "rb"))
	[Pf2, Pd2] = pickle.load(open("ROCsExp2%s.txt"%suffix, "rb"))
	[Pf3, Pd3] = pickle.load(open("ROCsExp3%s.txt"%suffix, "rb"))
	p1s = {}
	for i in range(0, len(Pf1)):
		p1s[Pf1[i]] = Pd1[i]
	Pf1 = sorted(p1s)
	Pd1 = [p1s[i] for i in Pf1]
	print "Pf1 = %s; Pd1 = %s; Pf2 = %s; Pd2 = %s;Pf3 = %s; Pd3 = %s;"%(Pf1, Pd1, Pf2, Pd2, Pf3, Pd3)
	print "axis square;";
	print "plot(Pf1, Pd1, Pf2, Pd2, Pf3, Pd3);"
	print "legend('No Uncertainty', 'Uncertainty and Modeled', 'Uncertainty But Not Modeled');"
	print "title('15Ghz ROC Curves');"
	print "xlabel('Pf');"
	print "ylabel('Pd');"
	print "axis([0, 1, 0.5, 1]);"
	ROC1 = pylab.plot(Pf1, Pd1)
	ROC2 = pylab.plot(Pf2, Pd2)
	ROC3 = pylab.plot(Pf3, Pd3)
	pylab.legend([ROC1, ROC2, ROC3], ['BoxDims Known Exactly', 'Uncertainty Accounted For', 'Uncertainty not accounted for'])
	pylab.show()
