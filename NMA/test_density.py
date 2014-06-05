import numpy as np
import GromCorr.gro2atom as g2a
import GromCorr.atom2density as a2d


def getEvenDist():
	den_sum = np.zeros((1000))
	N = 10
	for i in xrange(1000):
		den = a2d.atomCoarsen([i%N+ .5, i/N%N + .5, i/N/N + .5], 10, 10, 1.5, a2d.erfContribution)
		den_sum = den_sum + den
	return den_sum

	

def getDiffDelta(time_slice):
	Xi, Vi, size, type_i = g2a.loadOneAtomTimeFromFile(time_slice);
	
	d1 = a2d.atomToDeltaDensity(Xi, size, 10, type_i)
	d2 = a2d.atomToDeltaDensityV2(Xi, size, 10, type_i)

	return d1 - d2;
