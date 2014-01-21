import numpy as np
import sys


def pdShiftD(x,y,N,D):
	pd_index = 0
	factor = 1
	pd_index = pd_index + (x%N + y%N) % N * factor
	for i in xrange(D-1):
		x = x / N
		y = y / N
		factor = factor * N
		pd_index = pd_index + (x%N + y%N) % N * factor
    	return pd_index


def densityToCorr(density_field,N):
    	""" 
	Performs spatial autocorrelation, taking in an arbitrarily represented NxNxN ndarray and returning an (NxNxN, NxNxN) autocorrelation.
	The result is unnormalized, and periodic boundaries are applied. Symmetry should be applied afterwards as is appropriate.
    	"""

	print "\t\tComputing correlation...",
	sys.stdout.flush()
    	Corr = np.ndarray(shape=(N*N*N,N*N*N))
    	density_field.shape = (N*N*N)

	
	for x in xrange(N*N*N):
    		for dx in xrange(N*N*N):
        		Corr[x, dx] = density_field[x] * density_field[pdShiftD(x,dx,N,3)]
	print "Done."
	
    	#print "Corr size:", Corr.size
    	#print N*N*N
    	#print "X-direction from 0 in correlation tensor: ", Corr[0,0,:]
    	#print "Shape of Correlation: ", Corr.shape
    	return Corr



if __name__ == '__main__':
	N_times = 100
	N = 10 					#boxes per spatial dimension
	density_av = np.zeros(shape=(N*N*N))
	correla_av = np.zeros(shape=(N*N*N))

	printFileInfo()
	for time_index in xrange(N_times):
		Xi, Vi, system_size, type_i = loadOneAtomTimeFromFile(time_index)
		density = atomToDensity(Xi, system_size, N, type_i, erfContribution)
		corr    = densityToCorr(density, N)

		density_av = map(sum, zip(density_av, density))
		correla_av = map(sum, zip(correla_av, corr))
