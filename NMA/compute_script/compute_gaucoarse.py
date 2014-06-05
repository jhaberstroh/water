import cPickle
from atom2density import *

N_times = 500
N = 10 					#boxes per spatial dimension
density_av = np.zeros(shape=(N*N*N))
correla_av = np.zeros(shape=(N*N*N))
	
printFileInfo()
for time_index in xrange(N_times):
	Xi, Vi, system_size, type_i = loadOneAtomTimeFromFile(time_index)
	density = atomToDensity(Xi, system_size, N, type_i, gaussianContribution)
	corr    = densityToCorr(density, N)

	density_av = density_av + density
	correla_av = correla_av + corr

correla_av = correla_av / N_times
density_av = density_av / N_times

cPickle.dump( (correla_av, density_av), open('gauss_corr.pkl','wb'))
print "Saved gauss_corr.pkl."
