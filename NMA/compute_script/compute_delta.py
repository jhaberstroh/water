import cPickle
from atom2density import *

N_times = 500
N = 15 					#boxes per spatial dimension
density_av = np.zeros(shape=(N*N*N))
correla_av = np.zeros(shape=(N*N*N))
	
printFileInfo()
for time_index in xrange(N_times):
	Xi, Vi, system_size, type_i = loadOneAtomTimeFromFile(time_index)
	density = atomToDeltaDensity(Xi, system_size, N, type_i)
	corr    = densityToCorr(density, N)

	density_av = density_av + density
	correla_av = correla_av + corr

correla_av = correla_av / N_times
density_av = density_av / N_times

fname = 'delta_corr15.pkl'
cPickle.dump( (correla_av, density_av), open(fname,'wb'))
print "Saved "+fname+"!"
