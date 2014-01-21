import math
from operator import mul
import itertools as itt
import numpy as np
import os
import time
import sys



def erfContribution(dx, sigma, N):
    	"""
    	Solves for the distance to the closest periodic replica of x1 to x2. Integrates a gaussian of stdev sigma.
    	x1 and x2 are in reduced units, and N is the number of reduced units in the cubic system.
	Quantities in dx must be in [-N, N].
	Uses the same distance computation as GaussianContribution.
    	"""
	half = float(N)/2.0
    	dist_vec = [ abs(abs(abs(x)-half)-half) for x in dx]
	half_box = .5/sigma
    	my_erf = lambda x : - math.erf(-half_box - 2*half_box*x) + math.erf(half_box - 2*half_box*x)
    	contribution = reduce( mul, [ my_erf(x) for x in dist_vec] ) / 8
    	return contribution

def gaussianContribution(dx, sigma, N):
    	"""
    	Solves for the distance to the closest periodic replica of x1 to x2. Computes a gaussian of stdev sigma.
    	x1 and x2 are in reduced units, and N is the number of reduced units in the cubic system.
	Quantities in dx must be in [-N, N].
	Uses the same distance computation as ErfContribution.
    	"""
	half = float(N)/2.0
    	dist_vec = [ abs(abs(abs(x)-half)-half) for x in dx]
    	dist_scaled = sum([x*x for x in dist_vec]) / (2 * sigma * sigma)
    	return math.exp(-dist_scaled)

def distance(dx, sigma, N):
	return 

def atomCoarsen(x, system_len, N, sigma, method):
	"""
	Computes the coarse-grained density for the atom using method.
	Method should be either erfContribution or gaussianContribution.
	The benefit of this method:
		Easy checking of method for smoothness & sensibility
		Extension to velocity binning

	/input 	x 		natural-unit 3-tuple of position x,y,z, 
					Domain [0,system_len]
		system_len 	a scalar of unit natural unit describing the box size
		N 		the number of boxes that comprises one edge of the cubic density field
		sigma 		natural-unit the "width" of the smearing of the particle
		method 		a function with args (dx, sigma, N) that returns the contribution of a particle to a bin, given:
			dx		reduced-unit distance
						Domain [-N, N]
			sigma		reduced-unit sigma
			N		reduced-unit number of boxes, as described above

	/return density_field, a np.ndarray((N*N*N))

	/throws None
	"""
	scale = float(N) / system_len
	# Shift by .5 box length in all directions to compute the distance from the center of the box... ?
	density_field = [method( ((x[0]- box_pos%N)   *scale - .5,
				  (x[1]-(box_pos/N)%N)*scale - .5,
				  (x[2]- box_pos/N /N)*scale - .5), sigma/scale, N) for box_pos in xrange(N*N*N)]
	return np.array(density_field)


def atomToDeltaDensity(time_slice, system_len, N, type_i):
    	density_field = np.zeros((N,N,N))
    	weight_dict = {'OW': 15.9994, 'HW1': 1.008, 'HW2':1.008} 	#amu
	avg_weight = ( 15.9994 + 1.008 * 2 ) / 3
    	weight_i = [weight_dict[t] for t in type_i]
    	
      	histbins = ( np.linspace(0, system_len, N+1), ) * 3		# Tuple ctor, in nm
    	this_Xti_clean = time_slice
    	this_Xti_clean = (this_Xti_clean + system_len) %system_len  	# In nm
  
	print "\tComputing density field...",
	sys.stdout.flush()
    	density_field,edges = np.histogramdd(this_Xti_clean, bins = histbins, weights = weight_i)

	print "Binned mass: ", np.sum(density_field)
	print "System mass: ", len(this_Xti_clean) * avg_weight
	assert abs(np.sum(density_field) - len(this_Xti_clean) * avg_weight) < .1
	
    	V = system_len * system_len * system_len / N / N / N 		# Volume of one cell in nm^3
    	density_field = np.array( density_field / V * 1.66054 ) 	# conversion for amu/nm^3 -> kg/m^3, as well as to np.array
        	
	print "Done."
	density_field.shape = (N*N*N)
    	return density_field						#kg/m^3



def atomToDensity(time_slice, system_len, N, type_i, smear_fn):
    	density_field = np.zeros((N*N*N))
    	weight_dict = {'OW': 15.9994, 'HW1': 1.008, 'HW2':1.008} 					#amu
    	size_dict = {'OW': .3 * N /system_len, 'HW1': .2 * N /system_len, 'HW2': .2 * N /system_len}   	#nm -> boxes
    	weight_i = [weight_dict[t] for t in type_i]
    	size_i = [size_dict[t] for t in type_i]
    	size_sq_i = [size_dict[t] * size_dict[t] for t in type_i]
    	
    	this_Xti_clean = time_slice
    	this_Xti_clean = ((this_Xti_clean + system_len) %system_len) /system_len *N  # In reduced units
    	
	# VERY slow iterator, should probably change this out...
    	den_iter = np.nditer(density_field, flags=['multi_index'], op_flags=['writeonly'])
    	count = 0
    	clock_0 = time.clock()
    	last = 0
	print "\tComputing density field...",
	sys.stdout.flush()

	for den_iter in xrange(N*N*N):
        	count = count + 1
        	if int(time.clock() - clock_0) / 3 > last:
        		print 't =', time.clock() - clock_0, "s :", "{0:.2f} %".format(float(count) / float(density_field.size) * 100)
            		last = int(time.clock() - clock_0) / 3

        	neg_pos1 = [den_iter % N, (den_iter / N) % N, den_iter / N / N]
		#negpos1 = [ x, y, z]
        	
        	# Map the particle iterator & the weights to a single iterator via imap and lambda.
        	smear_lambda = lambda particle, size, weight: smear_fn( map(sum, zip(neg_pos1, particle)),  size, N) * weight
        	density_field[den_iter] = sum( itt.imap( smear_lambda, this_Xti_clean, size_i, weight_i ) )

	print "Done."
    	return density_field						# ?? / reduced_units



def atomToDensityV2(time_slice, system_len, N, type_i, smear_fn):
    	density_field = np.zeros((N*N*N))
	# Clean the positions
    	this_Xti_clean = time_slice
    	this_Xti_clean = ((this_Xti_clean + system_len) %system_len) /system_len *N  # In reduced units

	# Prepare the molecular details
    	weight_dict = {'OW': 15.9994, 'HW1': 1.008, 'HW2':1.008} 					#amu
    	size_dict = {'OW': .3 * N /system_len, 'HW1': .2 * N /system_len, 'HW2': .2 * N /system_len}   	#nm -> boxes
    	weight_i = [weight_dict[t] for t in type_i]
    	size_i = [size_dict[t] for t in type_i]
    	size_sq_i = [size_dict[t] * size_dict[t] for t in type_i]
    	
	# Preparing timing functions
    	count = 0
    	clock_0 = time.clock()
    	last = 0
	print "\tComputing density field...",
	sys.stdout.flush()

	for i, pos in enumerate(this_Xti_clean):
        	count = count + 1
        	if int(time.clock() - clock_0) / 3 > last:
        		print 't =', time.clock() - clock_0, "s :", "{0:.2f} %".format(float(count) / float(this_Xti_clean.size) * 3. * 100)
            		last = int(time.clock() - clock_0) / 3

		# Computing the acutal data
		atom_coarse = atomCoarsen(pos, system_len, N, size_i[i], smear_fn);
        	density_field = density_field + atom_coarse;

	print "Done."
    	return density_field						# ?? / reduced_units



def atomToDeltaDensityV2(time_slice, system_len, N, type_i):
    	density_field = np.zeros((N*N*N))
	# Clean the positions
    	this_Xti_clean = time_slice
    	this_Xti_clean = ((this_Xti_clean + system_len) %system_len) /system_len *N  # In reduced units

	# Prepare the molecular details
    	weight_dict = {'OW': 15.9994, 'HW1': 1.008, 'HW2':1.008} 					#amu
    	size_dict = {'OW': .0001 * N /system_len, 'HW1': .0001 * N /system_len, 'HW2': .0001 * N /system_len}   	#nm -> boxes
    	weight_i = [weight_dict[t] for t in type_i]
    	size_i = [size_dict[t] for t in type_i]
    	size_sq_i = [size_dict[t] * size_dict[t] for t in type_i]
    	
	# Preparing timing functions
    	count = 0
    	clock_0 = time.clock()
    	last = 0
	print "\tComputing density field...",
	sys.stdout.flush()

	for i, pos in enumerate(this_Xti_clean):
        	count = count + 1
        	if int(time.clock() - clock_0) / 3 > last:
        		print 't =', time.clock() - clock_0, "s :", "{0:.2f} %".format(float(count) / float(this_Xti_clean.size) * 3. * 100)
            		last = int(time.clock() - clock_0) / 3

		# Computing the acutal data
		atom_coarse = atomCoarsen(pos, system_len, N, size_i[i], erfContribution) * weight_i[i];
        	density_field = density_field + atom_coarse;

	print "Done."
    	return density_field						# ?? / reduced_units



def atomToDensityErf(time_slice, system_len, N, type_i):
    	density_field = np.zeros((N*N*N))
    	weight_dict = {'OW': 15.9994, 'HW1': 1.008, 'HW2':1.008} 					#amu
    	size_dict = {'OW': .3 * N /system_len, 'HW1': .2 * N /system_len, 'HW2': .2 * N /system_len}   	#nm -> boxes
    	weight_i = [weight_dict[t] for t in type_i]
    	size_i = [size_dict[t] for t in type_i]
    	size_sq_i = [size_dict[t] * size_dict[t] for t in type_i]
    	
    	this_Xti_clean = time_slice
    	this_Xti_clean = ((this_Xti_clean + system_len) %system_len) /system_len *N  # In reduced units
	print this_Xti_clean.shape
    	
	# VERY slow iterator, should probably change this out...
    	den_iter = np.nditer(density_field, flags=['multi_index'], op_flags=['writeonly'])
    	count = 0
    	clock_0 = time.clock()
    	last = 0
	print "\tComputing density field...",
	sys.stdout.flush()

	for den_iter in xrange(N*N*N):
        	count = count + 1
        	if int(time.clock() - clock_0) / 3 > last:
        		print 't =', time.clock() - clock_0, "s :", "{0:.2f} %".format(float(count) / float(density_field.size) * 100)
            		last = int(time.clock() - clock_0) / 3

        	neg_pos_box = [den_iter % N, (den_iter / N) % N, den_iter / N / N]
		half_box = float(N)/2.
		#negpos1 = [ x, y, z]
        	
        	# Map the particle iterator & the weights to a single iterator via imap and lambda.
        	#smear_lambda = lambda particle, size, weight: smear_fn( map(sum, zip(neg_pos1, particle)),  size, N) * weight
        	#density_field[den_iter] = sum( itt.imap( smear_lambda, this_Xti_clean, size_i, weight_i ) )

		for (particle_iter, xyz) in enumerate(this_Xti_clean):
			for (j, dist) in enumerate(xyz):
				dx = neg_pos_box[j] - dist
				dx = abs(dx - int(dx/half_box)*half_box) / size_i[particle_iter]

				density_field[den_iter] = density_field[den_iter] + ((math.erf(-half_box - 2 * half_box * dx) + math.erf(half_box - 2*half_box*dx)) * weight_i[particle_iter] / 2)
	print "Done."
    	return density_field						# ?? / reduced_units




