# Computing the 5-d integrals and solving the matrix equations.
import numpy as np
import numpy.polynomial.legendre as LG
import numpy.linalg as la
import matplotlib.pylab as plt

def leg_l(l_index, numgrid):
	leg_access = np.zeros(l_index + 1)
	leg_access[l_index] = 1
	#print leg_access
	x = np.linspace(0,1,numgrid)
	y = LG.legval(x, leg_access)
	return x,y



def rn(n_index, R_sph, numgrid):
	r = np.linspace(0,R_sph,numgrid)
	func = np.power((r - R_sph)/ R_sph, n_index) 
	#print r
	return r, func

C1G = 30  
C2G = 30
R1G = 20
R2G = 20
PG =  10

RAD_SPHERE = 3

def Jnlmj_R(n,l,m,j, R, chi_x, chi_y):
	return 2 * Jh11nlmj_R(n,l,m,j, R, chi_x, chi_y) \
		+ 2 * Jh12nlmj_R(n,l,m,j, R, chi_x, chi_y) \
		+ Pnl_Pmj(n,l,m,j)

def Pnl_Pmj(n,l,m,j):
	c1, l1 = leg_l(l,C1G)
	c2, l2 = leg_l(j,C2G)
	r1, rn1 = rn(n, RAD_SPHERE, R1G)
	r2, rn2 = rn(m, RAD_SPHERE, R2G)
	f1 = np.outer(l1,rn1).flatten()
	f2 = np.outer(l2,rn2).flatten()

	volume = (c1[1] - c1[0]) \
		* (c2[1] - c2[0]) \
		* (r1[1] - r1[0]) \
		* (r1[1] - r1[0]) 

	return np.sum(np.outer(f1,f2)) * (2 * np.pi)**2 * volume


def Jh11nlmj_R(n,l,m,j, R, chi_x, chi_y):
	return Jh12nlmj_R(n,l,m,j,0,chi_x,chi_y)

def Jh12nlmj_R(n,l,m,j, R, chi_x, chi_y):
	c1, l1 = leg_l(l,C1G)
	c2, l2 = leg_l(j,C2G)
	s1 = np.sqrt( 1 - np.square(c1) )
	s2 = np.sqrt( 1 - np.square(c2) )
	r1, rn1 = rn(n, RAD_SPHERE, R1G)
	r2, rn2 = rn(m, RAD_SPHERE, R2G)
	phi = np.linspace(0, 2*np.pi, PG, endpoint=False)

	volume = (phi[1] - phi[0])  \
		* (c1[1] - c1[0]) \
		* (c2[1] - c2[0]) \
		* (r1[1] - r1[0]) \
		* (r1[1] - r1[0]) 
		
	#print "Volume: ", volume
	#print "dPhi: ", phi[1] - phi[0]
	
	cosphi = np.cos(phi)

	#r1r2 = np.outer(r1,r2)
	#r1r2.shape = (R1G,  1,R2G,  1,  1)
	#r1c1 = np.outer(r1,c1)
	#r1c1.shape = (R1G,C1G,  1,  1,  1)
	#r2c2 = np.outer(r2,c2)
	#r2c2.shape = (  1,  1,R2G,C2G,  1)

	r1.shape     = (R1G,  1,  1,  1,  1)
	r2.shape     = (  1,  1,R2G,  1,  1)
	c1.shape     = (  1,C1G,  1,  1,  1)
	c2.shape     = (  1,  1,  1,C2G,  1)
	s1.shape     = (  1,C1G,  1,  1,  1)
	s2.shape     = (  1,  1,  1,C2G,  1)
	cosphi.shape = (  1,  1,  1,  1, PG)

	good_zero = np.zeros((R1G,C1G,R2G,C2G, PG))
	r1 = good_zero + r1
	r2 = good_zero + r2
	c1 = good_zero + c1
	c2 = good_zero + c2
	s1 = good_zero + s1
	s2 = good_zero + s2

	print good_zero.shape
	cosphi = good_zero + cosphi

	dist = np.square(r1) + np.square(r2) + R**2 \
			+ 2 * r1 * r2 * c1 * c2 \
			- 2 * r1 * c1 * R \
			- 2 * r2 * c2 * R \
			- 2 * r1 * r2 * s1 * s2

	chi_dx = chi_x[1] - chi_x[0]

	dist_ind = np.floor(dist/chi_dx).astype(int)
	max_el = len(chi_y) - 1
	print "repairing indexes...",
	bool_ind = dist_ind > max_el
	#print np.histogram(bool_ind)
	dist_ind[bool_ind] = max_el
	
	#disthist,edges = np.histogram(dist)
	print "computing distance responses...",
	chi_dist = chi_y[dist_ind]
	
	print "done!",
	vec1 = np.outer(rn1, l1)
	#lv.shape = (100,)
	vec2 = np.outer(rn2, l2)
	#print "rnv shape:",rnv.shape
	#rnv.shape = (100,)

	
	res = np.einsum('ij,kl',vec1, vec2)
	print res.shape
	res.shape = (R1G,C1G,R2G,C2G,  1)
	res = good_zero + res

	res = res * chi_dist
	
	#res = np.outer(np.outer(t1,t2), np.outer(r1,r2))
	return np.sum(res) * volume

def bnl(n,l):
	c1, l1 = leg_l(l,C1G)
	r1, rn1 = rn(n, RAD_SPHERE, R1G)
	volume = c1[1] - c1[0] *\
		r1[1] - r1[0]
	return np.sum(np.outer(l1,rn1)) * volume

def loadData(fname="Sk.dat"):
	with open(fname) as f:
		x = []
		for l in f:
			x.append(np.array([float(i) for i in str.split(l)]))
	x = np.array(x)
	#print x.shape
	#print x[:,0]
	#print x[:,1]
	return x[:,0], x[:,1]

def loadChi(fname="Sk.dat"):
	k, Sk = loadData(fname)

	dx = np.pi / (k[-1] - k[0])

	x = np.linspace(0,np.pi/(k[1]-k[0]),len(k))
	print x[0]
	print dx
	rhobar = .033

	hww = np.fft.irfft(-1.0j * Sk * k)[0:201] / (4 * np.pi * np.pi) * np.reciprocal(x)
	hww /= rhobar
	hww /= rhobar
	#hww += 1
	hww[0] = hww[1]
	##hww /= hww[0] * -1
	##hww += 1
	##hww *=  200
	##hww[0] += 1/dx
	##plt.plot(Sk)
	#hww = np.ones((201,))
	#plt.plot(x, hww[0:201])
	#plt.xlabel("Distance, angstroms")
	#plt.ylabel("hww or something")
	#plt.show()

	return x, hww[0:201]



if __name__ == "__main__":
	print "Running text of matrix element solver"
	
	x,hww = loadChi()
	#integral, dist = Jnlmj_R(1,1,1,1,8, x, hww)

	n_n = 5
	n_l = 6
	N,L = np.meshgrid(range(n_n),range(n_l))

	nl_index = np.array(zip(N.flatten(), L.flatten()))


	Jnlmj_mtx = np.zeros((n_n*n_l, n_n*n_l))
	bmj_vec = np.zeros((n_n*n_l,))

	print "COMPUTING MATRIX ELEMENTS"
	for i,nl in enumerate(nl_index):
		bmj_vec[i] = bnl(nl[0],nl[1])
		print "b_mj: ", bmj_vec[i]
		for j,mj in enumerate(nl_index[i:]):
			print "nl mj:",nl[0],nl[1],mj[0],mj[1]
			ij_val = Jnlmj_R(nl[0],nl[1],mj[0],mj[1],8, x, hww)
			Jnlmj_mtx[i,j] = ij_val
			Jnlmj_mtx[j,i] = ij_val
			print "Jnlmj:", Jnlmj_mtx[i,j]

	print "J Matrix:"
	print Jnlmj_matrix
	print "B Vector:"
	print bmj_vec
	
	x = la.solve(Jnlmj_mtx, bmj_vec)
	print x

	#print "RESULT: ", integral

	#plt.hist(dist.flatten(), normed=True, bins=np.linspace(0,199,200))
	#plt.show()
	#plt.plot(np.zeros(201))
	#plt.ylim([-.2,1.5])
	#plt.show()

	
