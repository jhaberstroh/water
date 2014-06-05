def DensityLine():
	import numpy as np
	import csv
	with open('Density5000CorrAvLin.csv','r') as f:
		read = csv.reader(f)

		corr = []
		for pos, csvline in enumerate(read):
			if (pos == 0):
				den = [float(x) for x in csvline]
			else:
				corr.append( np.array([float(x) for x in csvline]) )

	return np.vstack(np.array(corr)[0:1000]), np.array(den)


def OpenGr():
	import numpy as np
	import csv
	with open('GrCorrelation.csv','r') as f:
		read = csv.reader(f)

		corr = []
		for pos, csvline in enumerate(read):
			if (pos == 0):
				den  = np.array([float(x) for x in csvline])
			else
				corr = np.array([float(x) for x in csvline])

	return corr, den

def PlotOneSlice(Img3d, zslice):
	import numpy as np
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	plt.ion()
	fig = plt.figure();
	ax = fig.add_subplot(111,projection='3d')

	xv, yv = np.meshgrid(xrange(0,10), xrange(0,10))
	zv = np.array([ Img3d[x,y,zslice] for x,y in zip(xv,yv) ] )
	print xv.shape, yv.shape, zv.shape
	lines  = ax.plot_wireframe(xv,yv,zv)
	plt.plot()
	lines.remove()

def PlotSlices(Img3d, size):
	import numpy as np
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from time import sleep

	plt.ion()
	fig = plt.figure();
	ax = fig.add_subplot(111,projection='3d')

	plt.show()

	xv, yv = np.meshgrid(xrange(0,10), xrange(0,10))

	for z in xrange(size):
		zv = np.array([ Img3d[x,y,z] for x,y in zip(xv,yv) ] )
		print xv.shape, yv.shape, zv.shape
		lines  = ax.plot_wireframe(xv,yv,zv)
		plt.draw()
		sleep(1)
		lines.remove()



def MAIN1_Corr2D():
	corr, den = DensityLine()
	print "loaded"
	import numpy.linalg as LA
	w,v = LA.eig(corr)


	e0 = v[0]
	print "eig"
	e0.shape = (10,10,10)
	#print "Plotting Slices: "
	#PlotSlices(e0, 10)
	import numpy.fft as FFT
	ft0 = np.array([abs(x) for x in FFT.fftn(e0)])
	access = np.array([[[max(i,j,k) <= 5 for k in xrange(10)] for j in xrange(10)] for i in xrange(10)])
	ftpos = ft0[access]


def MAIN2_Corr1D():
	corr, den = OpenGr()
	PlotOneSlice(corr, 1)
