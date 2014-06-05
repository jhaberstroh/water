from matplotlib import pyplot
import pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import random

def main():
	with open('cover.txt','r') as f:
		nums = []
		dim = 0
		for l in f:
			if dim == 3:
				dim = 0
				nums.append(np.array([x,y,z]))
			if dim == 0:
				x = float(l)
			if dim == 1:
				y = float(l)
			if dim == 2:
				z = float(l)
			dim += 1
		print nums

	nums = np.array(nums)

#	fig = pylab.figure()
#	ax = Axes3D(fig)
#
#	ax.scatter(nums[:,0], nums[:,1], nums[:,2])
#	pyplot.show()





if __name__ == "__main__":
	main()
