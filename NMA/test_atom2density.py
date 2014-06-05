import atom2density
import unittest
import math


class Testatom2density(unittest.TestCase):

	def setUp(self):
		# Load in some pickle file with Xti data
		self.magicalvalue = 101


	def test_test(self):
		assert self.magicalvalue == 101

	def test_ErfContribution(self):
		val = atom2density.erfContribution([0,0,0], .1, 10)
		self.assertAlmostEqual(val, 1, delta = .0001)

		val = atom2density.erfContribution([.5,0,0], .1, 10)
		self.assertAlmostEqual(val, .5, delta = .0001)

		val = atom2density.erfContribution([-.5,0,0], .1, 10)
		self.assertAlmostEqual(val, .5, delta = .0001)

		val = atom2density.erfContribution([9.5,.5,0], .1, 10)
		self.assertAlmostEqual(val, .25, delta = .0001)

		val = atom2density.erfContribution([10.5,.5,0], .1, 11)
		self.assertAlmostEqual(val, .25, delta = .0001)

		val = atom2density.erfContribution([11,0,0], .1, 11)
		self.assertAlmostEqual(val, 1, delta = .0001)

	def test_GaussianContribution(self):
		val = atom2density.gaussianContribution([9,0,0],.5,10)
		self.assertAlmostEqual(val, math.exp(-2), delta = .000001)

		val = atom2density.gaussianContribution([.5,0,0],.5,10)
		self.assertAlmostEqual(val, math.exp(-.5), delta = .000001)

		val = atom2density.gaussianContribution([10,0,0],.5,11)
		self.assertAlmostEqual(val, math.exp(-2), delta = .000001)

		val = atom2density.gaussianContribution([-11,0,0],.5,11)
		self.assertAlmostEqual(val, math.exp(0), delta = .000001)
		
		val = atom2density.gaussianContribution([9.5,0,0],.5,10)
		self.assertAlmostEqual(val, math.exp(-2), delta = .000001)


if __name__ == '__main__':
	unittest.main()
