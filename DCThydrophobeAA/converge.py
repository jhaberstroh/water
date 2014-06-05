import numpy as np
import matplotlib.pylab as plt

import integrals as itg

def SetMesh(c1,c2,r1,r2,phi):
	itg.C1G = c1
	itg.C2G = c2
	itg.R1G = r1
	itg.R2G = r2
	itg.PG  = phi


def ConvergeCos(R, chi_x, chi_y):
	mesh_vals = np.arange(5,40)
	x = np.zeros(mesh_vals.shape)
	for i, cos in enumerate(mesh_vals):
		print "Mesh ",cos,"...",
		SetMesh(cos, cos, 12, 12, 6)
		x[i] = itg.Jnlmj_R(2,2,2,2,R,chi_x, chi_y)
		print "value =",x[i]
	return mesh_vals, x


R_C = 25
R_P = 6
R_MIN = 5
R_MAX = 30

def ConvergeR(R, chi_x, chi_y):
	mesh_vals = np.arange(R_MIN,R_MAX)
	x = np.zeros(mesh_vals.shape)
	for i, rmesh in enumerate(mesh_vals):
		print "Mesh ",rmesh,"...",
		SetMesh(R_C, R_C, rmesh, rmesh, R_P)
		x[i] = itg.Jnlmj_R(2,2,2,2,R,chi_x, chi_y)
		print "value =",x[i]
	return mesh_vals, x

def ConvergePhi(R, chi_x, chi_y):
	mesh_vals = np.arange(5,20)
	x = np.zeros(mesh_vals.shape)
	for i, rmesh in enumerate(mesh_vals):
		print "Mesh ",rmesh,"...",
		SetMesh(20, 20, 20, 20, rmesh)
		x[i] = itg.Jnlmj_R(2,2,2,2,R,chi_x, chi_y)
		print "value =",x[i]
	return mesh_vals, x



if __name__ == "__main__":
	chi_x, chi_y = itg.loadChi()

	cmesh, cmtx_el = ConvergeCos(8, chi_x, chi_y)
	plt.plot(cmesh,cmtx_el)
	plt.title("Converging cos mesh")
	plt.show()

	#rmesh, rmtx_el = ConvergeR(8, chi_x, chi_y)
	#plt.plot(rmesh,rmtx_el)
	#plt.title("Converging radial mesh")
	#plt.show()

	#pmesh, pmtx_el = ConvergePhi(8, chi_x, chi_y)
	#plt.plot(pmesh,pmtx_el)
	#plt.title("Converging phi mesh")
	#plt.show()


