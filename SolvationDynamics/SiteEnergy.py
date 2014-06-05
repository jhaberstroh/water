#!/usr/bin/env python2.7
import shlex
import numpy as np

f = open('fixed.gro', 'r')
smycharge = 5.5;
tip4pHcharge = 0.52;
tip4pMcharge = -1.04;
tip4pOcharge = 0;

PLUSpos = 0;
MINpos = 0;
elist = list() 
count = 0

while (f.readline()):
	print count
	num = int(f.readline())
	if (num == 0 or count > 80):
		break
	mpos = list()
	hpos = list()

	for _ in range(num):
		line = f.readline()
		tokens = shlex.split(line);
		tokens[3] = float(tokens[3]);
		tokens[4] = float(tokens[4]);
		tokens[5] = float(tokens[5]);
		if (tokens[1] == 'POS'):
			PLUSpos = np.array(tokens[3:6])
			print tokens
		if (tokens[1] == 'NEG'):
			MINpos = np.array(tokens[3:6])
			print tokens
		if (tokens[1] == 'MW'):
			mpos.append(np.array(tokens[3:6]))
		if (tokens[1] == 'HW1' or tokens[1] == 'HW2'):
			hpos.append(np.array(tokens[3:6]))

	energy = 0;
	for hp in hpos:
		Rplus = np.linalg.norm(PLUSpos - hp)
		Rminus = np.linalg.norm(MINpos - hp)
		energy += tip4pHcharge * smycharge * (1.0/Rplus - 1.0/Rminus)
	
	for mp in mpos:
		Rplus = np.linalg.norm(PLUSpos - mp)
		Rminus = np.linalg.norm(MINpos - mp)
		energy += tip4pMcharge * smycharge * (1.0/Rplus - 1.0/Rminus)
	elist.append(energy);
	count = count + 1 
	f.readline();

print elist
