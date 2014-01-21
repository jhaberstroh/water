from operator import mul
import numpy as np
import os
import time
import sys
import re

def printFileInfo(fpath = "data/ProductionMD/md_traj.gro"):
	fsize = os.path.getsize(fpath)
	print "File size:", fsize, "bytes (~", fsize / 1000000 , "MB )"
	
	N = 0
	with open(fpath, 'rb') as f:
    		f.readline()
    		N = int(f.readline())
	
	print "Number of particles:", N

	print "[L]: nm"
	print "[T]: ps"


def loadAtomFromFile(n_time_slice, fpath = "data/ProductionMD/md_traj.gro"):
	Xti = []
	Vti = []
	first_line = True
	first_pass = True
	
	fsize = os.path.getsize(fpath)
	
	N = 0
	with open(fpath, 'rb') as f:
    		f.readline()
    		N = int(f.readline())
	
	x_i = np.ndarray(shape=(N,3), dtype=float)
	v_i = np.ndarray(shape=(N,3), dtype=float)
	type_i = []
	sys_size_t = []
	times = []
	i = 0
	
	print "Loading Xti, Vti, type_i, and times into memory, assuming fixed number of particles...",
	sys.stdout.flush()
	with open(fpath, 'rb') as f:
    		header_flag = 0
    		N_particles = 0
    		system_size = 0
    		count = 0
    		for l in f:
        		
        		count = count+1
        		if count % (4234 * 20) == 0:
            			print float(f.tell()) / float(fsize) * 100, "% progress"
            		
        		if re.match(".*SOL.*", l):
            			#print l.split()
            			x = []
            			v = []
            			if first_pass:
                			type_i.append(l.split()[1])
                			
            			for val in l.split()[3:6]:
                			x.append(float(val))
            			for val in l.split()[6:]:
                			v.append(float(val))
            			x_i[i][:] = x
            			v_i[i][:] = v
            			i = i+1
			elif re.match(".*t=.*", l):
            			#print float(re.split("t=",l)[1])
            			i = 0
            			header_flag = 1
	    			print len(Xti), len(times)
				assert len(Xti) + 1 == len(times) or first_line
	    			if len(times) > n_time_slice:
	   	 			break
            			times.append(float(re.split("t=",l)[1]))
                     		
            			if not first_line:
                			first_pass = False
                			Xti.append(x_i)
                			Vti.append(v_i)
                			sys_size_t.append(system_size)
				elif first_line:
					first_line = False
			elif header_flag == 1:
				header_flag = 2
            			N_particles = int(l)
			elif header_flag == 2:
            			header_flag = 0
            			system_size = float(l.split()[0])
        		else:
            			print "Warning: bad line"
            			print l,
            			print header_flag
                           		
	Xti.append(x_i)
	Vti.append(v_i)
	sys_size_t.append(system_size)
		
	assert len(Xti) == len(Vti)
	assert len(Xti) == len(times)
	assert len(Xti[0]) == len(type_i)
	
	print "Done."
	print "Number of time points:", len(Xti)
	
	Xti = np.array(Xti)
	Vti = np.array(Vti)
	sys_size_t = np.array(sys_size_t)
	
	return Xti, Vti, sys_size_t, type_i


def loadOneAtomTimeFromFile(time_slice_index, fpath = "data/ProductionMD/md_traj.gro"):
	"""
	Grabs one timeslice from a .gro file, as zero-based indexed by time_slice_count.
	"""
	
	fsize = os.path.getsize(fpath)
	
	N = 0
	with open(fpath, 'rb') as f:
    		f.readline()
    		N = int(f.readline())
	
	x_i = np.ndarray(shape=(N,3), dtype=float)
	v_i = np.ndarray(shape=(N,3), dtype=float)
	times = []
	i = 0
	
	print "Loading a single Xi, Vi, type_i, into memory, using timeslice "+ str(time_slice_index)+ "...",
	sys.stdout.flush()
	with open(fpath, 'rb') as f:
    		header_flag = 0
    		N_particles = 0
    		system_size = 0
		type_i = []
    		linecount = 0
    		for l in f:
        		linecount = linecount+1
            		
        		if re.match(".*SOL.*", l) and len(times)-1 == time_slice_index:
				# This section only runs through once, when the time_slice_index matches up.
            			x = []
            			v = []
                		type_i.append(l.split()[1])
                			
            			for val in l.split()[3:6]:
                			x.append(float(val))
            			for val in l.split()[6:]:
                			v.append(float(val))
            			x_i[i][:] = x
            			v_i[i][:] = v
            			i = i+1
			elif re.match(".*t=.*", l):
            			#print float(re.split("t=",l)[1])
            			i = 0
            			header_flag = 1

				# Check that the program did not skip past...
				assert len(times)-1 <= time_slice_index
	    			if len(times)-1 == time_slice_index:
	   	 			break
				# Check the length above, because incrementation happens below.
            			times.append(float(re.split("t=",l)[1]))

			elif header_flag == 1 and len(times)-1 == time_slice_index:
				header_flag = 2
            			N_particles = int(l)
			elif header_flag == 2 and len(times)-1 == time_slice_index:
            			header_flag = 0
            			system_size = float(l.split()[0])
        		else:
				if len(times)-1 == time_slice_index:
            				print "Warning: bad line encountered"
            				print l,
            				print header_flag

	if (len(times)-1 != time_slice_index):
		raise ValueError("time_slice_index "+str(time_slice_index)+" was not found. It was probably too large.")
	assert len(x_i) == len(type_i)
	
	x_i = np.array(x_i)
	v_i = np.array(v_i)
	
	print "Time slice loaded."
	return x_i, v_i, system_size, type_i
