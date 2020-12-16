def join_logs(f1, f2, n=200000):
	"will produce a file with the n first lines of f1 and print f2 after"
	with open(f1+"_COMPLETE", "w") as w:
		with open(f1, "r") as r:
			for i, line in enumerate(r):
				if i<n:
					w.write(line)
				else:
					break
		with open(f2,"r") as r:
			for line in r:
				w.write(line)


def return_list_of_snapshots_and_atom_list(xyzfile):
	"""
	PRE: Takes in an xyz file
	POST: Returns a list of atoms and a list of snapshots of the trajectory as arrays
	"""
	natoms = -1
	snapshot_list=[]
	with open(xyzfile, 'rb') as r:
		tempxyz=[]
		for i, line in enumerate(r):
# =============================================================================
# 			if i==1000000: break
# =============================================================================
			if i==0:
				natoms=line.strip()
			if line.strip()==natoms and tempxyz!=[]:
				if i<1000:
					atm_list=[x.strip().split()[0] for x in tempxyz[2:]]
				xyz_array=parse_xyz_block(tempxyz)
				snapshot_list.append(xyz_array)
				tempxyz=[]
				nsnaps=len(snapshot_list)
				if nsnaps%10000==0: print "{}th snapshot processed".format(nsnaps)
			tempxyz.append(line)

	return atm_list, np.array(snapshot_list)[:]

def parse_xyz_block(xyz_block):
	"""
	PRE: Takes in an xyz block as a list of lines
	POST: Will return an array of array as np.array(np.array(), np.array(),...) each sub array contains [x,y,z] coord of a specific atom
	"""
	# let go of the atom names, the order is supposed to keep that information
	coords = np.array([np.array([float(y) for y in x.strip().split()[1:]]) for x in xyz_block[2:]], dtype=np.float64) # skip the two first lines that hold atom number and energy info
	return coords

def join_xyz(f1, f2, n=2000):
	"will produce a file with the n first lines of f1 and print f2 after"
	atm_list1, snaps1 = return_list_of_snapshots_and_atom_list(f1)
	atm_list2, snaps2 = return_list_of_snapshots_and_atom_list(f2)
	assert atm_list1==atm_list2
	with open(f1+"_COMPLETE.xyz","w") as w:
		for i, els in enumerate(snaps1[:n]):
			w.write(str(len(atm_list1))+"\n")
			w.write("i = {}\n".format(i*100))
			for at, xyz in zip(atm_list1, els):
				w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\n".format(at, xyz[0], xyz[1], xyz[2]))

		for i, els in enumerate(snaps2[:]):
			w.write(str(len(atm_list1))+"\n")
			w.write("i = {}\n".format(200000+i*100))
			for at, xyz in zip(atm_list1, els):
				w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\n".format(at, xyz[0], xyz[1], xyz[2]))



def get_real_time(metadynfile):
	"Will compute the accelerated time such as in plot coords"
	pass


def compute_energy_of_snap(atm_list, snap, working_directory="/home/macenrola/heptylamine/test_good_geom/", charge=0):
	"""
	will write out the snap, compute its energy using xtb
	"""
	tmpfile = "TMP.xyz"
	with open(working_directory+tmpfile, "w") as w:
		w.write(str(len(atm_list))+"\n\n")
		for at, xyz in zip(atm_list, snap):			
			w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\n".format(at, xyz[0], xyz[1], xyz[2]))

	os.chdir(working_directory)
	output = subprocess.check_output("/home/macenrola/xtb-6.3.3/bin/xtb {0} -c {1} -u 0 --namespace {0} --alpb water -P 4".format(tmpfile, charge).split())
	
	for line in output.split("\n"):
		if "TOTAL ENERGY" in line:
			energy = float(line.split()[3])
			return energy
	


def get_RMSD_and_bias(xyzfname, charge):
	"""
	will take an xyz file and compute the RMSD for each frame
	the xyz file needs to be pre 
	"""
	atm_list1, snaps1 = return_list_of_snapshots_and_atom_list(xyzfname)
	rmsds = []
	energy_and_bias = []
	energy_comp = []
	with open(xyzfname, "r") as r:
		for line in r:
			if "energy:" in line:
				energy_and_bias.append(float(line.split()[1]))
	for i, els in enumerate(snaps1[:]):
		rmsds.append( np.average(np.sum((snaps1[0][:151]-els[:151])**2) )**.5)
		if rmsds[-1]<100:
			energy_comp.append(compute_energy_of_snap(atm_list1, els, charge=charge))
		else: 
			break
	plt.plot(rmsds, label="RMSD (Angstrom)")
	plt.plot([(x[0]-x[1])*627.5 for x in zip(energy_and_bias, energy_comp)], label="BIAS (kcal/mol)")
	plt.xlabel("FRAME (#, every 50 fs)")
	plt.legend()
	plt.savefig(xyzfname+"_IMG.eps")
	plt.close()
	with open(xyzfname+"_DATA", "w") as w:
		for i, els in enumerate(zip(rmsds, energy_and_bias, energy_comp)):
			w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\t{4: 3.8f}\n".format(i, els[0], els[1], els[2], (els[1]-els[2])*627.5))
if __name__ == "__main__":
	import numpy as np
	import sys, os
	import glob
	import subprocess
	from matplotlib import pyplot as plt

	get_RMSD_and_bias(sys.argv[1], sys.argv[2])

#	for f in glob.glob("/media/macenrola/HUGO_+6590843639/Chemts/MTD_XTB/best_mols_50_trajs/CB7_ADA.xtb/aligned-*.xyz"):
#		try:
#			get_RMSD_and_bias(f, 0)
#		except: pass

#	for f in glob.glob("/media/macenrola/HUGO_+6590843639/Chemts/MTD_XTB/best_mols_50_trajs/cb7-ada-trop/aligned-*.xyz"):
#		try:
#			get_RMSD_and_bias(f,2)
#		except: pass

#	for f in glob.glob("/media/macenrola/HUGO_+6590843639/Chemts/MTD_XTB/best_mols_50_trajs/fluoropyridinium/aligned-*.xyz"):
#		try:
#			get_RMSD_and_bias(f,2)
#		except: pass

#	for f in glob.glob("/media/macenrola/HUGO_+6590843639/Chemts/MTD_XTB/best_mols_50_trajs/HsZNTKkSocNKX9ucNJBB9I3Ne7JTBHwD-CB-0/aligned-*.xyz"):
#		try:
#			get_RMSD_and_bias(f,2)
#		except: pass

#	for f in glob.glob("/media/macenrola/HUGO_+6590843639/Chemts/MTD_XTB/best_mols_50_trajs/mS9F1cTE4JQM7ya5kPTvEWRxtqWkwlqR-CB-0/aligned-*.xyz"):
#		try:
#			get_RMSD_and_bias(f,2)
#		except: pass
	#for f in glob.glob("/media/macenrola/HUGO_+6590843639/HEPTYLAMINE/MTD_SHUT_FLIP/*00_OLD/*-pos-1.xyz")[:]:
	#	print f
		#join_xyz(f, f.replace("_OLD", ""))

	#for f in glob.glob("/media/macenrola/HUGO_+6590843639/HEPTYLAMINE/MTD_SHUT_FLIP/*00_OLD/*-COLVAR.metadynLog")[:]:
	#	print f
	#	try:
	#		join_logs(f, f.replace("_OLD", ""))
	#	except:
	#		print "ERROR : {}".format(f)
	

