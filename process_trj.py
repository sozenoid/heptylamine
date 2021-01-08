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
				try:
					xyz_array=parse_xyz_block(tempxyz)
					snapshot_list.append(xyz_array)
					tempxyz=[]
					nsnaps=len(snapshot_list)
					if nsnaps%10000==0: print "{}th snapshot processed".format(nsnaps)
				except: break

			tempxyz.append(line)
	print len(snapshot_list)
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



def get_real_time(metadynfile, nzeros=10, avg_window=30, shake_number=10):
	"""Will compute the accelerated time such as in plot coords
	0	 0.00000000	-484.96818530	-484.96893524	 0.47058574
	step rmsd listed_energy computed_energy kcal/mol_bias
	"""
	short_name = metadynfile.split('/')[-1]
	with open(metadynfile, 'r') as r:
		cum_time = 0
		consecutive_below_1 = 0 # consecutive bias number below 1 kcal/mol, at least nzeros consecutive
					 # times of low bias means the adamantanone is flying off
		bias_history = []
		
		for line in r:
			temp_step, temp_rmsd, _,_,temp_bias = line.strip().split()
			temp_rmsd, temp_bias = float(temp_rmsd), float(temp_bias)
			bias_history.append(temp_bias)
			if temp_bias<1:
				consecutive_below_1+=1
			else:
				consecutive_below_1=0
			if consecutive_below_1==nzeros and temp_rmsd>20:
			#if temp_rmsd>50:
				print "CONVERGED: {0} at step {1} with {2:e}".format(short_name, temp_step, cum_time)
				return cum_time
			avg_bias = sum(sorted(bias_history[-avg_window:])[:-shake_number])/(avg_window-shake_number)
			#avg_bias = bias_history[-1]
			temp_time = np.exp(avg_bias/1.987E-3/300)*50e-15
			cum_time += temp_time
			#print line.strip(), "{0:2.2f}".format(np.exp(temp_bias/1.987E-3/300)*50e-15)
		#print "NOT CONVERGED: {0} at step {1} with {2:e}".format(short_name, temp_step, cum_time)
		return None
			

def compute_energy_of_snap(atm_list, snap, tmpfile, charge=0):
	"""
	will write out the snap, compute its energy using xtb
	"""
	with open(tmpfile, "w") as w:
		w.write(str(len(atm_list))+"\n\n")
		for at, xyz in zip(atm_list, snap):
			w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\n".format(at, xyz[0], xyz[1], xyz[2]))

	output = subprocess.check_output("{2} {0} -c {1} -u 0 --namespace {0} --alpb water -P 1".format(tmpfile, charge, xtb_path).split())

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
		rmsds.append(np.average(np.sum((snaps1[0][:151]-els[:151])**2) )**.5)
		if rmsds[-1]<100:
			print i
			energy_comp.append(compute_energy_of_snap(atm_list1, els, xyzfname+'-TMP.xyz',charge=charge))
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
			w.write("{0}\t{1: 3.8f}\t{2: 3.8f}\t{3: 3.8f}\t{4: 3.8f}\n".format(i, els[0],
			els[1], els[2], (els[1]-els[2])*627.5))
			
			
			
def compute_avg_time_from_fit(timedist):
	"""
	PRE: Takes in a range of times assumed to belong to an exponential distribution
	POST: fits the experimental cdf of the timedist to 1-exp(-t/tau) to find tau
	"""
	unfilled_markers = [m for m, func in Line2D.markers.items()
                    if func != 'nothing' and m not in Line2D.filled_markers]
	sortedist = sorted(timedist)
	p = 1. * np.arange(len(sortedist)) / (len(sortedist) - 1)
	#
	def f_to_fit(x, scale):
		return 1-np.exp(-x/scale)
	#
	xn = np.geomspace(sortedist[0], sortedist[-1], 200)
	p0 = np.median(sortedist)
	popt, pcov = optimize.curve_fit(f_to_fit, sortedist, p, p0=p0)
	#
	plt.semilogx(sortedist, p, marker = random.choice(unfilled_markers))
	#plt.semilogx(xn, f_to_fit(xn, popt))
	plt.xlabel('Time to escape [s]')
	plt.ylabel('Empirical CDF [-]')
	

	return popt, sortedist, p

def KS_test(timedist, scale):
	"""
	PRE: a sequence of escape times from the botton energy well (computed from the MTD time and the bias potential)
	POST: Returns wether this time distribution follows a Poisson distribution AKA the law of rare events
	"""
	print stats.kstest(rvs=timedist, cdf='expon', args=(0,scale), N=len(timedist))
	
def plot_a_time_dist_file(fname):
	with open(fname, 'r') as r:
		time_dist = []
		for line in r:
			time_dist.append(float(line.strip().split()[-1]))
		#shuffle(time_dist)
		time_dist=sorted(time_dist)[:]
		popt, sortedist, p = compute_avg_time_from_fit(time_dist)
		print popt
		KS_test(time_dist, popt)
		return popt[0]

def plot_a_bunch(flist):
	times = []
	for f in flist:
		print f
		times.append(float(plot_a_time_dist_file(f)))
	cas_dic = {"CB7_ADA.xtb-all":"No cap", "cb7-ada-trop-all":"Tropylium", "HsZNTKkSocNKX9ucNJBB9I3Ne7JTBHwD-CB-0-all":"CAS_53485-21-5", "fluoropyridinium-all":"CAS_140623-89-8", "mS9F1cTE4JQM7ya5kPTvEWRxtqWkwlqR-CB-0-all":"CAS_273-81-4", "iopamidol-all":"CAS_60166-93-0", "iodoxamic_slow-MINUS-all":"CAS_31127-82-9", "pentabromofluorobenzene-all":"CAS_827-05-4", "acetrizoic_acid-MINUS-all":"CAS_85-36-9", "hexabromobenzene-all":"CAS_87-82-1","hexaiodobenzene-all":"CAS_608-74-2", "sodium_ditrizoate-MINUS-all":"CAS_737-31-5"}
	#plt.legend([r'{0:<16} $\tau$={1:} s'.format(cas_dic[x.split('/')[-1]],"{0:.3e}".format(y)) for x, y in zip(flist,times)])
	plt.legend([r'{0:<16} $\tau$={1:} s'.format(x.split('/')[-1],"{0:.3e}".format(y)) for x, y in zip(flist,times)])
	matplotlib.rc('xtick', labelsize=30) 
	matplotlib.rc('ytick', labelsize=30) 
	plt.tight_layout()
	plt.show()
		
if __name__ == "__main__":
	import numpy as np
	import sys, os
	import glob
	import subprocess
	import matplotlib 
	#matplotlib.use('Agg')
	from matplotlib import pyplot as plt
	from scipy import stats, optimize
	import scipy
	import random 
	from matplotlib.lines import Line2D



	xtb_path='/home/macenrola/xtb-6.3.3/bin/xtb'
	xtb_path='/home/uccahcl/xtb-6.3.3/bin/xtb'
	#get_RMSD_and_bias(sys.argv[1], sys.argv[2])
	#get_real_time(sys.argv[1])
	#plot_a_time_dist_file(sys.argv[1])
	#plot_a_bunch(['/home/macenrola/Desktop/test/ALL_STEPS/TIMES/{}'.format(f) for f in ["CB7_ADA.xtb-all", "cb7-ada-trop-all", "HsZNTKkSocNKX9ucNJBB9I3Ne7JTBHwD-CB-0-all", "fluoropyridinium-all", "mS9F1cTE4JQM7ya5kPTvEWRxtqWkwlqR-CB-0-all", "iopamidol-all", "iodoxamic_slow-MINUS-all", "pentabromofluorobenzene-all", "acetrizoic_acid-MINUS-all", "hexabromobenzene-all","hexaiodobenzene-all", "sodium_ditrizoate-MINUS-all"]])
	plot_a_bunch(['/home/macenrola/Desktop/test/controls/EPS_DATA_OUT/CONTROLS/{}'.format(f) for f in ["C6_H11_NH3p-all","C6_H11_NMe3p-all","C6_H11_NMe2Hp-all","C6_H11_NMeH2p-all","cyclobutylmethylamine-cb6-all","cyclopentylmethylamine-cb6-all","cyclohexylmethylamine-cb6-all","cyclopropylmethylamine-cb6-all","oxiplatin_CB7-all"]])



