def extract_confs(xyzfile):
	"""
	Takes in a regular xyz file corresponding to a sequence of same size configurations
	returns all the conformations as a list
	"""
	with open(xyzfile, "r") as r:
		lines = r.readlines()
	all_confs = list(chunks(lines, int(lines[0])+2))
	
	return random.choices(all_confs,k=1000)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_partial_charges_for_conf(conf):
	"""
	Uses xtb to compute the partial charges of a conformation
	"""
	with open("tmp.xyz", "w") as w:
		w.writelines(conf)
	cmd = "/home/macenrola/xtb-6.3.3/bin/xtb tmp.xyz -c 1 -u 0 -P 1"
	res=subprocess.check_output(cmd.split())
	with open("tmp.out", "wb") as w:
		w.write(res) 
	with open("tmp.out", "r") as r:
		lines = r.readlines()
		for i,line in enumerate(lines):
			if "number of atoms" in line:
				n = int(line.split()[-1])
			if "#   Z          covCN         q      C6AA      α(0)" in line:
				charges = [x.split()[4] for x in lines[i+1:i+1+n]]
				return charges
			
				

			

def write_to_xlsx(confs, name):
	masses = {'C':12, 'N':14, 'O':16, 'H':1}
	radii = {'C':1.7, 'N':1.55, 'O':1.52, 'H':1.1}
	sym, x, y, z, r, q, overview, m = 0, 1, 2, 3, 4, 5, 6, 7
	totmass = 0
	workbook = xlsxwriter.Workbook(name)
	for k,conf in enumerate(confs):
		print("{}/{}".format(k, len(confs)))
		worksheet = workbook.add_worksheet()
		charges = get_partial_charges_for_conf(conf)
		totmass=0
		for i in range(len(conf[2:])):
			parts = conf[2+i].strip().split()
			worksheet.write(i, sym, parts[0])
			worksheet.write(i, x, float(parts[1]))
			worksheet.write(i, y, float(parts[2]))
			worksheet.write(i, z, float(parts[3]))
			worksheet.write(i, r, radii[parts[0]])
			worksheet.write(i, q, float(charges[i]))
			worksheet.write(i, m, masses[parts[0]])
			totmass += masses[parts[0]]
		worksheet.write(0, overview, "TOTAL z")
		worksheet.write(1, overview, 1)
		worksheet.write(2, overview, "Totalmass")
		worksheet.write(3, overview, totmass)
		worksheet.write(4, overview, "Crosssection")
	workbook.close()

def write_confs(confs, name):
	with open(name, "w") as w:
		for c in confs:
			w.writelines(c)

	
if __name__ == "__main__":
	import random
	import xlsxwriter
	import subprocess
	import glob
	flist = glob.glob("/media/macenrola/HUGO_+6590843639/HEPTYLAMINE/*xyz")
	for i,f in enumerate(flist):
		print("{}/{} is {}".format(i, len(flist),f))
		confs = extract_confs(f)
		write_confs(confs, f+"-rdconfs.xyz")
		write_to_xlsx(confs, f+"-1k.xlsx")
	#print(get_partial_charges_for_conf("tmp.xyz"))
