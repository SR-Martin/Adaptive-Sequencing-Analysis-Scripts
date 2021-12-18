import sys, getopt, errno

try:
	opts, args = getopt.getopt(sys.argv[1:],"hi:qg",["ifile="])
except getopt.GetoptError:
	print "Option not recognised."
	print "GetCoverageBySAM.py -i <inputPAFfile>"
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print "GetCoverageBySAM.py -i <inputPAFfile>"
		sys.exit()
	elif opt in ("-i", "--ifile"):
		inputfile = arg

if inputfile == '':
	print "You must specify -i"
	print "GetCoverageBySAM.py -i <inputPAFfile>"
	sys.exit(2)

lengths = {	"Escherichia_coli_chromosome" : 4765434,
			"E.coli.b2207.final.genome"	: 5111512,
			"Escherichia_coli_B3008" : 4739263,
			"E.coli.B766.final.genome" : 5062632,
			"E.coli.JM109.final.genome" : 4497410
}

coverages = dict()	

try:
	with open(inputfile, 'r') as infile:
		for line in infile:
			fields = line.split()
			ref_name = fields[2]
			if ref_name in lengths.keys():

				ref_length = lengths[ref_name]
				ref_start = int(fields[3])
				ref_end = ref_start + len(fields[9])
				if ref_end > ref_length:
					ref_end = ref_length

				if ref_name not in coverages.keys():
					coverages[ref_name] = [0] * ref_length
					print "New ref " + ref_name + " of length " + str(ref_length)
				
				#print "Adding coverage for (" + str(ref_start) + ", " + str(ref_end) + ") on " + ref_name + "of length " + str(ref_length)
				for i in range(ref_start, ref_end):
					coverages[ref_name][i] = 1

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not find file " + inputfile
	sys.exit(2)

for key in coverages.keys():
	totalCoveredBases = 0
	for i in coverages[key]:
		if i == 1:
			totalCoveredBases += 1
	totalBases = len(coverages[key])
	coverage = float(totalCoveredBases) * 100 / totalBases

	print key + ": " + str(coverage) + "%"