#!/usr/bin/env python

import sys, getopt, errno, array

readMappingCount = dict()
refName = ""

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetTimeHist.py -s <sequencing summary file> -c <channels> -r <reference name>")
	print("python GetTimeHist.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetTimeHist.py -s <sequencing summary file> -c <channels>")
		print("Gets total sequences and total sequence length of sequences mapping to <reference name> over time.")
		print("-s <sequencing summary file>\t\t Name of file containing sequencing summary.")
		print("-c <channels> Channels to split results by e.g. 1-256.")
		print("-r <reference name> Only consider reads mapping to this reference.")
		sys.exit()
	elif opt in ("-s"):
		sequencingSummaryFilename = arg
	elif opt in ("-c"):
		channels=arg
		fields=arg.split("-")
		min_channel = int(fields[0])
		max_channel = int(fields[1])
	elif opt in ("-r"):
		refName = arg
		refName = refName.lower()

bin_size = 3600
numBins = 72 * 60 * 60 / bin_size	#3 days!
counts_target = array.array('i',(0 for i in range(0,numBins)))
lengths_target = array.array('i',(0 for i in range(0,numBins)))
counts_non_target = array.array('i',(0 for i in range(0,numBins)))
lengths_non_target = array.array('i',(0 for i in range(0,numBins)))

try:
	with open(sequencingSummaryFilename, 'r') as infile:
		count = 0
		for line in infile:
			if count > 0:
				fields = line.split()
				channel = int(fields[5])
				startTime = fields[7]
				length = int(fields[14])
				ref = fields[23].lower()
				bin_num = int(float(startTime) / bin_size)
				if refName == "" or ref == refName:
					if channel >= min_channel and channel <= max_channel:
						counts_target[bin_num] += 1
						lengths_target[bin_num] += length
					else:
						counts_non_target[bin_num] += 1
						lengths_non_target[bin_num] += length
			count += 1
			
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + sequencingSummaryFilename)
	print(e)
	sys.exit(2)

outfile_target_name = "target_times.hist"
outfile_non_target_name = "non_target_times.hist"
if refName != "":
	outfile_target_name = "target_times_" + refName + ".hist"
	outfile_non_target_name = "non_target_times_" + refName + ".hist"

with open(outfile_target_name, 'w') as outfile:
	outfile.write("time\tfrequency\tlength\n")
	for i in range(0, numBins):
		outfile.write(str(i) + "\t" + str(counts_target[i]) + "\t" + str(lengths_target[i]) + "\n")
with open (outfile_non_target_name, 'w') as outfile:
	outfile.write("time\tfrequency\tlength\n")
	for i in range(0, numBins):
		outfile.write(str(i) + "\t" + str(counts_non_target[i]) + "\t" + str(lengths_non_target[i]) + "\n")	

