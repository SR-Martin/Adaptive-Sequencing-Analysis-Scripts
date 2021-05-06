#!/usr/bin/env python

import sys, getopt, errno, array

readMappingCount = dict()
refName = ""

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:p:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetTimeHist.py -s <sequencing summary file> -c <channels> -r <reference name> -p <outfile prefix>")
	print("python GetTimeHist.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetTimeHist.py -s <sequencing summary file> -c <channels>")
		print("Gets the average over all query sequences of the number of mappings to reference.")
		print("-s <sequencing summary file>\t\t Name of file containing sequencing summary.")
		print("-c <channels> Channels to split results by e.g. 1-256.")
		print("-r <reference name> Only consider reads mapping to this reference.")
		print("-p <outfile prefix> Prefix for outfile names.")
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
	elif opt in ("-p"):
		prefix = arg

bin_size = 90 * 60
numBins = (6 * 3600) / (15 * 60)	#bin every 15 minutes

bin_sizes = [15,30,45,60,75,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90]
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
				startTime = float(fields[7])
				length = int(fields[14])
				ref = fields[23].lower()
				for i in range(numBins):
					binEnd = ((i+1) * 15) * 60
					binStart = binEnd - (bin_sizes[i] * 60)
					if startTime >= binStart and startTime < binEnd:				
						bin_num = i
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

outfile_target_name = prefix + "_target_times.hist"
outfile_non_target_name = prefix + "_non_target_times.hist"
if refName != "":
	outfile_target_name = "target_times_" + refName + ".hist"
	outfile_non_target_name = "non_target_times_" + refName + ".hist"

with open(outfile_target_name, 'w') as outfile:
	outfile.write("time\tfrequency\tmbph\n")
	for i in range(0, numBins):
		time = i * 0.25 # time in hours
		mbph = (lengths_target[i] / (bin_sizes[i]*60)) * 3600. / 1000000
		outfile.write(str(time) + "\t" + str(counts_target[i]) + "\t" + str(mbph) + "\n")
with open (outfile_non_target_name, 'w') as outfile:
	outfile.write("time\tfrequency\tmbph\n")
	for i in range(0, numBins):
		time = i * 0.25
		mbph = (lengths_non_target[i] / (bin_sizes[i]*60)) * 3600. /1000000
		outfile.write(str(time) + "\t" + str(counts_non_target[i]) + "\t" + str(mbph) + "\n")	

