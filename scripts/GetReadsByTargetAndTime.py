#!/usr/bin/env python

import sys, getopt, errno, array

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetReadsByTargetAndTime.py -s <sequencing summary file> -c <channels> -r <reference name>")
	print("python GetReadsByTargetAndTime.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetReadsByTargetAndTime.py -s <sequencing summary file> -c <channels>")
		print("Gets lists of readnames for target <reference name>, for each cumulative hour, and by channel groups -c.")
		print("-s <sequencing summary file>\t\t Name of file containing sequencing summary.")
		print("-c <channels> Channels to split results by e.g. 1-256.")
		print("-r <reference name> Name of target sequence.")
		sys.exit()
	elif opt in ("-s"):
		sequencingSummaryFilename = arg
	elif opt in ("-c"):
		channels=arg
		fields=arg.split("-")
		min_channel = int(fields[0])
		max_channel = int(fields[1])
	elif opt in ("-r"):
		refName = arg.lower()


target_reads_enrich = dict()
target_reads_control = dict()
for i in [1,2,3,4,5,6]:
	target_reads_enrich[i] = list()
	target_reads_control[i] = list()


try:
	with open(sequencingSummaryFilename, 'r') as infile:
		count = 0
		for line in infile:
			if count > 0:
				fields = line.split()
				readname = fields[3]
				channel = int(fields[5])
				length = int(fields[14])
				reference = fields[23].lower()
				startTime = float(fields[7])
				if reference == refName:
					if channel >= min_channel and channel <= max_channel:
						for i in [1,2,3,4,5,6,]:
							if startTime <= 60 * 60 * i:
								target_reads_enrich[i].append(readname)
					else:
						for i in [1,2,3,4,5,6]:
							if startTime <= 60 * 60 * i:
								target_reads_control[i].append(readname)
			count += 1		
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + sequencingSummaryFilename)
	print(e)
	sys.exit(2)


for i in [1,2,3,4,5,6]:
	with open("target_reads_enrich" + str(i) + "hr.txt", 'w') as outfile:
		for read in target_reads_enrich[i]:
			outfile.write(read + "\n")
	with open("target_reads_control" + str(i) + "hr.txt", 'w') as outfile:
		for read in target_reads_control[i]:
			outfile.write(read + "\n")









