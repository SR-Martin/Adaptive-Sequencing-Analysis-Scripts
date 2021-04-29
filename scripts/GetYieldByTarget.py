#!/usr/bin/env python

import sys, getopt, errno, array

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetYieldByTarget.py -s <sequencing summary file> -c <channels> -r <reference name>")
	print("python GetYieldByTarget.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetYieldByTarget.py -s <sequencing summary file> -c <channels>")
		print("Gets the yield per hours, split by target <reference name> and total, and by channel groups -c.")
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
		refNames = arg.split(",")
		for ref in refNames:
			ref = ref.lower()

numBins = 12 #every 6 hours
binInterval = 60 * 60 * 6

target_channels_reference_yield = dict()
non_target_channels_reference_yield = dict()
for ref in refNames:
	target_channels_reference_yield[ref] = [0] * numBins
	non_target_channels_reference_yield[ref] = [0] * numBins

target_channels_total_yield = [0] * numBins
non_target_channels_total_yield = [0] * numBins

try:
	with open(sequencingSummaryFilename, 'r') as infile:
		count = 0
		for line in infile:
			if count > 0:
				fields = line.split()
				channel = int(fields[5])
				length = int(fields[14])
				reference = fields[23].lower()
				startTime = float(fields[7])
				bin_num = int(startTime / binInterval)
				if channel >= min_channel and channel <= max_channel:
					if reference in refNames:
						(target_channels_reference_yield[reference])[bin_num] += length
					else:
						target_channels_total_yield[bin_num] += length					
				else:
					if reference in refNames :
						(non_target_channels_reference_yield[reference])[bin_num] += length
					else:
						non_target_channels_total_yield[bin_num] += length
			count += 1
			
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + sequencingSummaryFilename)
	print(e)
	sys.exit(2)

outfile_target_name = "yield_target_channels.hist"
outfile_non_target_name = "yield_non_target_channels.hist"

cumulative_target_channels_reference_yield = dict()
cumulative_non_target_channels_reference_yield = dict()
for ref in refNames:
	cumulative_target_channels_reference_yield[ref] = [0] * numBins
	cumulative_non_target_channels_reference_yield[ref] =  [0] * numBins
	cumulative_target_channels_reference_yield[ref][0] = target_channels_reference_yield[ref][0]
	cumulative_non_target_channels_reference_yield[ref][0] = non_target_channels_reference_yield[ref][0]

cumulative_target_channels_total_yield = [0] * numBins
cumulative_non_target_channels_total_yield = [0] * numBins
cumulative_target_channels_total_yield[0] = target_channels_total_yield[0]
cumulative_non_target_channels_total_yield[0] = non_target_channels_total_yield[0]

for i in range(1,numBins):
	for ref in refNames:
		cumulative_target_channels_reference_yield[ref][i] = cumulative_target_channels_reference_yield[ref][i-1] + target_channels_reference_yield[ref][i]
		cumulative_non_target_channels_reference_yield[ref][i] = cumulative_non_target_channels_reference_yield[ref][i-1] + non_target_channels_reference_yield[ref][i]
	cumulative_target_channels_total_yield[i] = cumulative_target_channels_total_yield[i-1] + target_channels_total_yield[i]
	cumulative_non_target_channels_total_yield[i] = cumulative_non_target_channels_total_yield[i-1] + non_target_channels_total_yield[i]

with open(outfile_target_name, 'w') as outfile:
	header = "time\t"
	for ref in refNames:
		header += (ref + "\t")
	header += "remaining\n"
	outfile.write(header)
	for i in range(0, numBins):
		time = i 
		outString = str(time) + "\t"
		for ref in refNames:
			outString += (str(cumulative_target_channels_reference_yield[ref][i]) + "\t")
		outString += str(cumulative_target_channels_total_yield[i]) +"\n"
		outfile.write(outString)

with open(outfile_non_target_name, 'w') as outfile:
	header = "time\t"
	for ref in refNames:
		header += (ref + "\t")
	header += "remaining\n"
	outfile.write(header)
	for i in range(0, numBins):
		time = i 
		outString = str(time) + "\t"
		for ref in refNames:
			outString += (str(cumulative_non_target_channels_reference_yield[ref][i]) + "\t")
		outString += str(cumulative_non_target_channels_total_yield[i]) +"\n"
		outfile.write(outString)