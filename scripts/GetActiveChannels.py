#!/usr/bin/env python

import sys, getopt, errno, array

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetActiveChannels.py -s <sequencing summary file> -c <channels> ")
	print("python GetActiveChannels.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetActiveChannels.py -s <sequencing summary file> -c <channels>")
		print("Gets the number of active channels in each time window, split by channel groups -c.")
		print("-s <sequencing summary file>\t\t Name of file containing sequencing summary.")
		print("-c <channels> Channels to split results by e.g. 1-256.")
		sys.exit()
	elif opt in ("-s"):
		sequencingSummaryFilename = arg
	elif opt in ("-c"):
		channels=arg
		fields=arg.split("-")
		min_channel = int(fields[0])
		max_channel = int(fields[1])

binInterval = 3600  # 1 hour
numBins = 72
target_channels = [ [] for _ in range(numBins) ]
non_target_channels =  [ [] for _ in range(numBins) ]

try:
	with open(sequencingSummaryFilename, 'r') as infile:
		count = 0
		for line in infile:
			if count > 0:
				fields = line.split()
				channel = int(fields[5])
				startTime = float(fields[7])
				for i in range(numBins):
					if startTime >= (i * binInterval):
						bin_num = i
						if channel >= min_channel and channel <= max_channel:
							if channel not in target_channels[bin_num]:
								target_channels[bin_num].append(channel)
						else:
							if channel not in non_target_channels[bin_num]:
								non_target_channels[bin_num].append(channel)
			count += 1
			
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + sequencingSummaryFilename)
	print(e)
	sys.exit(2)

outfile_target_name = "target_active_channels.hist"
outfile_non_target_name = "non_target_active_channels.hist"

with open(outfile_target_name, 'w') as outfile:
	outfile.write("time\tnum_channels\n")
	for i in range(0, numBins):
		time = i * binInterval / 3600
		outfile.write(str(time) + "\t" + str(len(target_channels[i])) + "\n")
with open (outfile_non_target_name, 'w') as outfile:
	outfile.write("time\tnum_channels\n")
	for i in range(0, numBins):
		time = i * binInterval / 3600
		outfile.write(str(time) + "\t" + str(len(non_target_channels[i])) + "\n")

