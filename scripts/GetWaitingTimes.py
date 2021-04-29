#!/usr/bin/env python

import sys, getopt, errno, array

try:
	opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python GetTimeBetweenChannels.py -s <sorted sequencing summary file> -c <channels> -r <refName>")
	print("python GetTimeBetweenChannels.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("GetTimeBetweenChannels.py -s <sequencing summary file> -c <channels>")
		print("Gets list of waiting time between reads matching <refName> on channels <channels>, and list of waiting times between any two reads on all other channels.")
		print("-s <sequencing summary file>\t\t Name of file containing sequencing summary. Must be sorted by start time and no header.")
		print("-c <channels> Channels to split results by e.g. 1-256.")
		print("-r <refName> Name of reference sequence")
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

num_channels = 512
channel_times = [0] * num_channels
#for i in range(1,numBins):
#	target_channels[i] = array.array()
#	non_target_target_channels = array.array()
target_channels = []
non_target_channels = []

try:
	with open(sequencingSummaryFilename, 'r') as infile:
		for line in infile:
			fields = line.split()
			channel = int(fields[5]) - 1
			startTime = float(fields[7])
			endTime = startTime + float(fields[8])
			reference = fields[23].lower()
			if reference == refName:
				if channel_times[channel] != 0:
					diff = startTime - channel_times[channel]
					if channel >= min_channel and channel <= max_channel:
						target_channels.append(diff)
					else:
						non_target_channels.append(diff)
				channel_times[channel] = endTime

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + sequencingSummaryFilename)
	print(e)
	sys.exit(2)

count = 0



outfile_target_name = "target_waiting_times.hist"
outfile_non_target_name = "non_target_waiting_times.hist"

with open(outfile_target_name, 'w') as outfile:
	for time in target_channels:
		outfile.write(str(time) + "\n")
with open (outfile_non_target_name, 'w') as outfile:
	for time in non_target_channels:
		outfile.write(str(time) + "\n")