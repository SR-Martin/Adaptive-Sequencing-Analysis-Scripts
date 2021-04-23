#!/usr/bin/env python

import sys, getopt, errno

target_channel_reads = list()
non_target_channel_reads = list()
read_lengths = dict()

class ReferenceStats:
	def __init__(self, reference):
		self.totalReads = 0
		self.totalLength = 0
		self.reference = reference


def AnalyseForChannels(_readsForChannels, _samFilename):
	ReferenceStatDict = dict()
	try:
		with open(_samFilename, 'r') as infile:
			for line in infile:
				fields = line.split()
				if fields[0] in _readsForChannels:
					reference = fields[2]
					sequenceLength = read_lengths[fields[0]]
					if reference not in ReferenceStatDict.keys():
							ReferenceStatDict[reference] = ReferenceStats(reference)
					stats = ReferenceStatDict[reference]
					stats.totalReads += 1
					stats.totalLength += sequenceLength
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print("Could not find file " + infile)
		print(e)
		sys.exit(2)

	return ReferenceStatDict

try:
	opts, args = getopt.getopt(sys.argv[1:],"hf:s:c:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python analyse_RU.py -f <fastq file> -s <sam file> -c <channels>")
	print("python analyse_RU.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("analyse_RU.py -i <inputfile>")
		print("-f <fastq file>\t\t\t\t Name of fastq file containing all reads.")
		print("-s <sam file>\t\t\t Name of sam file containing all mappings to references.")
		print("-c <channels>\t channels to separate, in the form \"a-b\".")
		sys.exit()
	elif opt in ("-f"):
		fastqFile = arg
	elif opt in("-s"):
		samFile = arg
	elif opt in ("-c"):
		channels=arg
		fields=arg.split("-")
		min_channel = int(fields[0])
		max_channel = int(fields[1])

if fastqFile == '' or samFile == '' or channels == '':
	print("You must specify -f and -s and -c")
	sys.exit(2)

#split the reads by channel
try:
	with open(fastqFile, 'r') as infile:
		count = 0
		for line in infile:
			if count % 4 == 0:
				assert(line[0]=='@')
				fields = line.split()
				readName = fields[0][1:]
				channel = int(fields[3].split("=")[1])
				assert(fields[3].split("=")[0] == "ch")
				assert(channel >= 0 and channel <= 512)
				if channel >= min_channel and channel <= max_channel:
					target_channel_reads.append(readName)
				else:
					non_target_channel_reads.append(readName)
			if count % 4 == 1:
				length = len(line.strip())
				read_lengths[readName] = length
			count += 1
		
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + fastqFile)
	print(e)
	sys.exit(2)

dictionary = AnalyseForChannels(target_channel_reads, samFile)
print("Reference stats for channels " + str(channels) + ": ")
totalMapped = 0
for ref in dictionary.values():
		print( ref.reference + "\t" + str(ref.totalReads) + "\t" + str(ref.totalLength) )
		totalMapped += ref.totalReads
print("Total number of reads mapped: " + str(totalMapped) + "/" + str(len(target_channel_reads)))


dictionary2 = AnalyseForChannels(non_target_channel_reads, samFile)
print("\nReference stats for all other channels: ")
totalMapped = 0
for ref in dictionary2.values():
		print( ref.reference + "\t" + str(ref.totalReads) + "\t" + str(ref.totalLength) )
		totalMapped += ref.totalReads
print("Total number of reads mapped: " + str(totalMapped) + "/" + str(len(non_target_channel_reads)))


