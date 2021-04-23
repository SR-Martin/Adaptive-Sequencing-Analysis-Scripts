#!/usr/bin/python

import sys, getopt, random
import time
import errno

class SpeciesStats:
	def __init__(self, name):
		self.name 			= name
		self.unblock 		= 0
		self.no_decision 	= 0
		self.stop_receiving = 0
		self.other 			= 0
		self.duplicates		= 0


	def PrintStats(self):
		total = self.unblock + self.no_decision + self.stop_receiving + self.other
		unblock_pc = 100 * self.unblock/float(total)
		no_decision_pc = 100 * self.no_decision/float(total)
		stop_receiving_pc = 100 * self.stop_receiving/float(total)
		other_pc = 100 * self.other/float(total)
		print "Stats for " + self.name + ":"
		print "Unblock:\t" + str(self.unblock) + "\t" + str(unblock_pc) + "%"
		print "No decision:\t" + str(self.no_decision) + "\t" + str(no_decision_pc) + "%"
		print "Stop Receiving:\t" + str(self.stop_receiving) + "\t" + str(stop_receiving_pc) + "%"
		print "Other:\t" + str(self.other) + "\t" + str(other_pc) + "%"
		print "Duplicates:\t" + str(self.duplicates) + "\n"

def PrintHelp():
	print "RU_decision_stats.py -c <RU CSV file> -m <all mappings SAM file> -r <reference name for unblocked reads>"

try:
	opts, args = getopt.getopt(sys.argv[1:],"hc:m:r:")
except getopt.GetoptError:
	print "Option not recognised."
	PrintHelp()
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		PrintHelp()
		sys.exit()
	elif opt in ("-c"):
		csv_filename = arg
	elif opt in ("-m"):
		mapping_filename = arg
	elif opt in ("-r"):
		reference_name = arg

if(not csv_filename or not mapping_filename):
	PrintHelp()
	sys.exit()

read_ref_dict = dict()
ref_stats = dict()

try:
	with open(mapping_filename, 'r') as mappings:
		for line in mappings:
			if line[0] != '@':
				fields = line.split()
				reference = fields[2].strip()
				read_ref_dict[fields[0].strip()] = reference
				if reference not in ref_stats.keys():
					ref_stats[reference] = SpeciesStats(reference)

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not open file " + mapping_filename
		sys.exit(2)

num_reads_not_found = 0
reads_seen = list()
reference_unblocked = list()

try:
	with open(csv_filename, 'r') as csv_file:
		for line in csv_file:
			if line[0] != ",":
				fields = line.split(",")
				read_ID = fields[5]
				decision = fields[7].strip()
				if read_ID in read_ref_dict.keys():
					reference = read_ref_dict[read_ID]
					stats = ref_stats[reference]
					if decision == "unblock":
						stats.unblock += 1
						if reference_name and reference_name == reference:
							reference_unblocked.append(read_ID)

					elif decision == "no_decision":
						stats.no_decision += 1
					elif decision == "stop_receiving":
						stats.stop_receiving += 1
					else:
						stats.other += 1
						print "Unknown decision for read " + read_ID + ": " + decision

					if read_ID in reads_seen:
						stats.duplicates += 1
					else:
						reads_seen.append(read_ID)
				else:
					num_reads_not_found += 1
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print "Could not open file " + csv_filename
		sys.exit(2)

if reference_name:
	try:
		with open(reference_name + "_unblocked.txt", 'w') as outfile:
			for read_ID in reference_unblocked:
				outfile.write(read_ID + "\n")
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print "Could not open file " + reference + "_unblocked.txt"
			sys.exit(2)

print "Number of reads in CSV not found in SAM: " + str(num_reads_not_found)
for key in ref_stats.keys():
	ref_stats[key].PrintStats()

