#!/usr/bin/env python

import sys
from itertools import groupby
from operator import itemgetter

def read_mapper_output(file_object):
	for line in file_object:
		yield line.rstrip().split('\t',1)

def bin_counts():
	counts = read_mapper_output(sys.stdin)
	for current_bin,group in groupby(counts, itemgetter(0)):
		try:
			total_count = sum(int(count) for b,count in group)
			print '%s\t%d' % (current_bin,total_count)
		except ValueError:
			pass

def read_bins():
	read_bin_pairs = read_mapper_output(sys.stdin)
	for current_read,group in groupby(read_bin_pairs, itemgetter(0)):
		try:
			bin_list = set([int(b) for read,b in group])
			print '%s\t%s' % (current_read,','.join([str(b) for b in bin_list]))
		except ValueError:
			pass

if __name__ == "__main__":
	#bin_counts()
	read_bins()