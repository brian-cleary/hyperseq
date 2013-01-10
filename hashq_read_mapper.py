#!/usr/bin/env python

from numpy import *
import glob,os
import sys

def read_generator(file_object):
	read_strings = []
	line = file_object.readline().strip()
	while line != '':
		if line[0] == '@':
			if len(read_strings) == 5:
				try:
					I = 'NEWLINE'.join(read_strings[:-1])
					B = [int(c) for c in read_strings[-1][10:-1].split(',')]
					yield (I,B[0],B[1:])
				except Exception,err:
					pass
			read_strings = []
		read_strings.append(line)
		line = file_object.readline().strip()

def membership_generator(H,Hkeys,GW,match_thresh=0.4):
	for a in read_generator(sys.stdin):
		try:
			read_set = set(a[2])
			read_match_thresh = GW[a[2]].sum()*match_thresh
			for h in range(len(H)):
				sect_sum = GW[list(read_set & H[h])].sum()
				if sect_sum > read_match_thresh:
					print '%s\t%s' % (Hkeys[h],a[0])
		except:
			pass

if __name__ == "__main__":
	global_weights = load('global_weights.npy')
	cluster_files = glob.glob(os.path.join('./cluster_vectors','*.npy'))
	cluster_weights = []
	cluster_keys = []
	for cf in cluster_files:
		cluster_weights.append(set(load(cf)))
		cluster_keys.append(cf[cf.rfind('/')+1:cf.rfind('.')])
	membership_generator(cluster_weights,cluster_keys,global_weights)