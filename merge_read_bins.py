import glob,os
from fastq_reader import read_until_new
from ctypes import c_uint8
from collections import defaultdict

def merge_files(bin_path,read_path,out_prefix='/mnt/'):
	hash_prefix = 'k, bins: '
	bin_files = glob.glob(os.path.join(bin_path,'part-*'))
	read_files = glob.glob(os.path.join(read_path,'*.fastq'))
	for bf in bin_files:
		Read_Bin_Index = index_bin_file(bf)
		for rf in read_files:
			f = open(rf,'r')
			g = open(out_prefix+rf[rf.rfind('/')+1:rf.rfind('.')]+'.hashq','a')
			last = None
			while last != f.tell():
				last = f.tell()
				read_strings = read_until_new(f)
				if len(read_strings) > 1:
					if read_strings[0].strip().split()[0] in Read_Bin_Index:
						g.writelines(read_strings)
						# SUPER DUMB: HARDCODE WRITING KMER SIZE
						g.write(hash_prefix+'[35,'+Read_Bin_Index[read_strings[0].strip().split()[0]][:-1]+']\n')
			f.close()
			g.close()
			print bf,rf

def merge_into_hash(read_prefix,bin_prefix,s,out_path):
	H = (c_uint8*2**s)()
	bin_files = glob.glob(os.path.join(bin_path,'part-*'))
	read_files = [read_prefix+'.short.singleton.fastq',read_prefix+'.short.pairs.fastq']
	for bf in bin_files:
		Read_Bin_Index = index_bin_file(bf)
		for rf in read_files:
			f = open(rf,'r')
			last = None
			while last != f.tell():
				last = f.tell()
				read_strings = read_until_new(f)
				if read_strings:
					if read_strings[0].strip() in Read_Bin_Index:
						read_bins = [int(rb) for rb in Read_Bin_Index[read_strings[0].strip()].split(',')]
						for rb in read_bins:
							H[rb] = min(255,H[rb]+1)
			f.close()
	f0 = open(out_path,'wb')
	f0.write(H)
	f0.close()
	return H

def index_bin_file(fp):
	f = open(fp,'r')
	#last = None
	Index = defaultdict(str)
	#while last != f.tell():
	#	last = f.tell()
	#	line = f.readline().strip().split()
	#	if line:
	#		Index[line[0]] = line[1]
	Lines = f.readlines()
	f.close()
	for l in Lines:
		line = l.strip().split()
		if line:
			Index[line[0]] += line[1]+','
	Lines = None
	return Index