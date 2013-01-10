#from vdb_reader import *
from fastq_reader import *
from hyper_sequences import affine_hull,letters_to_coords,generator_to_bins,get_wheels
import cPickle
from bitarray import bitarray
import math
from ctypes import c_uint8
from collections import defaultdict
import tempfile
from multiprocessing import Pool

def set_wheels(f,num_dimensions,spokes=41,wheels=5,out_path='/mnt/'):
	Wheels = []
	for w in xrange(wheels):
		Wheels += one_wheel(f,w,spokes,num_dimensions)
	Wheels.sort()
	f = open(out_path+'Wheels.txt','w')
	cPickle.dump(Wheels,f)
	f.close()

def one_wheel(f,w,spokes,dims):
	S = []
	for s in range(spokes):
		L = rand_kmer_coords(f,dims)
		P = affine_hull(L.values())
		C = P.pop()
		S.append((w,s,P,C))
	return S

def rand_kmer_coords(f,nodes):
	N = {}
	while len(N) < nodes:
		rk = rand_kmer(f,nodes)
		if rk['q'].count(2) < 3:
			N[len(N)] = list(letters_to_coords(rk))
	return N

def write_hashed_reads(read_file,out_file,s,w=1,wheel_path='/mnt/Wheels.txt',block_size=10000):
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=w)
	k = len(W[0]['p'])
	hash_prefix = 'k, bins: '
	read_file.seek(0)
	last = None
	while last != read_file.tell():
		last = read_file.tell()
		try:
			A,B = generator_to_bins(read_generator(read_file,max_reads=block_size,verbose_ids=True,kmer_size=k),W)
			# WRITING JUST ONE WHEEL HERE, ASSUMING SORTED BY READ
			B0 = []
			last_a = None
			for a in range(len(A)):
				if A[a] != last_a:
					if B0:
						out_file.write(last_a+hash_prefix+str([k] + B0)+'\n')
						B0 = []
					last_a = A[a]
				B0.append(B[0][a])
		except Exception, err:
			print Exception,str(err)
		print read_file.tell()

def create_kmer_hash(f,s,w=1,wheel_path='/mnt/Wheels.txt',block_size=10000,out_path='/mnt/Kmer_Hash.txt',reverse_compliments=True):
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=w)
	k = len(W[0]['p'])
	f.seek(0)
	H = bitarray(2**s)
	# APPARENTLY BITARRAY IS NOT GUARANTEED TO INITIALIZE EMPTY
	H.setall(False)
	last = None
	while last != f.tell():
		last = f.tell()
		try:
			A,B = generator_to_bins(read_generator(f,max_reads=block_size,kmer_size=k),W,rc=reverse_compliments)
			for b in range(len(B)):
				for a in range(len(A)):
					H[B[b][a]] = True
		except:
			pass
		print f.tell()
	fo = open(out_path,'wb')
	H.tofile(fo)
	fo.close()
	return H

def hash_counts_from_hashq(s,hashq_path,out_path='/mnt/Kmer_Hash_Counts.txt'):
	H = (c_uint8*2**s)()
	f = open(hashq_path,'r')
	last = None
	while last != f.tell():
		last = f.tell()
		for a in hashq_read_generator(f):
			for b in a[2]:
				H[b] = min(255,H[b]+1)
	f.close()
	f0 = open(out_path,'wb')
	f0.write(H)
	f0.close()
	return H

def create_kmer_hash_counts_fasta(f,s,w=1,wheel_path='/mnt/Wheels.txt',block_size=1,out_path='/mnt/Kmer_Hash_Counts.txt'):
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=w)
	k = len(W[0]['p'])
	H = (c_uint8*2**s)()
	last = None
	f.seek(0)
	while last != f.tell():
		last = f.tell()
		try:
			A,B = generator_to_bins(read_generator(f,max_reads=block_size,kmer_size=k),W)
			for b in range(len(B)):
				for a in range(len(A)):
					H[B[b][a]] = min(255,H[B[b][a]]+1)
		except Exception,err:
			print str(err)
	f0 = open(out_path,'wb')
	f0.write(H)
	f0.close()
	return H

def create_kmer_hash_counts(f,s,w=1,wheel_path='/mnt/Wheels.txt',out_path='/mnt/Kmer_Hash_Counts.txt',temp_file_size=5*10**5):
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=w)
	kmer_size = len(W[0]['p'])
	H = (c_uint8*2**s)()
	block = []
	last = None
	f.seek(0)
	while last != f.tell():
		last = f.tell()
		f0 = []
		f0 += f.readlines(temp_file_size)
		f0 += read_until_new(f)
		block.append((f0,kmer_size,W))
		if len(block) > 500:
			pool = Pool()
			results = pool.map(hash_count_part,block,max(1,len(block)/8))
			for result in results:
				for k,v in result.items():
					# would be great to have overflow checking
					H[k] = min(255,H[k]+v)
			pool.close()
			pool.join()
			block = []
	if block:
		pool = Pool()
		results = pool.map(hash_count_part,block,max(1,len(block)/8))
		for result in results:
			for k,v in result.items():
				# would be great to have overflow checking
				H[k] = min(255,H[k]+v)
		pool.close()
		pool.join()
		block = []
	f0 = open(out_path,'wb')
	f0.write(H)
	f0.close()
	return H
	
def hash_count_part(args):
	read_lines,k,W = args
	H = defaultdict(int)
	try:
		A,B = generator_to_bins(reads_from_string(read_lines,kmersize=k),W,rc=True)
		for b in range(len(B)):
			for a in range(len(A)):
				H[B[b][a]] += 1
	except Exception,err:
		print str(err)
	return H

def open_count_hash(file_path,s):
	f = open(file_path,'rb')
	H = (c_uint8*2**s)()
	f.readinto(H)
	f.close()
	return H

def membership_generator(f,H,GW,block_size=10000,match_thresh=0.75):
	f.seek(0)
	last = None
	while last != f.tell():
		last = f.tell()
		try:
			for a in read_generator(f,max_reads=block_size,verbose_ids=True,from_hashq=True):
				sect_sums = H[:,a[2]].sum(1)
				read_match_thresh = GW[a[2]].sum()*match_thresh
				for h in range(H.shape[0]):
					if sect_sums[h] > read_match_thresh:
						yield h,a[0]
		except Exception, err:
			print str(err)
		print f.tell()

def bitarray_from_array(A,s):
	H = bitarray(2**s)
	H.setall(False)
	for a in A:
		H[a] = True
	return H