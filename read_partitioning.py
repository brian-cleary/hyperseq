#from vdb_reader import *
from fastq_reader import *
from hyper_sequences import affine_hull,letters_to_coords,generator_to_bins,get_wheels
import cPickle
from bitarray import bitarray
import math
from ctypes import c_uint8
from collections import defaultdict

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

def create_kmer_hash(f,s,w=1,wheel_path='/mnt/Wheels.txt',block_size=10000,out_path='/mnt/Kmer_Hash.txt'):
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
			A,B = generator_to_bins(read_generator(f,max_reads=block_size,kmer_size=k),W)
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

def create_kmer_hash_counts(f,s,w=1,wheel_path='/mnt/Wheels.txt',block_size=10000,out_path='/mnt/Kmer_Hash_Counts.txt'):
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=w)
	k = len(W[0]['p'])
	f.seek(0)
	H = (c_uint8*2**s)()
	last = None
	while last != f.tell():
		last = f.tell()
		try:
			A,B = generator_to_bins(read_generator(f,max_reads=block_size,kmer_size=k),W)
			for b in range(len(B)):
				for a in range(len(A)):
					# would be great to have overflow checking
					H[B[b][a]] = min(255,H[B[b][a]]+1)
		except Exception,err:
			print str(err)
		print f.tell()
	fo = open(out_path,'wb')
	fo.write(H)
	fo.close()
	return H

def open_count_hash(file_path,s):
	f = open(file_path,'rb')
	H = (c_uint8*2**s)()
	f.readinto(H)
	f.close()
	return H

# DUMB: match_thresh SHOULD PROBABLY BE DETERMINED FROM READS
def membership_generator(f,H,wheel_path='/mnt/Wheels.txt',block_size=10000,match_thresh=35):
	s = int(math.log(len(H[0]),2))
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=1)
	f.seek(0)
	last = None
	while last != f.tell():
		last = f.tell()
		try:
			for a in read_generator(f,max_reads=block_size,verbose_ids=True,from_hashq=True):
				read_set = set(a[2])
				for h in range(len(H)):
					if len(read_set & H[h]) >= match_thresh:
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