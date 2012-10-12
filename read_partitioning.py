#from vdb_reader import *
from fastq_reader import *
from hyper_sequences import affine_hull,letters_to_coords,generator_to_bins,get_wheels
import cPickle
from bitarray import bitarray
import math
from ctypes import c_uint8

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

# NOT THE MOST EFFICIENT THING EVER...PASSING OVER AND HASHING ALL DATA AGAIN
def membership_generator(f,H,wheel_path='/mnt/Wheels.txt',block_size=10000):
	s = int(math.log(len(H),2))
	W = get_wheels(wheel_path,spoke_limit=s,wheel_limit=1)
	k = len(W[0]['p'])
	f.seek(0)
	last = None
	while last != f.tell():
		last = f.tell()
		try:
			A,B = generator_to_bins(read_generator(f,max_reads=block_size,verbose_ids=True,kmer_size=k),W)
			reads_found = {}
			for b in range(len(B)):
				for a in range(len(A)):
					if H[B[b][a]]:
						# A SIMPLETONS METHOD FOR ONLY WRITING EACH READ ONCE
						if A[a] not in reads_found:
							reads_found[A[a]] = True
							yield A[a]
		except Exception, err:
			print str(err)
		print f.tell()