from vdb_reader import *
from hyper_sequences import affine_hull,letters_to_coords,generator_to_bins,get_wheels
import cPickle
from bitarray import bitarray

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
	last = None
	while last != f.tell():
		last = f.tell()
		A,B = generator_to_bins(read_generator(f,max_reads=block_size,kmer_size=k),W)
		for b in range(len(B)):
			for a in range(len(A)):
				H[B[b][a]] = True
		print f.tell()
	fo = open(out_path,'wb')
	H.tofile(fo)
	fo.close()
	return H