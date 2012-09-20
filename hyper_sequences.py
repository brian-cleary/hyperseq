from numpy import *
import random
from pymongo import ASCENDING,DESCENDING,Connection
import cPickle
from operator import itemgetter
import math
from index_tools import push_index

conn = Connection()
GENOME = 'arabidopsis_thaliana'

def add_sequence_file(file_path):
	# STICK ALL THE 35MERS IN MONGO
	db = conn[GENOME]
	docs = db.kmer_ranges.find().sort([("n",DESCENDING)]).limit(1)
	doc = [d for d in docs]
	if doc:
		n = doc[0]['n'] + 1
	else:
		n = 0
	n0 = n
	f = open(file_path)
	info = f.readline().split('|')
	last = None
	line = ''
	while f.tell() != last:
		last = f.tell()
		lines = f.readlines(100)
		line += ''.join([l.strip() for l in lines])
		i = 0
		while i < max(0,len(line)-35):
			new_seq = line[i:i+35]
			if new_seq.count('N') > 2:
				i += new_seq.rfind('N') + 1
			elif len(new_seq) == 35:
				db.kmers.insert({"_id": n+i,"s": new_seq})
				i += 1
		line = line[i:]
		n += i
	db.kmer_ranges.insert({"n0": n0,"n": n,"info": info})

def generator_to_coords(sequence_generator):
	for s in sequence_generator:
		coords = letters_to_coords(s)
		yield s['_id'],coords

# THIS EQUATES A/T AND C/G MISMATCHES WITH NEGATIVE MATCHES, AND ALL OTHERS AS NON-MATCHES
ltc_q = {'A': complex(-1,0),'T': complex(1,0),'C': complex(0,-1),'G': complex(0,1)}
p_correct = 1-1./2000
ltc_no_q = {'A': complex(-1,0)*p_correct,'T': complex(1,0)*p_correct,'C': complex(0,-1)*p_correct,'G': complex(0,1)*p_correct}
def letters_to_coords(S):
	if 'q' in S:
		return array([ltc_q.get(l,complex(0,0)) for l in S['s']])*quality_to_prob(S['q'])
	else:
		return [ltc_no_q.get(l,complex(0,0)) for l in S['s']]

def quality_to_prob(Q):
	return array([1-10**(-q/10.) for q in Q])

def collision_probs(num_spokes=28,num_wheels=36,seq_len=35):
	for i in range(1,6):
		c = (seq_len - i*4./3)/seq_len
		p0 = 1 - math.acos(c)/math.pi
		p1 = p0**num_spokes
		p2 = 1 - (1 - p1)**num_wheels
		print 'sequences with',i,'mismatch(es); approx. collision prob:',p2

def set_wheels(num_dimensions,realm=None,spokes=41,wheels=200,out_path='/mnt/'):
	Wheels = []
	for w in xrange(wheels):
		Wheels += one_wheel(w,spokes,num_dimensions,realm)
	Wheels.sort()
	f = open(out_path+'Wheels.txt','w')
	cPickle.dump(Wheels,f)
	f.close()

def get_wheels(wheel_path,spoke_limit=999,wheel_limit=999999):
	f = open(wheel_path)
	Wheels = cPickle.load(f)
	f.close()
	Wheels = [{'w': x[0],'s': x[1],'p': x[2],'c': x[3]} for x in Wheels if (x[0] < wheel_limit) and (x[1] < spoke_limit)]
	return Wheels

def hyperplane_wheel(max_spokes,max_wheels,operation_size=100000,wheel_path='Wheels.txt'):
	db = conn[GENOME]
	Wheels = get_wheels(wheel_path,spoke_limit=max_spokes,wheel_limit=max_wheels)
	num_wheels = Wheels[-1]['w'] + 1
	num_spokes = Wheels[-1]['s'] + 1
	docs = db.kmers.find().sort([("_id",DESCENDING)]).limit(1)
	largest_id = [d for d in docs][0]['_id']
	n = 0
	while n < largest_id:
		n = bin_one_block(n,operation_size,Wheels)

def bin_one_block(start_id,size,Wheels):
	db = conn[GENOME]
	docs = db.kmers.find({"_id": {'$gte': start_id}}).sort([("_id",ASCENDING)]).limit(size)
	A,L = generator_to_bins(docs,Wheels)
	num_wheels = Wheels[-1]['w'] + 1
	num_spokes = Wheels[-1]['s'] + 1
	for w in range(num_wheels):
		D = []
		for a in xrange(len(A)):
			D.append(L[w][a])
		push_index(D,w)
	return A[-1]

def coords_to_bins(A,C,Wheels,reverse_compliments=False):
	num_wheels = Wheels[-1]['w'] + 1
	num_spokes = Wheels[-1]['s'] + 1
	pow2 = [2**j for j in range(num_spokes)]
	L = dot(C,transpose([w['p'] for w in Wheels]).conjugate())
	L -= [w['c'] for w in Wheels]
	L = int_((sign(L) + 1)/2)
	B = [dot(L[:,ws:ws+num_spokes],pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
	if reverse_compliments:
		L = dot(array(C)[:,::-1]*-1,transpose([w['p'] for w in Wheels]).conjugate())
		L -= [w['c'] for w in Wheels]
		L = int_((sign(L) + 1)/2)
		B2 = [dot(L[:,ws:ws+num_spokes],pow2) for ws in range(0,num_wheels*num_spokes,num_spokes)]
		return A,B,B2
	else:
		return A,B

def generator_to_bins(sequence_generator,Wheels,rc=False,return_terminals=False):
	C = []
	A = []
	for a,c in generator_to_coords(sequence_generator):
		if return_terminals:
			a = (a,c[0],c[-1])
		A.append(a)
		C.append(c)
	return coords_to_bins(A,C,Wheels,reverse_compliments=rc)

def one_wheel(w,spokes,dims,realm=None):
	if realm:
		db = conn[realm]
	else:
		db = conn[GENOME]
	N = db.kmers.count()
	S = []
	for s in range(spokes):
		L = pick_leaf_noloc(dims,db,N)
		P = affine_hull(L.values())
		C = P.pop()
		S.append((w,s,P,C))
	return S

def pick_leaf_noloc(nodes,db,n0):
	new_leaf = {}
	while len(new_leaf) < nodes:
		docs = db.kmers.find({"_id": {'$gt': random.randint(0,n0)}}).sort([("_id",ASCENDING)]).limit(1)
		nl = [_ for _ in generator_to_coords(docs)]
		if nl:
			new_leaf[len(new_leaf)] = list(nl[0][1])
	return new_leaf

def affine_hull(linear_system):
	# linear_system: d dimensions of n docs in this hyperplane
	for row in linear_system:
		row.append(-1)
	linear_system.append([0]*len(linear_system[0]))
	linear_system = array(linear_system)
	U,W,V = linalg.svd(linear_system)
	return list(V[-1,:])