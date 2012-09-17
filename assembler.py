from hyper_sequences import set_wheels,get_wheels,generator_to_bins
from pymongo import Connection
from random import randint
import itertools
from collections import defaultdict
from operator import itemgetter

conn = Connection()

# PROBABLY FINE FOR TESTING, NOT OPTIMIZED FOR PRODUCTION PERFORMANCE
# RUN LIKE THIS:
# >>> S = generate_test_genome(n=100000,read_size=75,k=31)
# >>> W = do_wheels(41)
# choose spokes based on memory requirements
# >>> W = get_wheels('/mnt/test_Wheels.txt',spoke_limit=30)
# >>> H = hash_test_kmers(W)
# >>> PR = read_many_sequences(H,W,31)

def generate_test_genome(n=10000,read_size=50,repeats=10,k=20):
	L = {0: 'A',1: 'T',2: 'C',3: 'G'}
	S = ''.join([L[randint(0,3)] for _ in xrange(n)])
	db = conn['test_genome']
	db.kmers.drop()
	T = 0
	for c in range(repeats):
		i = 0
		while i < len(S):
			r = randint(-read_size/5,read_size/5) + read_size
			s = S[i:i+r]
			db.reads.insert({"_id": i,"s": s})
			i += r
			l = 0
			while l < len(s)-k:
				db.kmers.insert({'_id': str(T)+s[l]+s[l+k-1],'s': s[l:l+k]})
				T += 1
				l += 1
	return S

def do_wheels(s):
	db = conn['test_genome']
	d = len(db.kmers.find_one()['s'])
	set_wheels(d,realm='test_genome',spokes=s,wheels=1,out_path='/mnt/test_')
	W = get_wheels('/mnt/test_Wheels.txt')
	return W

def hash_test_kmers(W):
	db = conn['test_genome']
	docs = db.kmers.find({},timeout=False)
	A,B = generator_to_bins(docs,W)
	H = defaultdict()
	for i in xrange(len(A)):
		H[B[0][i]] = (A[i][-2],A[i][-1])
	return H

def extend_kmer(s,forward=True):
	A = ['A','T','C','G']
	for a in A:
		if forward:
			yield {'_id': a,'s': s[1:] + a}
		else:
			yield {'_id': a,'s': a + s[:-1]}

def sequence_terminals(s,k):
	yield {'_id': k,'s': s[:k]}
	yield {'_id': k+1,'s': s[1:k+1]}
	yield {'_id': -k-1,'s': s[-k-1:-1]}
	yield {'_id': -k,'s': s[-k:]}

def update_known_paths(Hp,P,Pb,forward=True):
	if forward:
		for i in range(len(P)-3):
			if Pb[i:i+3] not in Hp:
				Hp[Pb[i:i+3]] = (P[i+3:],Pb[i+3:])
	else:
		for i in range(len(P),3,-1):
			if Pb[i-3:i] not in Hp:
				Hp[Pb[i-3:i]] = (P[:i-3],Pb[:i-3])
	return Hp

# Known_Paths may occupy more memory than anything else...need to make this more efficient
def read_many_sequences(H,W,kmer_size):
	db = conn['test_genome']
	Read_Terminals = defaultdict(list)
	Known_Paths = {}
	docs = db.reads.find({},timeout=False)
	for doc in docs:
		if len(doc['s']) > kmer_size:
			f,b = read_sequence(doc['s'],H,Known_Paths,W,k=kmer_size)
			Known_Paths = update_known_paths(Known_Paths,f[0],f[1])
			Known_Paths = update_known_paths(Known_Paths,b[0],b[1])
			for x in f[1][-2:] + b[1][:2]:
				Read_Terminals[x].append(doc['_id'])
	Read_Terminals = [(len(v),k,v) for k,v in Read_Terminals.iteritems()]
	Read_Terminals.sort(reverse=True)
	Partitioned_Reads = {}
	for rt in Read_Terminals:
		partition_counts = defaultdict(int)
		for r in rt[2]:
			partition_counts[Partitioned_Reads.get(r,None)] += 1
		partition_counts[None] = 0
		partition_counts = sorted(partition_counts.iteritems(),key=itemgetter(1),reverse=True)
		if partition_counts[0][1] > rt[0]*.25:
			p0 = partition_counts[0][0]
		else:
			p0 = rt[1]
		for r in rt[2]:
			if r not in Partitioned_Reads:
				Partitioned_Reads[r] = p0
	return Partitioned_Reads

# Hp = {(b0,b1,b2): ('atcg',[b,b,b,b])}
def read_sequence(initial_sequence,H,Hp,W,k=None):
	forward_extension = ''
	backward_extension = ''
	forward_extension_bins = ()
	backward_extension_bins = ()
	if not k:
		k = len(initial_sequence)
	while True:
		extended_sequence = initial_sequence + forward_extension
		extension = longest_path(extended_sequence[-k:],H,W,True)
		if extension:
			forward_extension += extension[0]
			forward_extension_bins += extension[1]
			if forward_extension_bins[-3:] in Hp:
				forward_extension += Hp[forward_extension_bins[-3:]][0]
				forward_extension_bins += Hp[forward_extension_bins[-3:]][1]
				break
		else:
			break
	while True:
		extended_sequence = backward_extension + initial_sequence
		extension = longest_path(extended_sequence[:k],H,W,False)
		if extension:
			backward_extension = extension[0] + backward_extension
			backward_extension_bins = extension[1] + backward_extension_bins
			if backward_extension_bins[:3] in Hp:
				backward_extension = Hp[backward_extension_bins[:3]][0] + backward_extension
				backward_extension_bins = Hp[backward_extension_bins[:3]][1] + backward_extension_bins
				break
		else:
			break
	return (forward_extension,forward_extension_bins),(backward_extension,backward_extension_bins)

def extend_path(s0,H,W,fb):
	k = len(s0)
	A,B = generator_to_bins(extend_kmer(s0,forward=fb),W)
	E = []
	for a in range(len(A)):
		if H.get(B[0][a],(None,None))[int(fb)] == A[a]:
			E.append((A[a],(B[0][a],)))
	return E

def longest_path(s0,H,W,fb):
	k = len(s0)
	P = extend_path(s0,H,W,fb)
	while len(P) > 1:
		did_extension = False
		for i in range(len(P)):
			p = P[i][0]
			if fb:
				sp = s0[len(p):] + p[-k:]
			else:
				sp = p[:k] + s0[:max(0,k-len(p))]
			for e in extend_path(sp,H,W,fb):
				if fb:
					P.append((p+e[0],P[i][1]+e[1]))
				else:
					P.append((e[0]+p,e[1]+P[i][1]))
				did_extension = True
		P = [(len(p[0]),p) for p in P]
		P.sort(reverse=True)
		P = [p[1] for p in P if p[0]==P[0][0]]
		if not did_extension:
			P = P[:1]
	if P:
		return P[0]
	else:
		return None

def test_false_positives(S,H,W,k=20):
	i = 0
	fp = 0
	fn = 0
	while i < len(S)-k:
		Ef = extend_path(S[i:i+k],H,W,True)
		if len(Ef) == 0:
			fn += 1
		elif len(Ef) > 1:
			fp += 1
		Eb = extend_path(S[i:i+k],H,W,False)
		if len(Eb) == 0:
			fn += 1
		elif len(Eb) > 1:
			fp += 1
		i += 1
	return fp,fn