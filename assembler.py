from hyper_sequences import set_wheels,get_wheels,generator_to_bins
from pymongo import Connection
from random import randint
import itertools
from collections import defaultdict

conn = Connection()

# PROBABLY FINE FOR TESTING, NOT OPTIMIZED FOR PRODUCTION PERFORMANCE

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

def read_sequence(initial_sequence,H,W):
	full_sequence = initial_sequence
	k = len(full_sequence)
	while True:
		extension = longest_path(full_sequence[-k:],H,W,True)
		if extension:
			full_sequence += extension
		else:
			break
	while True:
		extension = longest_path(full_sequence[:k],H,W,False)
		if extension:
			full_sequence = extension + full_sequence
		else:
			break
	return full_sequence

def extend_path(s0,H,W,fb):
	k = len(s0)
	A,B = generator_to_bins(extend_kmer(s0,forward=fb),W)
	E = []
	for a in range(len(A)):
		if H.get(B[0][a],(None,None))[int(fb)] == A[a]:
			E.append(A[a])
	return E

def longest_path(s0,H,W,fb):
	k = len(s0)
	P = extend_path(s0,H,W,fb)
	while len(P) > 1:
		did_extension = False
		for i in range(len(P)):
			p = P[i]
			if fb:
				sp = s0[len(p):] + p[-k:]
			else:
				sp = p[:k] + s0[:max(0,k-len(p))]
			for e in extend_path(sp,H,W,fb):
				if fb:
					P.append(p+e)
				else:
					P.append(e+p)
				did_extension = True
		P = [(len(p),p) for p in P]
		P.sort(reverse=True)
		P = [p[1] for p in P if p[0]==P[0][0]]
		if not did_extension:
			P = P[:1]
	if P:
		return P[0]
	else:
		return ''

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