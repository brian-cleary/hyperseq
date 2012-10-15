from numpy import *
from scipy.sparse.linalg import svd as svds
from scipy.sparse import csr_matrix
from read_partitioning import open_count_hash, membership_generator
from random import randint
import math
from itertools import combinations
from multiprocessing import Pool

def matrix_from_file_paths(path_list,s):
	M = []
	for p in path_list:
		H = array(open_count_hash(p,s),dtype=float32)
		if len(M):
			M = concatenate((M,[H]))
		else:
			M = [H]
	return M

def condition_matrix(M,block_size=1000000):
	rows,cols = M.shape
	log_rows = math.log(rows)
	j = 0
	block = []
	while j < cols:
		block.append((j,M[:,j],log_rows))
		j += 1
		if len(block) > block_size:
			pool = Pool()
			updates = pool.map(pooled_log_entropy,block,max(1,len(block)/10))
			for c in updates:
				for r in c:
					M[r[0],r[1]] = r[2]
			pool.close()
			pool.join()
			updates = None
			block = []
			print j
	if block:
		pool = Pool()
		updates = pool.map(pooled_log_entropy,block,max(1,len(block)/10))
		for c in updates:
			for r in c:
				M[r[0],r[1]] = r[2]
		pool.close()
		pool.join()
		updates = None
		block = []
	return M

def pooled_log_entropy(args):
	col_index,col,log_rows = args
	gf = sum(col)
	rc_updates = []
	if gf > 0:
		g_weight = 1 + sum([tf/gf*math.log(tf/gf)/log_rows for tf in col if tf>0])
		for row in range(len(col)):
			if col[row] != 0:
				rc_updates.append((row,col_index,g_weight*math.log(col[row] + 1)))
	return rc_updates

def eigenkmers(M,num_dims=10):
	U,W,Vt = svds(M,num_dims)
	return W,Vt

def random_cols(M,n):
	num_cols = M.shape[1]
	RC = []
	while len(RC) < n:
		r = randint(0,num_cols-1)
		if abs(sum(M[:,r])) > 0:
			RC.append(M[:,r])
	return RC

def kmer_clusters(M,initial_clusters=25,cluster_thresh=0.9,cluster_iters=2):
	Clusters = {}
	for v in random_cols(M,initial_clusters):
		Clusters[len(Clusters)] = {'center': v,'members': []}
	num_cols = M.shape[1]
	for _ in range(cluster_iters):
		Clusters = cluster_centers(Clusters,M)
		for i in xrange(num_cols):
			if abs(sum(M[:,i])) > 0:
				found_cluster = False
				for c in Clusters:
					center = Clusters[c]['center']
					if cos_sim(M[:,i],center) > cluster_thresh:
						Clusters[c]['members'].append(i)
						found_cluster = True
				if not found_cluster:
					Clusters[len(Clusters)] = {'center': M[:,i],'members': [i]}
	return Clusters

def cos_sim(v1,v2):
	return dot(v1,v2)/(dot(v1,v1)**.5*dot(v2,v2)**.5)

def cluster_centers(C,M,combine_thresh=.98):
	for k,v in C.items():
		if len(v['members']) > 0:
			sampled_members = array(random_cols(M[:,v['members']],min(500,len(v['members']))))
			C[k]['center'] = sampled_members.sum(axis=0)/len(v['members'])
	for c in combinations(C.keys(),2):
		if (c[0] in C) and (c[1] in C):
			if cos_sim(C[c[0]]['center'],C[c[1]]['center']) > combine_thresh:
				if len(C[c[0]]['members']) >= len(C[c[1]]['members']):
					del C[c[1]]
				else:
					del C[c[0]]
	for c in C:
		C[c]['members'] = []
	return C