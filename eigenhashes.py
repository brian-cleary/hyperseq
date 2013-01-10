from numpy import *
from scipy.sparse.linalg import svd as svds
from scipy.sparse import csr_matrix
from read_partitioning import open_count_hash, membership_generator, bitarray_from_array
from random import randint
import math
from itertools import combinations
from multiprocessing import Pool
from scipy.spatial import distance

def matrix_from_file_paths(path_list,s):
	M = []
	for p in path_list:
		H = array(open_count_hash(p,s),dtype=float32)
		if len(M):
			M = concatenate((M,[H]))
		else:
			M = [H]
	return M

def abundance_to_conditioned_nonzeros(M,out_prefix='/mnt/'):
	total_cols = M.shape[1]
	NZ = nonzero_cols(M)
	save(out_prefix+'nonzero_indices.npy',NZ)
	M = M[:,NZ]
	GWnz = global_entropy_weights(M)
	save(out_prefix+'global_nonzero_weights.npy',GWnz)
	GW = zeros(total_cols,dtype=float32)
	for i in xrange(len(NZ)):
		GW[NZ[i]] = GW[i]
	save(out_prefix+'global_weights.npy',GW)
	M = log(M + 1)*GWnz
	save(out_prefix+'conditioned_nonzeros.npy',M)
	return M

def nonzero_cols(M):
	i = 0
	Nonzeros = []
	while i < M.shape[1]:
		m = M[:,i:i+1000000].sum(0)
		for s in range(len(m)):
			if m[s] != 0:
				Nonzeros.append(i+s)
		i += 1000000
		print i
	return Nonzeros

def global_entropy_weights(M,block_size=1000000):
	rows,cols = M.shape
	log_rows = math.log(rows)
	Global_Weights = zeros(cols)
	j = 0
	block = []
	while j < cols:
		block.append((j,M[:,j],log_rows))
		j += 1
		if len(block) > block_size:
			pool = Pool()
			updates = pool.map(pooled_entropy,block,max(1,len(block)/10))
			for c in updates:
				if c:
					Global_Weights[c[0]] = c[1]
			pool.close()
			pool.join()
			updates = None
			block = []
			print j
	if block:
		pool = Pool()
		updates = pool.map(pooled_entropy,block,max(1,len(block)/10))
		for c in updates:
			if c:
				Global_Weights[c[0]] = c[1]
		pool.close()
		pool.join()
		updates = None
		block = []
	return Global_Weights

def pooled_entropy(args):
	col_index,col,log_rows = args
	gf = sum(col)
	rc_updates = []
	if gf > 0:
		g_weight = 1 + sum([tf/gf*math.log(tf/gf)/log_rows for tf in col if tf>0])
		return (col_index,g_weight)

def eigenkmers(M,num_dims=10,out_prefix='/mnt/'):
	L,V,R = svds(M,num_dims)
	save(out_prefix+'eigenleft.npy',L)
	save(out_prefix+'eigenvalues.npy',V)
	save(out_prefix+'eigenright.npy',R)
	return transpose(transpose(R)*V)

def random_cols(M,n):
	num_cols = M.shape[1]
	RC = []
	while len(RC) < n:
		r = randint(0,num_cols-1)
		if abs(sum(M[:,r])) > 0:
			RC.append(M[:,r])
	return RC

def kmer_clusters(M,initial_clusters=200,cluster_thresh=0.8,cluster_iters=3,block_size=100):
	Clusters = []
	Centers = []
	for v in random_cols(M,initial_clusters):
		Clusters.append([])
		if len(Centers) > 0:
			Centers = concatenate((Centers,[v]))
		else:
			Centers = [v]
	num_cols = M.shape[1]
	for _ in range(cluster_iters):
		Clusters,Centers = cluster_centers(Clusters,Centers,M)
		i = 0
		while i < num_cols:
			Clusters,Centers = distance_block((i,i+block_size),M,Clusters,Centers,cluster_thresh)
			i += block_size
			if i%100000 == 0:
				print _,i,len(Clusters)
	return Clusters

def distance_block(indices,M,C,Cm,ct):
	D = distance.cdist(transpose(M[:,indices[0]:indices[1]]),Cm,'cosine')
	D = D < (1 - ct)
	indices = range(indices[0],indices[1])
	for i in xrange(D.shape[0]):
		found_cluster = False
		for j in xrange(D.shape[1]):
			if D[i,j]:
				C[j].append(indices[i])
				found_cluster = True
		if not found_cluster:
			C.append([indices[i]])
			Cm = concatenate((Cm,[M[:,indices[i]]]))
	return C,Cm

def cluster_centers(C,Cm,M,combine_thresh=.8):
	for k in range(len(C)):
		v = C[k]
		if len(v) > 0:
			sampled_members = array(random_cols(M[:,v],min(75000,len(v))))
			Cm[k,:] = sampled_members.sum(axis=0)/len(sampled_members)
	remove_clusters = {}
	D = distance.pdist(Cm,'cosine')
	D = D < (1 - combine_thresh)
	i = 0
	j = 1
	for d in D:
		if j >= Cm.shape[0]:
			i += 1
			j = i+1
		if d:
			if len(C[i]) >= len(C[j]):
				remove_clusters[j] = True
			else:
				remove_clusters[i] = True
		j += 1
	return [[] for _ in range(len(C)-len(remove_clusters))],Cm[[i for i in range(len(C)) if i not in remove_clusters],:]

def save_clusters(Clusters,Nonzeros,out_prefix='/mnt/')
	for i in range(len(Clusters)):
		c = [Nonzeros[x] for x in Clusters[i]]
		save(out_prefix+str(i)+'.npy',c)

def write_reads_from_clusters(read_files,cluster_files,s,Weights,out_prefix='/mnt',out_type='.fastq',max_hash=5,kmer_match_thresh=0.4):
	c = 0
	while c < len(cluster_files):
		Cluster_Hashes = []
		out_files = []
		i = 0
		while (i < max_hash) and (c < len(cluster_files)):
			cluster_indices = load(cluster_files[c])
			new_cluster_weights = zeros(2**28,dtype=float32)
			for ci in cluster_indices:
				new_cluster_weights[ci] = Weights[ci]
			Cluster_Hashes.append(new_cluster_weights)
			new_cluster_weights = None
			out_files.append(open(out_prefix+str(c)+out_type,'w'))
			c += 1
			i += 1
		Cluster_Hashes = array(Cluster_Hashes)
		for r in read_files:
			F = open(r,'r')
			for cluster,read in membership_generator(F,Cluster_Hashes,Weights,block_size=20000,match_thresh=kmer_match_thresh):
				out_files[cluster].write(read)
			F.close()
		for f in out_files:
			f.close()