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

def kmer_clusters(M,initial_clusters=25,cluster_thresh=0.9,cluster_iters=2,block_size=1000):
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
		block = []
		for i in xrange(num_cols):
			#if abs(sum(M[:,i])) > 0:
			block.append(i)
			if len(block) > block_size:
				Clusters,Centers = distance_block(block,M,Clusters,Centers,cluster_thresh)
				block = []
			if i%1000000 == 0:
				print _,i,len(Clusters)
		if len(block) > 0:
			Clusters,Centers = distance_block(block,M,Clusters,Centers,cluster_thresh)
	return Clusters

def distance_block(indices,M,C,Cm,ct):
	D = distance.cdist(transpose(M[:,indices]),Cm,'cosine')
	for i in range(D.shape[0]):
		found_cluster = False
		for j in range(D.shape[1]):
			if D[i,j] < 1-ct:
				C[j].append(indices[i])
				found_cluster = True
		if not found_cluster:
			C.append([indices[i]])
			Cm = concatenate((Cm,[M[:,indices[i]]]))
	return C,Cm

def cluster_centers(C,Cm,M,combine_thresh=.9):
	for k in range(len(C)):
		v = C[k]
		if len(v) > 0:
			sampled_members = array(random_cols(M[:,v],min(25000,len(v))))
			Cm[k,:] = sampled_members.sum(axis=0)/len(sampled_members)
	remove_clusters = {}
	D = distance.pdist(Cm,'cosine')
	i = 0
	j = 1
	for d in D:
		if j >= Cm.shape[0]:
			i += 1
			j = i+1
		if d < 1-combine_thresh:
			if len(C[i]) >= len(C[j]):
				remove_clusters[j] = True
			else:
				remove_clusters[i] = True
		j += 1
	return [[] for _ in range(len(C)-len(remove_clusters))],Cm[[i for i in range(len(C)) if i not in remove_clusters],:]

# WILL OPEN TOO MANY FILES AND BREAK WITH LOTS (THOUSANDS) OF CLUSTERS
def write_reads_from_clusters(read_files,cluster_files,s,out_prefix='/mnt',out_type='.fastq',max_hash=250000000):
	c = 0
	while c < len(cluster_files):
		Cluster_Hashes = []
		out_files = []
		i = 0
		while (i < max_hash) and (c < len(cluster_files)):
			Cluster_Hashes.append(set(load(cluster_files[c])))
			out_files.append(open(out_prefix+str(c)+out_type,'w'))
			c += 1
			i += len(Cluster_Hashes[-1])
		for r in read_files:
			F = open(r,'r')
			for cluster,read in membership_generator(F,Cluster_Hashes):
				out_files[cluster].write(read)
			F.close()
		for f in out_files:
			f.close()