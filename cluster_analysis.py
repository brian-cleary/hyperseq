import pylab
from numpy import *
from read_partitioning import create_kmer_hash_counts,open_count_hash
from hcluster import pdist, linkage, dendrogram
import glob,os
from collections import defaultdict


def condition_and_map_clusters(eigenvectors,nonzeros,global_weights,file_paths,completed_count_files,s,out_prefix='/mnt/'):
	M = []
	for fp in file_paths:
		cluster_id = fp[fp[:fp.index('_velvet')].rfind('/')+1:fp.index('_velvet')]
		outfile = out_prefix+cluster_id+'.txt'
		if outfile in completed_count_files:
			H = open_count_hash(outfile,s)
		else:
			f = open(fp,'r')
			#H = create_kmer_hash_counts(f,s,block_size=15000,out_path=outfile,temp_file_size=5*10**5)
			H = create_kmer_hash_counts_fasta(f,s,block_size=1,out_path=outfile)
		H = array(H,dtype=float32)[nonzeros]
		# THIS IS ASKING FOR BAD MAMAJU - NOT USING THE ORIGINAL FUNCTION TO CONDITION...
		H = log(H + 1)*global_weights
		H = dot(H,eigenvectors)
		if len(M):
			M = concatenate((M,[H]))
		else:
			M = [H]
		print cluster_id
	return M

def plot_cluster_tree(cluster_coords,Labels=None,link_method='single',color_thresh=.25,fontsize=8):
	D = pdist(cluster_coords,'cosine')
	# SEEMS THERE MAY SOMETIME BE VERY SMALL NEGATIVE DISTANCES ie -2*10**-16
	D = abs(D)
	L = linkage(D,method=link_method,metric='cosine')
	if Labels:
		dendrogram(L,labels=Labels,orientation='left',color_threshold=color_thresh)
	else:
		dendrogram(L,orientation='left',color_threshold=color_thresh)
	pylab.title('HMP Buccal Mucosa - Latent Strain Analysis')
	pylab.xlabel('Cosine Distance')
	pylab.ylabel('Strain with the Most Alignments to Each Cluster')
	pylab.rcParams.update({'font.size': fontsize})
	pylab.show()

def do_clusters(cluster_coords,Labels=None,link_method='single',d=0.2):
	D = pdist(cluster_coords,'cosine')
	# SEEMS THERE MAY SOMETIME BE VERY SMALL NEGATIVE DISTANCES ie -2*10**-16
	D = abs(D)
	L = linkage(D,method=link_method,metric='cosine')
	F = fcluster(L,d,'distance','cosine')
	C = defaultdict(list)
	for i in range(len(F)):
		if Labels:
			C[F[i]].append(Labels[i])
		else:
			C[F[i]].append(i)
	return C

#C = [n for n in Ns if n[1][0] == 'Haemophilus parainfluenzae ATCC']
#Contigs = {}
#for c in C:
#	s = sum([v[1] for v in Contigs.values()])
#	Contigs[c[0]] = (s,c[1][1])
def get_all_contigs(name_path):
	name_files = glob.glob(os.path.join(name_path,'summary.*.txt'))
	Contigs = defaultdict(dict)
	for fp in name_files:
		f = open(fp,'r')
		lines = f.readlines()
		for n in lines:
			if n[:8] == 'Sequence':
				ns = n.strip().split('\t')
				nsw = ns[1].split()
				name = ' '.join(nsw[1:4])
				Contigs[name][nsw[0]] = (sum([v[1] for v in Contigs[name].values()]),int(ns[2]))
		f.close()
	return Contigs

def cluster_genome_seqs(cluster,contigs,start,stop):
	genome_prefix = contigs.keys()[0][:9]
	ltc = {'A': 0,'T': 1,'C': 2,'G': 3,'N': 4}
	A = {}
	for i in range(start,stop):
		A[i] = [0,0,0,0,0]
	F = glob.glob(os.path.join('/mnt2/bwa/alignments',str(cluster)+'.*.sam'))
	for cluster_path in F:
		f = open(cluster_path,'r')
		last = None
		while last != f.tell():
			last = f.tell()
			line = f.readline().split()
			if line:
				if (len(line) == 20) and (line[2][:9] == genome_prefix):
					s = 0
					while (contigs[line[2]][0] + int(line[3]) + s >= start) and (contigs[line[2]][0] + int(line[3]) + s < stop) and (s < len(line[9])):
						A[contigs[line[2]][0] + int(line[3]) + s][ltc[line[9][s]]] += 1
						s += 1
		f.close()
	return A

def position_diffs(A,start,stop,d):
	Diffs = defaultdict(list)
	for k in range(start,stop):
		m = array([A[i][k] for i in range(len(A))])
		m_sum = m.sum(0)
		if m_sum > 0:
			for i in range(len(A)):
				m_sum_i = m_sum-m[i,:]
				if (sum(A[i][k]) > 3) and (sum(m_sum_i) > 0):
					if distance.cosine((m_sum_i)/(len(A)-1),m[i,:]) > d:
						Diffs[i].append(k,m)
	return Diffs

def cluster_genome_alignment_coords(cluster,contigs):
	genome_prefix = contigs.keys()[0][:9]
	A = []
	F = glob.glob(os.path.join('/mnt/bwa/alignments',str(cluster)+'.*.sam'))
	for cluster_path in F:
		f = open(cluster_path,'r')
		last = None
		while last != f.tell():
			last = f.tell()
			line = f.readline().split()
			if line:
				if (len(line) == 20) and (line[2][:9] == genome_prefix):
					A.append(contigs[line[2]][0] + int(line[3]))
		f.close()
	return A

def plot_genome_alignments(g_length,Alignments,al_length=100):
	num_aligns = float(len(set([al[0] for al in Alignments])))
	x = []
	last = None
	i = 0
	for al in Alignments:
		if al[0] != last:
			if len(x) > 0:
				x = list(set(x))
				x.sort()
				pylab.polar(x,[1*(1 - i/num_aligns)]*len(x))
			i += 1
			x = []
			last = al[0]
			xi = int(al[1]) + al_length
		if int(al[1]) < xi:
			x.append(float(al[1])/g_length*2*pi)
		else:
			x = list(set(x))
			x.sort()
			pylab.polar(x,[1*(1 - i/num_aligns)]*len(x))
			x = []
		xi = int(al[1]) + al_length
	x = list(set(x))
	x.sort()
	pylab.polar(x,[1*(1 - i/num_aligns)]*len(x))
	all_aligns = [int(al[1]) for al in Alignments]
	all_aligns.sort()
	x = []
	xi = all_aligns[0] + al_length
	for al in all_aligns:
		if al < xi:
			x.append(float(al)/g_length*2*pi)
		else:
			x = list(set(x))
			x.sort()
			pylab.polar(x,[1.1]*len(x))
			x = []
		xi = al + al_length
	x = list(set(x))
	x.sort()
	pylab.polar(x,[1.1]*len(x))
	pylab.grid(False)
	pylab.show()

