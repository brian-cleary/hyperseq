from fasta_reader import read_until_new
from operator import itemgetter
from numpy import *
import glob,os
from collections import defaultdict

def fastq_from_mr_output(mr_result_dir,out_prefix='/mnt/'):
	PF = glob.glob(os.path.join(mr_result_dir,'part*'))
	F = {}
	for pf in PF:
		f = open(pf)
		L = f.readlines()
		for l in L:
			l_id,l_info = l.strip().split('\t')
			if l_id not in F:
				F[l_id] = open(out_prefix+l_id+'.short.fastq','w')
			F[l_id].writelines([s+'\n' for s in l_info.split('NEWLINE')])
		f.close()
	for f in F.values():
		f.close()

assembly_suffix = '_velvet'
def process_cluster(c,readprefix,outprefix):
	c = str(c)
	total_reads = sort_read_pairs(readprefix+c+'.short.fastq',out_prefix=readprefix)
	assemble_cluster(c,readprefix,outprefix,exp_cov=max(total_reads*100/2000000,20),cov_cutoff=max(total_reads*100/2000000/15,2))
	fragment_long_reads(outprefix+c+assembly_suffix+'/contigs.fa',out_prefix=outprefix+c+assembly_suffix+'/',cov=1,frag_length=3000)
	align_assembly(c,outprefix,outprefix)

def cluster_summary_stats(c,read_prefix,post_prefix):
	c = str(c)
	Cluster_Stats = {}
	Cluster_Stats['Reads'] = read_counts(c,read_prefix)
	Cluster_Stats['Assembly'] = assembly_stats(read_prefix+c+'.fa',post_prefix+c+assembly_suffix+'/contigs.fa')
	Cluster_Stats['Alignments'] = top_alignments(post_prefix+c+'_alignments.txt')
	return Cluster_Stats

def assembly_stats(long_path,contig_path):
	AS = {}
	try:
		f = open(long_path,'r')
		line = f.readline()
	except:
		print 'long path not found:',long_path
		line = ''
	L = []
	l = 0
	while line != '':
		if line[0] == '>':
			if l > 0:
				L.append(l)
				l = 0
		else:
			l += len(line) - 1
		line = f.readline()
	if L:
		L.sort(reverse=True)
		AS['pre assembly'] = {'median contig': L[len(L)/2],'total contigs': len(L),'total bp': sum(L),'top contigs': L[:5]}
	else:
		AS['pre assembly'] = None
	L = []
	try:
		f = open(contig_path,'r')
		line = f.readline()
	except:
		print 'contig path not found:',contig_path
		line = ''
	while line != '':
		if line[0] == '>':
			L.append(int(line.split('_')[3]))
		line = f.readline()
	if L:
		L.sort(reverse=True)
		AS['post assembly'] = {'median contig': L[len(L)/2],'total contigs': len(L),'total bp': sum(L),'top contigs': L[:5]}
	else:
		AS['post assembly'] = None
	return AS

def read_counts(c,read_prefix):
	RC = {}
	RC['pairs (two per pair)'] = read_count(read_prefix+c+'.short.pairs.fastq',startchar='@')
	RC['singleton'] = read_count(read_prefix+c+'.short.singleton.fastq',startchar='@')
	return RC

def read_count(fp,startchar='>'):
	f = open(fp,'r')
	r = 0
	line = f.readline()
	while line != '':
		if line[0] == startchar:
			r += 1
		line = f.readline()
	f.close()
	return r

def read_count_dict(fp,startchar='>'):
	f = open(fp,'r')
	R = defaultdict(int)
	line = f.readline()
	while line != '':
		if line[0] == startchar:
			R[line.split('=')[-1].strip()] += 1
		line = f.readline()
	f.close()
	return R

def assemble_cluster(c,read_prefix,out_path,exp_cov=50,cov_cutoff=15):
	os.system('mkdir '+out_path+c+assembly_suffix)
	#os.system('~/src/velvet_1.2.08/velveth '+out_path+c+assembly_suffix+'/ 31 -fastq -short '+read_prefix+c+'.short.singleton.fastq -shortPaired '+read_prefix+c+'.short.pairs.fastq -fasta -long '+read_prefix+c+'.fragmented.fa')
	os.system('~/src/velvet_1.2.08/velveth '+out_path+c+assembly_suffix+'/ 31 -fastq -short '+read_prefix+c+'.short.singleton.fastq -shortPaired '+read_prefix+c+'.short.pairs.fastq')
	os.system('~/src/velvet_1.2.08/velvetg '+out_path+c+assembly_suffix+'/ -exp_cov '+str(exp_cov)+' -cov_cutoff '+str(cov_cutoff)+' -min_contig_lgth 1000')

def align_assembly(c,contig_prefix,out_path,blastdb='microbial_genomes.fa'):
	os.system('/mnt2/ncbi-blast-2.2.27+/bin/blastn -query '+contig_prefix+c+assembly_suffix+'/contigs.fragmented.fa -db /mnt2/ncbi-blast-2.2.27+/db/'+blastdb+' -task megablast -max_target_seqs 100 -outfmt "7 qseqid qlen sseqid pident length evalue bitscore score" -out '+out_path+c+'_alignments.txt')

def top_alignments(alignment_path):
	f = open(alignment_path,'r')
	A = defaultdict(int)
	last_hit = None
	for l in f.readlines():
		if l[0] != '#':
			q,ql,h,p,hl,e,bs,s = l.split('\t')
			if (h[:9] != last_hit) and (float(e) < 0.01):
				A[h[:9]] += int(hl)
				last_hit = h[:9]
		else:
			last_hit = None
	A = sorted(A.items(),key=itemgetter(1),reverse=True)
	return A[:3]

def sort_read_pairs(read_path,out_prefix='/mnt/'):
	f = open(read_path,'r')
	fn = read_path[read_path.rfind('/')+1:read_path.rfind('.')]
	pair_file = open(out_prefix+fn+'.pairs.fastq','w')
	singleton_file = open(out_prefix+fn+'.singleton.fastq','w')
	all_reads = []
	Lines = f.readlines()
	l = 0
	while l < len(Lines):
		if Lines[l][0] == '@':
			read_id = Lines[l]
			read_info = Lines[l+1:l+4]
			if len(read_info) == 3:
				all_reads.append((read_id,read_info))
			l += 1 + len(read_info)
		else:
			l += 1
	Lines = None
	all_reads = sorted(all_reads,key=itemgetter(0))
	r = 0
	while r < len(all_reads):
		if (len(all_reads)-r > 2) and (all_reads[r][0].strip().split()[0][:-1] == all_reads[r+1][0].strip().split()[0][:-1]):
			pair_file.writelines([all_reads[r][0]]+all_reads[r][1])
			pair_file.writelines([all_reads[r+1][0]]+all_reads[r+1][1])
			r += 2
		else:
			singleton_file.writelines([all_reads[r][0]]+all_reads[r][1])
			r += 1
	total_reads = len(all_reads)
	all_reads = None
	f.close()
	pair_file.close()
	singleton_file.close()
	return total_reads

def fragment_long_reads(read_path,out_prefix='/mnt/',cov=8,frag_length=3000):
	f = open(read_path,'r')
	fn = read_path[read_path.rfind('/')+1:read_path.rfind('.')]
	frag_file = open(out_prefix+fn+'.fragmented.fa','w')
	read_strings = ['dummy','dummy']
	while len(read_strings) > 1:
		read_strings = read_until_new(f)
		if len(read_strings) > 1:
			read_id = read_strings[0].strip()
			read_seq = ''.join([s.strip() for s in read_strings[1:]])
			if cov > 1:
				for i in range(max(cov,len(read_seq)/frag_length*cov)):
					r = random.randint(0,len(read_seq))
					rand_seq = read_seq[max(0,r-frag_length/2):min(len(read_seq),r+frag_length/2)]
					frag_file.write(read_id+'_'+str(i)+'\n')
					rx = 0
					while rx < len(rand_seq):
						frag_file.write(rand_seq[rx:rx+60] + '\n')
						rx += 60
			else:
				for i in range(0,len(read_seq),frag_length):
					seq = read_seq[i:i+frag_length]
					frag_file.write(read_id+'_'+str(i)+'\n')
					rx = 0
					while rx < len(seq):
						frag_file.write(seq[rx:rx+60] + '\n')
						rx += 60
	f.close()
	frag_file.close()