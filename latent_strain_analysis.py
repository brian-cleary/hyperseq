# Latent Strain Analysis - Full Procedure

# HASHING THE RAW READS
	# via MapReduce (large data sets)
		# test like: $ head -n100 somefile.fastq | python fastq_read_mapper.py | sort -k1 | python read_reducer.py
		# upload all fastq files to S3 bucket data directory (eg s3://cleary-metagenomics/data)
		# create an elastic map reduce job:
			# input dir: s3n://cleary-metagenomics/data
			# output dir: s3n://cleary-metagenomics/hashed_reads
			# mapper: s3n://cleary-metagenomics/job/fastq_read_mapper.py
			# reducer: s3n://cleary-metagenomics/job/read_reducer.py
			# extra args: -cacheFile s3n://cleary-metagenomics/cache_files/Wheels.txt#Wheels.txt
		# with 20m reads this took 2h18m running with 19 m1.large
		# download results: s3cmd get s3://cleary-metagenomics/hashed_reads/part*
		from merge_read_bins import merge_files
		merge_files('/mnt2/mr_out/','/mnt/grinder_output/',out_prefix='/mnt/grinder_output/hashed_reads/')
		# merging 34 parts of ~600Mb each took couple hours
		from read_partitioning import hash_counts_from_hashq
		for hf in hashed_read_files:
			hash_counts_from_hashq(28,hf,out_path='/mnt/grinder_output/'+hf[hf.rfind('/')+1:hf.index('.')]+'Kmer_Hash_Counts.txt')
	# on a single node (small data sets)
		# read_partitioning: create_kmer_hash_counts (should write hashq!)

# CONDITION, SVD, AND CLUSTER KMER ABUNDANCE MATRIX
	from eigenhashes import *
	M = matrix_from_file_paths(Kmer_Hash_Count_Files,28)
	M = abundance_to_conditioned_nonzeros(M,out_prefix='/mnt/grinder_output/')
	ER = eigenkmers(M,num_dims=10,out_prefix='/mnt/grinder_output/')
	M = None
	C = kmer_clusters(ER)
	NZ = load('nonzeros.npy')
	save_clusters(C,NZ,out_prefix='/mnt/grinder_output/clustered_reads/')

# WRITE FULL READS INTO CLUSTERS
	# via MapReduce (large data sets)
		# upload all hashq files to S3 (eg s3://cleary-metagenomics/hashq-data)
		# upload global_weights to S3 (eg s3://cleary-metagenomics/cache_files/global_weights.npy')
		# create a mapreduce job for bundles of clusters:
			# tar -xcvf clusters0-4.tgz cluster0.npy cluster1.npy cluster2.npy cluster3.npy cluster4.npy
			# upload cluster tarball (eg s3://cleary-metagenomics/cache_files/clusters0-4.tgz)
			# input dir: s3n://cleary-metagenomics/hashq-data
			# output dir: s3n://cleary-metagenomics/clusters0-4
			# mapper: s3n://cleary-metagenomics/job/hashq_read_mapper.py
			# reducer: NONE
			# extra args: -cacheFile s3n://cleary-metagenomics/cache_files/global_weights.npy#global_weights.npy -cacheArchive s3n://cleary-metagenomics/cache_files/clusters0-4.tgz#cluster_vectors
			# Bootstrap Memory Intensive Configuration
			# download results: s3cmd get s3://cleary-metagenomics/clusters0-4/part*
			from cluster_assembly_and_alignment import fastq_from_mr_output
			fastq_from_mr_output('/mnt2/mr_output/',out_prefix='/mnt/grinder_output/clustered_reads/')
	# on a single node (small data sets)
		# eigenhashes: write_reads_from_clusters (will need hashq!)

# ANALYZE CLUSTER CONTENT
	from cluster_assembly_and_alignment import process_cluster, cluster_summary_stats
	for c in all_clusters:
		process_cluster(c,'/mnt/grinder_output/clustered_reads/','/mnt/grinder_output/clustered_reads/')
	cluster_summary_stats