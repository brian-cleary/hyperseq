from random import randint

def read_generator(file_object,max_reads=10**15,verbose_ids=False,from_hashq=False,kmer_size=None):
	last = None
	r = 0
	while (last != file_object.tell()) and (r < max_reads):
		last = file_object.tell()
		line = file_object.readline()
		if line[0] == '@':
			# ASSUMING READ PAIRS ARE SPLIT INTO THEIR OWN LINES
			verbose_id = line
			I = line.split()[1]
			I = I[I.index(':')+1:]
			verbose_id += file_object.readline()
			S = verbose_id.split('\n')[-2]
			verbose_id += file_object.readline()
			verbose_id += file_object.readline()
			if verbose_ids:
				I = verbose_id
			Q = quality_code_to_int(verbose_id.split('\n')[-2])
			if kmer_size:
				for kmer in kmers_from_read(S,Q,kmer_size):
					if kmer['q'].count(2) > 2:
						break
					if verbose_ids:
						kmer['_id'] = I
					yield kmer
			elif from_hashq:
				B = file_object.readline()
				if 'k, bins:' in B:
					B = [int(c) for c in B[10:-2].split(',')]
					yield (I,B[0],B[1:])
			elif (S) and (Q):
				yield {'_id': I,'s': S,'q': Q}
			r += 1

quality_codes = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
def quality_code_to_int(code_string):
	return [quality_codes.index(c) for c in code_string]

def kmers_from_read(s,q,k):
	i = 0
	while i < len(s)-k+1:
		yield {'_id': i,'s': s[i:i+k],'q': q[i:i+k]}
		i += 1

def rand_kmer(f,k,max_seek=10**9):
	while True:
		f.seek(randint(0,max_seek))
		rs = [_ for _ in read_generator(f,max_reads=1)]
		if len(rs) > 0:
			if len(rs[0]['s']) > k:
				break
	ri = randint(0,len(rs[0]['s'])-k)
	return {'s': rs[0]['s'][ri:ri+k],'q': rs[0]['q'][ri:ri+k]}