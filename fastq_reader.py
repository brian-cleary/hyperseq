from random import randint

def hashq_read_generator(file_object,max_reads=10**15):
	last = None
	r = 0
	while (last != file_object.tell()) and (r < max_reads):
		last = file_object.tell()
		read_strings = read_until_new(file_object)
		if len(read_strings) == 5:
			try:
				I = ''.join(read_strings[:-1])
				B = [int(c) for c in read_strings[-1][10:-2].split(',')]
				yield (I,B[0],B[1:])
			except Exception,err:
				print str(err)
			r += 1

def fastq_read_generator(file_object,max_reads=10**15,verbose_ids=False,kmer_size=None):
	r = 0
	read_strings = read_until_new(file_object)
	while read_strings != ['']:
		for rx in reads_from_string(read_strings,verboseids=verbose_ids,kmersize=kmer_size):
			yield rx
		r += 1

def reads_from_string(RS,verboseids=False,kmersize=None):
	l = 0
	while l < len(RS):
		line = RS[l]
		l += 1
		if line[0] == '@':
			# ASSUMING READ PAIRS ARE SPLIT INTO THEIR OWN LINES
			verbose_id = line
			I = line.split()[1]
			I = I[I.index(':')+1:]
			#I = line[line.index(':')+1:].strip()
			verbose_id += RS[l]
			l += 1
			S = verbose_id.split('\n')[-2]
			verbose_id += RS[l]
			l += 1
			verbose_id += RS[l]
			l += 1
			if verboseids:
				I = verbose_id
			Q = quality_code_to_int(verbose_id.split('\n')[-2])
			if kmersize:
				for kmer in kmers_from_read(S,Q,kmersize):
					if kmer['q'].count(2) > 2:
						break
					if verboseids:
						kmer['_id'] = I
					yield kmer
			elif (S) and (Q):
				yield {'_id': I,'s': S,'q': Q}

def read_until_new(file_object,read_start='@'):
	last = file_object.tell()
	line = file_object.readline()
	S = [line]
	line = 'dummyline'
	while (last != file_object.tell()) and (line[0] != read_start):
		last = file_object.tell()
		line = file_object.readline()
		S.append(line)
	if last != file_object.tell():
		file_object.seek(last)
		S = S[:-1]
	return S

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