from random import randint


def read_generator(file_object,max_reads=10**15,kmer_size=None):
	id_line = '          NAME: '
	quality_line = '       QUALITY: '
	sequence_line = '          READ: '
	read_split_line = '    READ_START: '
	last = None
	r = 0
	while (last != file_object.tell()) and (r < max_reads):
		last = file_object.tell()
		line = read_until_line_type(file_object,id_line)
		I = line[len(id_line):].strip()
		line = read_until_line_type(file_object,quality_line)
		if line:
			Q = [int(q) for q in line[len(quality_line):].split(',')]
		else:
			break
		line = read_until_line_type(file_object,sequence_line)
		if line:
			S = line[len(sequence_line):].strip()
		else:
			break
		line = read_until_line_type(file_object,read_split_line)
		if line:
			split = [int(x) for x in line[len(read_split_line):].split(',')][1]
		else:
			break
		if kmer_size:
			for kmer in kmers_from_read(S[:split],Q[:split],kmer_size):
				if kmer['q'].count(2) > 2:
					break
				yield kmer
			for kmer in kmers_from_read(S[split:],Q[split:],kmer_size):
				if kmer['q'].count(2) > 2:
					break
				yield kmer
		else:
			yield {'_id': I+'a','s': S[:split],'q': Q[:split]}
			yield {'_id': I+'b','s': S[split:],'q': Q[split:]}
		r += 1

def read_until_line_type(f,t):
	l = f.readline()
	last = None
	while (l[:len(t)] != t) and (last != f.tell()):
		last = f.tell()
		l = f.readline()
	if last == f.tell():
		return None
	return l

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
			break
	ri = randint(0,len(rs[0]['s'])-k)
	return {'s': rs[0]['s'][ri:ri+k],'q': rs[0]['q'][ri:ri+k]}