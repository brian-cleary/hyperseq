def read_generator(file_object,max_reads=10**15,verbose_ids=False,from_hasha=False,kmer_size=None):
	last = None
	r = 0
	while (last != file_object.tell()) and (r < max_reads):
		last = file_object.tell()
		Read = read_until_new(file_object)
		I = Read[0].strip()
		if from_hasha:
			if verbose_ids:
				I = ''.join(Read[:-1])
			B = Read[-1]
			if 'k, bins:' in B:
				B = [int(c) for c in B[10:-2].split(',')]
				yield (I,B[0],B[1:])
		else:
			S = ''.join([s.strip() for s in Read[1:]])
			if kmer_size:
				for kmer in kmers_from_read(S,kmer_size):
					if verbose_ids:
						kmer['_id'] = ''.join(Read)
					yield kmer
			elif (S):
				yield {'_id': I,'s': S}
		r += 1

def read_until_new(file_object):
	last = file_object.tell()
	line = file_object.readline()
	S = [line]
	line = 'dummyline'
	while (last != file_object.tell()) and (line[0] != '>'):
		last = file_object.tell()
		line = file_object.readline()
		S.append(line)
	if last != file_object.tell():
		file_object.seek(last)
		S = S[:-1]
	return S

def kmers_from_read(s,k):
	i = 0
	while i < len(s)-k+1:
		yield {'_id': i,'s': s[i:i+k]}
		i += 1