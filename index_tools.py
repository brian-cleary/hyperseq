import heapq
import tempfile
import array
from struct import *
from collections import defaultdict

wheel_index_prefix = '/mnt/index_files/'

# MAY NEED TO INCREASE THIS TO AVIOD OPEN FILE LIMITS
BYTE_READ_BUFFER = 800000

def merge_sort_bins(num_wheels):
	for wheel in xrange(num_wheels):
		f = open(wheel_index_prefix+'wheel_index_'+str(wheel)+'.txt','r')
		iters = []
		while True:
			I = sorted(bins_from_file(f))
			if len(I) == 0:
				break
			f0 = tempfile.TemporaryFile(dir='/mnt/tmp')
			flush_tmp(f0,I)
			f0.seek(0)
			iters.append(article_bins_from_file(f0))
		f.close()
		push_sorted(heapq.merge(*iters),wheel)

def bins_from_file(f):
	start_index = f.tell()/4
	a = array.array('I')
	a.fromstring(f.read(BYTE_READ_BUFFER))
	for element in a:
		start_index += 1
		yield (element,start_index)

def flush_tmp(f,I):
	T = array.array('I')
	for t in I:
		T.append(t[0])
		T.append(t[1])
	T.tofile(f)

def article_bins_from_file(f):
	while True:
		a = array.array('I')
		a.fromstring(f.read(BYTE_READ_BUFFER))
		if not a:
			break
		for i in range(0,len(a),2):
			yield (a[i],a[i+1])

def push_sorted(sorted_iter,wheel):
	f = open(wheel_index_prefix+'sorted_wheel_index_'+str(wheel)+'.txt','w')
	a = array.array('I')
	for bin_article in sorted_iter:
		a.append(bin_article[0])
		a.append(bin_article[1])
		if len(a) >= BYTE_READ_BUFFER/4:
			a.tofile(f)
			del a[:]
	if a:
		a.tofile(f)
	f.close()

# THIS METHOD ASSUMES ARTICLES ARE SORTED
# NOTE THAT THE 8TH article_id != 8
def push_index(I,wheel):
	f = open(wheel_index_prefix+'wheel_index_'+str(wheel)+'.txt','a')
	I = array.array('I',I)
	I.tofile(f)
	f.close()

# ASSUMES SORTED BINS
def wheel_bin_counts(f,article_bins):
	Bin_Lookups = defaultdict(list)
	for a,b in enumerate(article_bins):
		Bin_Lookups[b].append(a)
	Bin_Lookups = sorted(Bin_Lookups.iteritems())
	Matches = []
	bi = 0
	for b,a in article_bins_from_file(f):
		try:
			while Bin_Lookups[bi][0] < b:
				bi += 1
		except IndexError:
			break
		if b == Bin_Lookups[bi][0]:
			for m in Bin_Lookups[bi][1]:
				Matches.append((m,a))
	return Matches