hyperseq
========

Hyperplanes and so forth for sequence alignment

The goal of this project is to create a hashing structure in which kmers with very few mismatches collide. The probabilistic data srtucture defined herein will enable proper scaling of reference genome alignment tasks in which the number of permitted mismatches grows beyond 1 or 2.


BACKGROUND
Reference genome alignment tools, such as Bowtie, seem to use explicit match algorithms that inhibit scaling beyond queries with 1-2 mismatched letters. The reason for this, is that they cannot make any "fuzzy" queries - they simply query the exact sequence, then all allowed permutations until a match is found, or the permutations exhausted. As the number of permitted permutations (mismatches) grows, the query space explodes.


PROBABILISTIC QUERIES
A probabilistic setup reduces the task significantly. The query space grows with desired accuracy, and tends to be sub-linear with the number of permutations. The tradeoff is that existing matches that fit the desired permutations may not be found, since the results are probabilistic in nature.
Probabilistic data structures make use of hashing functions; each of the corpus elements are mapped (hashed) to one or more values, then the query is mapped using the same function(s), and collisions correspond to matches. However, with many hash functions collisions are random; exact matches would collide, but other collisions would not be indicative of similarity.


HYPERPLANES
Given a high-dimensional space, a series of n hyperplanes can be used to partition the space into 2^n bins. Elements can be mapped in n operations by determining if the element lies above or below each successive plane. The probability of collision grows with cosine similarity: P = 1-acos(similarity)/pi
Furthermore, the partitioning into 2^n bins can be repeated m times, with a different set of hyperplanes each time, in order to achieve a desired collision profile:

Example: n=28, m=36, k=35
sequences with 0 mismatches; approx. collision prob: 1.0
sequences with 1 mismatches; approx. collision prob: 0.940754862659
sequences with 2 mismatches; approx. collision prob: 0.578871949642
sequences with 3 mismatches; approx. collision prob: 0.2871101997
sequences with 4 mismatches; approx. collision prob: 0.138402124139
sequences with 5 mismatches; approx. collision prob: 0.0679560358736


COMPLEX ALPHABET
In order to find hyperplanes and partition the space, we need to map our alphabet to a number space so that we may perform SVDs, dot products, and so forth. I've chosen the following mapping: {'A': complex(-1,0),'T': complex(1,0),'C': complex(0,-1),'G': complex(0,1)}.
Why complex numbers? With real numbers each mapping would be a linear combination of any of the others, which is bad because we want them all to be independent. In our complex plane, A=-T and C=-G, but the other pairs are all independent, and we still get to use all our standard matrix math.
This choice also allows us to easily incorporate read probability errors (ie 'A': complex(-1,0)*1-error_prob), and to include sequences with SNPs (ie 'N': complex(0,0)).


RUNNIN THE STUFF (dirty)
>>>from hyper_sequences import *
>>>for f in <list of chr1.fa type files>:
>>>...	add_sequence_file(f)
>>>set_wheels()
>>>hyperplane_wheel(28,36)