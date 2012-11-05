hyperseq
========

Hyperplanes and so forth for sequence analysis

Currently, there are 3 arms in this project:
	Creating a probabilistic hashing structure via hyperplanes.
	Efficient comparative analysis by means of logical operations on multiple hashed data sets.
	Performing a "latent strain analysis" on multiple hashed samples.


PROBABILISTIC QUERIES
Probabilistic data structures herein make use of hashing functions; each of the corpus elements are mapped (hashed) to one or more values. With many hash functions collisions are random; exact matches would collide, but other collisions would not be indicative of similarity. A hash function in which similar sequences are more likely to collide is desirable.


HYPERPLANES
Given a high-dimensional space, a series of n hyperplanes can be used to partition the space into 2^n bins. Elements can be mapped in n operations by determining if the element lies above or below each successive plane. The probability of collision grows with cosine similarity: P = 1-acos(similarity)/pi


COMPLEX ALPHABET
In order to find hyperplanes and partition the space, we need to map our alphabet to a number space so that we may perform SVDs, dot products, and so forth. I've chosen the following mapping: {'A': complex(-1,0),'T': complex(1,0),'C': complex(0,-1),'G': complex(0,1)}.
Why complex numbers? With real numbers each mapping would be a linear combination of any of the others, which is bad because we want them all to be independent. In our complex plane, A=-T and C=-G, but the other pairs are all independent, and we still get to use all our standard matrix math.
This choice also allows us to easily incorporate read probability errors (ie 'A': complex(-1,0)*1-error_prob), and to include sequences with SNPs (ie 'N': complex(0,0)).


COMPARATIVE ANALYSIS
Sequence data that has been compressed with hyperplane hashing can be used directly in comparative analysis. Imagine the following hashed data sets:

	1.  0010110001011
	2.  1011010000101
	3.  0001010111010

Set (union, intersection, etc.) queries can be executed very efficiently with this data structure. For example, it is easy to see that ( (1 AND 2) NOT 3 ) is:

	    0010000000001

In practice, data are much larger, making this technique highly advantageous.


LATENT STRAIN ANALYSIS
Moving beyond logical operations, we can perform a matrix factorization on hashed data. In this case, the value at each bit in the hash represents the number of times a sequence mapped to that value. The final vector of hashed data, then, approximately represents kmer abundances.
We examine the eigenvectors (eigenkmers here) found by singular value decomposition of the sample,abundance matrix (s samples by k kmers). Clustering based on eigenkmer similarity allows us to identify kmers that belong to latent (unobserved) strains. Righteous!