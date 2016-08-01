# kate: syntax Python;
"""
Quality trimming.
"""

def quality_trim_index(str qualities, int cutoff_front, int cutoff_back, int base=33):
	"""
	Find the positions at which to trim low-quality ends from a nucleotide sequence.
	Return tuple (start, stop) that indicates the good-quality segment.

	Qualities are assumed to be ASCII-encoded as chr(qual + base).

	The algorithm is the same as the one used by BWA within the function
	'bwa_trim_read':
	- Subtract the cutoff value from all qualities.
	- Compute partial sums from all indices to the end of the sequence.
	- Trim sequence at the index at which the sum is minimal.
	"""
	cdef int s
	cdef int max_qual
	cdef int stop = len(qualities)
	cdef int start = 0
	cdef int i

	# find trim position for 5' end
	s = 0
	max_qual = 0
	for i in range(len(qualities)):
		s += cutoff_front - (ord(qualities[i]) - base)
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			start = i + 1

	# same for 3' end
	max_qual = 0
	s = 0
	for i in reversed(xrange(len(qualities))):
		s += cutoff_back - (ord(qualities[i]) - base)
		if s < 0:
			break
		if s > max_qual:
			max_qual = s
			stop = i
	if start >= stop:
		start, stop = 0, 0
	return (start, stop)
