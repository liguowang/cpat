import sys
class ORFFinder:
	"""
	Find the ORFs in a given sequence. 
	"""
	def __init__(self, seq, min_orf = 75):
		self.seq = seq.upper()
		self.result = []
		#self.winner = 0
		self.cutoff = min_orf	#minimum ORF size. 75 is the default value for NCBI ORF finder
  
	def _reverse_comp(self):
		"""
		Reverse complement DNA nucleotides.
		"""
		swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
		return "".join(swap[b] for b in self.seq)[::-1]
  
	def codons(self, frame):
		"""
		A generator that yields DNA in one codon blocks.
		"""
		start = frame
		while start + 3 <= len(self.seq):
			yield (self.seq[start:start+3], start)
			start += 3 

	def run_one(self, frame_number, direction, start_coden, stop_coden):
		"""
		Search in one reading frame
		
		Parameters
		----------
		frame_number : int
			Reading frame. Must be 0, 1, or 2.
		direction : str
			Sense ('+') or antisense ('-'). Must be '+' or '-'. 
		start_coden : list of string
			Start coden(s).
		stop_coden : list of string
			Stop coden(s).
		
		Return
		------
			List of ORF candidates ([[direction, frame, ORF_start, ORF_end, ORF_length, ORF_seq],...])
		"""
		codon_gen = self.codons(frame_number)  
		start_codens = start_coden
		stop_codens = stop_coden   
		while True:
			try: 
				c , index = next(codon_gen)
			except StopIteration:
				break 
			if c in start_codens or not start_codens and c not in stop_codens:
				orf_start = index  # we'll return the result as 0-indexed
				end = False
				while True:
					try: 
						c, index = next(codon_gen)
					except StopIteration:
						end = True
					if c in stop_codens:
						end = True
					if end:
						orf_end = index + 3 # because index is realitve to start of codon
						L = (orf_end - orf_start)
						if L > self.cutoff:
							self.result.append( [direction, frame_number+1, orf_start, orf_end, L, self.seq[orf_start:orf_end]])
						break
    
	def orf_candidates(self, start_coden=['ATG'], stop_coden=['TAG','TAA','TGA'], antisense = False, n_candidate = 3):
		"""
		Find ORF candidates.
		
		Parameters
		----------
		start_coden : list of string
			Start coden(s).
		stop_coden : list of string
			Stop coden(s).
		antisense : bool
			Whether to search ORFs from the antisense strand.
		n_candidate : int
			Number of candidate ORFs returned. 
		
		Return
		------
		List of list
			List of ORF candidates ([[direction, frame, ORF_start, ORF_end, ORF_length, ORF_seq],...])
		"""
    
		#search for ORFs from the '+' strand
		for frame in range(3):
			self.run_one(frame, '+', start_coden, stop_coden)
    	# Also search ORFs from the antisense strand
		if antisense:
			try:
				self.seq = self._reverse_comp()
			except:
				print (self.seq)
				sys.exit(0)	  
        
			for frame in range(3):
				self.run_one(frame, '-', start_coden, stop_coden)    
		candidates = self.result
		return sorted(candidates, key=lambda x: x[4], reverse=True)[:n_candidate]
      
#===================
def little_test():
  seq=''
  for line in open(sys.argv[1],'r'):
    line=line.strip('\n\r')
    if line.startswith('>'):
  	  continue
    seq	+= line
  tmp = ORFFinder(seq).orf_candidates()
  for orf in tmp:
  	print ('\t'.join([str(i) for i in orf]))
  
if __name__ == "__main__":
  little_test()
  