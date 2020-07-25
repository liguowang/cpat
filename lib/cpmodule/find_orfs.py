import sys
class ORFFinder:
	"""Find the longest ORF in a given sequence 
	"seq" is a string, if "start" is not provided any codon can be the start of 
	and ORF. If muliple ORFs have the longest length the first one encountered
	is printed 
	"""
	def __init__(self, seq, min_orf = 75):
		self.seq = seq.upper()
		self.result = []
		#self.winner = 0
		self.cutoff = min_orf	#minimum ORF size. 75 is the default value for NCBI ORF finder
  
	def _reverse_comp(self):
		swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
		return "".join(swap[b] for b in self.seq)[::-1]
  
	def codons(self, frame):
		""" A generator that yields DNA in one codon blocks 
		"frame" counts for 0. This function yelids a tuple (triplet, index) with 
		index relative to the original DNA sequence 
		"""
		start = frame
		while start + 3 <= len(self.seq):
			yield (self.seq[start:start+3], start)
			start += 3 

	def run_one(self, frame_number, direction,start_coden, stop_coden):
		""" Search in one reading frame """
		codon_gen = self.codons(frame_number)  
		start_codens = start_coden
		stop_codens = stop_coden   
		while True:
			try: 
				c , index = next(codon_gen)
			except StopIteration:
				break 
			# Lots of conditions here: checks if we care about looking for start 
			# codon then that codon is not a stop
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
    
		#search for '+' strand
		for frame in range(3):
			self.run_one(frame, '+', start_coden, stop_coden)
    
		if antisense:
			# Also search ORFs from antisense strand
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
  