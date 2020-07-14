#!/usr/bin/env python

'''
--------------------------------------
Find the longest Open Reading Frame
--------------------------------------
'''

import sys,os
import collections
from textwrap import wrap
from optparse import OptionParser
from cpmodule  import orf
from cpmodule  import fasta

__author__ = "Liguo Wang"
__contributor__="Liguo Wang"
__copyright__ = "Copyright 2020, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="2.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"

def reverse_comp(s):
	swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
	return "".join(swap[b] for b in s)[::-1]

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-s","--sequence",action="store",type="string", dest="sequence_file",help="DNA sequnence(s) in FASTA format.")
	parser.add_option("-o","--outfile",action="store",type="string", dest="out_file",help="Prefix of output files. Sequences in \"prefix.longest_orf.fa\" are always in 5' -> 3' direction. ORF nucleotides are in uppercase, everything else is in lowercase.")
	parser.add_option("--start",action="store",type="string", dest="start_codons",default='ATG',help="Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). default=%default")
	parser.add_option("--stop",action="store",type="string", dest="stop_codons",default='TAG,TAA,TGA',help="Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ','. default=%default")
	parser.add_option("--width",action="store",type="int", dest="line_width",default=100,help="Line width of output fasta file.  default=%default")
	(options,args)=parser.parse_args()

	for file in ([options.sequence_file,options.out_file]):
		if not (file):
			parser.print_help()
			sys.exit(0)
	
	
	start_codons = options.start_codons.replace(' ','').split(',')
	stop_codons = options.stop_codons.replace(' ','').split(',')
	
	OUT = open(options.out_file + '.longest_orf.fa', 'w')
	
	print ("Start codons used: [%s]" % ','.join(start_codons), file=sys.stderr)
	print ("Stop codons used: [%s]" % ','.join(stop_codons), file=sys.stderr)

	
	print ("Read sequence file: \"%s\" ..." % options.sequence_file,file=sys.stderr)
	f = fasta.Fasta(fastafile = options.sequence_file)
	seqs = f.seqs
	print ("%d sequences loaded" % len(seqs))
	
	#results = collections.defaultdict(dict)
	longest_orfs = collections.defaultdict(list)
	
	print ("Searching for longest ORFs ...")
	#print ("\t".join(["Seq_ID", "fwd_Frame", "fwd_Start", "fwd_End", "fwd+Size","rev_Frame", "rev_Start", "rev_End", "rev+Size" ]))
	for name,seq in seqs.items():
		
		
		#look for orf on "+" strand
		tmp1 = orf.ORFFinder(seq)
		for frame in range(3):
			tmp1.run_one(frame, '+', start_codons, stop_codons)
		f_results = tmp1.result[1:]	#frame, start, end, length
		
		#look for orf on "-" strand
		rseq = reverse_comp(seq)
		tmp2 = orf.ORFFinder(rseq)
		for frame in range(3):
			tmp2.run_one(frame, '-', start_codons, stop_codons)
		r_results = tmp2.result[1:]
		#print (name + '\t' + '\t'.join([str(i) for i in f_results]) + '\t' + '\t'.join([str(i) for i in r_results]))
		
		#compare orf size between forward and reverse strand	
		if f_results[3] > r_results[3]:
			new_name = name + (" ORF_Strand='+',ORF_Frame=%d,ORF_st=%d,ORF_end=%d,ORF_size=%d" % (f_results[0], f_results[1], f_results[2],f_results[3]))
			p5 = seq[0:f_results[1]].lower()
			cds = seq[f_results[1]:f_results[2]].upper()
			p3 = seq[f_results[2]:].lower()
			print (">" + new_name, file=OUT)
			print ('\n'.join(wrap(p5 + cds + p3, width=options.line_width)), file=OUT)
		elif f_results[3] < r_results[3]:
			new_name = name + (" ORF_Strand='-',ORF_Frame=%d,ORF_st=%d,ORF_end=%d,ORF_size=%d" % (r_results[0], r_results[1], r_results[2], r_results[3]))
			p5 = rseq[0:r_results[1]].lower()
			cds = rseq[r_results[1]:r_results[2]].upper()
			p3 = rseq[r_results[2]:].lower()
			print (">" + new_name, file=OUT)
			print ('\n'.join(wrap(p5 + cds + p3, width=options.line_width)), file=OUT)
		

if __name__=='__main__':
	main()	










	