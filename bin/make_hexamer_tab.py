#!/usr/bin/env python
'''---------------------------------------------------------------------------------------
Make Hexamer frequency table
------------------------------------------------------------------------------------------'''

import os,sys
import string
from optparse import OptionParser
import warnings
import string
from cpmodule.FrameKmer import kmer_freq_file
__author__ = "Liguo Wang"
__contributor__="Liguo Wang, Hyun Jung Park, Wei Li"
__copyright__ = "Copyright 2012, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="2.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu;wangliguo78@gmail.com"
__status__ = "Production"

def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-c","--cod",action="store",dest="coding_file",help="Coding sequence (must be CDS without UTR, i.e. from start coden to stop coden) in fasta format")
	parser.add_option("-n","--noncod",action="store",dest="noncoding_file",help="Noncoding sequences in fasta format")
	(options,args)=parser.parse_args()

	if not options.coding_file and not options.noncoding_file:
		parser.print_help()
		sys.exit(0)		
	cod = kmer_freq_file(fastafile = options.coding_file, word_size = 6, step_size = 3, frame = 0)
	noncod = kmer_freq_file(fastafile = options.noncoding_file, word_size = 6, step_size = 1, frame = 0)
	
	#for i,j in cod.items():
	#	print str(i) + '\t' + str(j)
	
	cod_sum = 0.0
	cod_sum += sum(cod.values())
	noncod_sum = 0.0
	noncod_sum += sum(noncod.values())
	
	print('hexamer' + '\t' + 'coding' + '\t' + 'noncoding')
	for kmer in cod:
		if 'N' in kmer:
			continue
		print(kmer + '\t' + str(float(cod[kmer]/cod_sum))  + '\t' + str(float(noncod[kmer]/noncod_sum))) 

if __name__ == '__main__':
	main()
