#!/usr/bin/env python

'''
--------------------------------------
Find the longest Open Reading Frame
--------------------------------------
'''
import sys
import logging
import signal
from textwrap import wrap
from optparse import OptionParser

from cpmodule import fickett
from cpmodule  import FrameKmer
from cpmodule  import ireader
from cpmodule import find_orfs
from cpmodule.utils import *

__author__ = "Liguo Wang"
__contributor__="Liguo Wang"
__copyright__ = "Copyright 2020, Mayo Clinic"
__credits__ = []
__license__ = "GPLv2"
__version__="3.0.1"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"		


def signal_handler(sig, frame):
	logging.info('\nYou pressed Ctrl+C. Exit!')
	sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

			
def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-g","--gene",action="store",type="string", dest="gene_file",help="Genomic sequnence(s) of RNA in FASTA or BED format. If this is a BED file, '-r/--ref' must also be specified. It is recommended to use *short* and *unique* sequence identifiers in this FASTA or BED file. The input  FASTA or BED file could be a regular text file or compressed file (*.gz, *.bz2) or accessible URL (http://, https://, ftp://).")
	parser.add_option("-o","--outfile",action="store",type="string", dest="out_file",help="The prefix of output files.")
	parser.add_option("-d","--logitModel",action="store",dest="logit_model",help="Prebuilt training model (Human, Mouse, Fly, Zebrafish). Run 'make_logitModel.py' to build logit model out of your own training datset.")
	parser.add_option("-x","--hex",action="store",dest="hexamer_dat",help="Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish). Run 'make_hexamer_tab.py' to make this table out of your own training dataset.")
	parser.add_option("-r","--ref",action="store",dest="ref_genome",help="Reference genome sequences in FASTA format. Ignore this option if FASTA file was provided to '-g/--gene'. Reference genome file will be indexed automatically (produce *.fai file along with the original *.fa file within the same directory) if hasn't been done.")
	parser.add_option("--antisense",action="store_true",dest="antisense",default=False,help="Also search for ORFs from the anti-sense strand. *Sense strand* (or coding strand) is DNA strand that carries the translatable code in the 5′ to 3′ direction. default=False (i.e. only search for ORFs from the sense strand)")
	parser.add_option("--start",action="store",type="string", dest="start_codons",default='ATG',help="Start codon used by ORFs. Use 'T' instead of 'U'. default=%default")
	parser.add_option("--stop",action="store",type="string", dest="stop_codons",default='TAG,TAA,TGA',help="Stop codons used by ORFs. Multiple stop codons should be separated by ','. Use 'T' instead of 'U'. default=%default")
	parser.add_option("--min-orf",action="store",type="int", dest="min_orf_len",default=75,help="Minimum ORF length in nucleotides.  default=%default")
	parser.add_option("--top-orf",action="store",type="int", dest="n_top_orf",default=5,help="Number of top ORF candidates. Many RNAs have dozens of possible ORFs, in most cases, the real ORF is ranked (by size) in the top several. To increase speed, we do not need to calculate \"Fickett\" score, \"Hexamer\" score and \"coding probability\" for all of them. default=%default")
	parser.add_option("--width",action="store",type="int", dest="line_width",default=100,help="Line width of output ORFs in FASTA format.  default=%default")
	parser.add_option("--log-file",action="store",type="string", dest="log_file",default='CPAT_run_info.log',help="Name of log file. default=\"%default\"")
	parser.add_option("--best-orf",action="store",type="string", dest="mode",default='p',help="Criteria to select the best ORF: \"l\"=length, selection according to the \"ORF length\"; \"p\"=probability, selection according to the \"coding probability\". The \"p\" mode usually gives more accurate prediction than the \"l\"mode. default=\"%default\"")
	(options,args)=parser.parse_args()
	
	for file in ([options.gene_file, options.hexamer_dat, options.logit_model, options.out_file]):
		if not (file):
			parser.print_help()
			sys.exit(0)
	
	if options.line_width < 1:
		sys.exit(0)
	
	if options.mode not in ["p", "l"]:
		print ("Please specifiy either \"p\" or \"l\" to --best-orf.", file=sys.stderr)
		sys.exit(0)	
	
	#logging to file
	logging.basicConfig(filename='%s' % options.log_file, filemode='w',format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	#logging to console
	logFormat = logging.Formatter("%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S')
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormat)
	logging.getLogger().addHandler(consoleHandler)
				  
		
	logging.info ("Running CPAT version %s..." % ( __version__))
	start_codons = options.start_codons.replace(' ','').split(',')
	stop_codons = options.stop_codons.replace(' ','').split(',')
	
	
	SEQOUT = open(options.out_file + '.ORF_seqs.fa', 'w')
	INFOUT = open(options.out_file + '.ORF_info.tsv', 'w')
	NOORF = open(options.out_file + '.no_ORF.txt', 'w')
	
	logging.info ("Start codons used: [%s]" % ','.join(start_codons))
	logging.info ("Stop codons used: [%s]" % ','.join(stop_codons))
	
	#build hexamer table from hexamer frequency file
	logging.info ("Reading %s" % options.hexamer_dat)
	coding = {}
	noncoding = {}
	for line in open(options.hexamer_dat):
		line = line.strip()
		fields = line.split()
		if fields[0] == 'hexamer':
			continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] =  float(fields[2])
                	
	
	count = 0
	logging.info ("Checking format of \"%s\"" % options.gene_file)
	file_format = bed_or_fasta(options.gene_file)
	if file_format == 'UNKNOWN':
		logging.error("Unknown file format:%s" % options.gene_file)
		sys.exit(0)
	
	elif file_format == 'FASTA':
		logging.info("Input gene file is in FASTA format")
		if options.ref_genome:
			logging.warning("\"%s\" is a sequence file. The reference genome file \"%s\" will be ignored." % (options.gene_file, options.ref_genome))

		logging.info ("Searching for ORFs ...")
		print ("\t".join(["ID", "mRNA","ORF_strand", "ORF_frame", "ORF_start", "ORF_end","ORF", "Fickett", "Hexamer" ]), file=INFOUT)	## do NOT change these labels, they are R variable names in the model. 
		for name,seq in FrameKmer.seq_generator(options.gene_file):
			count += 1
			RNA_len = len(seq)
			#ORF serial number, starting from 1
			orf_sn = 1
			tmp1 = find_orfs.ORFFinder(seq = seq, min_orf = options.min_orf_len)
			ORFs = tmp1.orf_candidates(antisense = options.antisense, n_candidate = options.n_top_orf)
			if len(ORFs) == 0:
				logging.warning("No ORFs found for %s" % name)
				print (name, file=NOORF)
				continue
			for orf in ORFs:
				# (direction, frame_number+1, orf_start, orf_end, L, sequence)
				orf_seq = orf[-1]
				if orf[0] == '+':
					orf[2] = orf[2] + 1	#change 0-based into 1-based to be consistent with NCBI ORFfinder output (https://www.ncbi.nlm.nih.gov/orffinder/)
				elif orf[0] == '-':
					orf[2] = RNA_len - (orf[2])
					orf[3] = RNA_len - orf[3] + 1
			
				orf_id = name + '_ORF_' + str(orf_sn) + '\t' + str(RNA_len) + '\t' + '\t'.join([str(i) for i in orf[:-1]])
			
				fickett_score = fickett.fickett_value(orf_seq)
				hexamer_score = FrameKmer.kmer_ratio(orf_seq,6,3,coding,noncoding)
				print (orf_id + '\t' + str(fickett_score) + '\t' + str(hexamer_score), file=INFOUT)
			
				print (">" + orf_id, file=SEQOUT)
				print ('\n'.join(wrap(orf_seq, width = options.line_width)), file=SEQOUT)
				orf_sn += 1
			print("%d sequences finished\r" % count, end=' ', file=sys.stderr)
		print ("\n", file=sys.stderr)
	
	
	elif file_format == 'BED':
		logging.info("Input gene file is in BED format")
		if not options.ref_genome:
			logging.error("Reference genome file (-r/--ref) must be provided.")
			parser.print_help()
			sys.exit(0)	
		
		logging.info ("Searching for ORFs ...")
		print ("\t".join(["ID", "mRNA","ORF_strand", "ORF_frame", "ORF_start", "ORF_end","ORF", "Fickett", "Hexamer" ]), file=INFOUT)	## do NOT change these labels, they are R variable names in the model. 
		
		index_fasta(options.ref_genome)	
		
		for line in ireader.reader(options.gene_file):
			count +=1
			if line.startswith('track'):continue
			if line.startswith('#'):continue
			if line.startswith('browser'):continue	
			name,seq = seq_from_bed(line, options.ref_genome)
			
			
			RNA_len = len(seq)
			#ORF serial number, starting from 1
			orf_sn = 1
			tmp1 = find_orfs.ORFFinder(seq = seq, min_orf = options.min_orf_len)
			ORFs = tmp1.orf_candidates(antisense = options.antisense, n_candidate = options.n_top_orf)
			if len(ORFs) == 0:
				logging.warning("No ORFs found for %s" % name)
				print (line, file=NOORF)
				continue
			for orf in ORFs:
				# (direction, frame_number+1, orf_start, orf_end, L, sequence)
				orf_seq = orf[-1]
				if orf[0] == '+':
					orf[2] = orf[2] + 1	#change 0-based into 1-based to be consistent with NCBI ORFfinder output (https://www.ncbi.nlm.nih.gov/orffinder/)
				elif orf[0] == '-':
					orf[2] = RNA_len - (orf[2])
					orf[3] = RNA_len - orf[3] + 1
			
				orf_id = name + '_ORF_' + str(orf_sn) + '\t' + str(RNA_len) + '\t' + '\t'.join([str(i) for i in orf[:-1]])
			
				fickett_score = fickett.fickett_value(orf_seq)
				hexamer_score = FrameKmer.kmer_ratio(orf_seq,6,3,coding,noncoding)
				print (orf_id + '\t' + str(fickett_score) + '\t' + str(hexamer_score), file=INFOUT)
			
				print (">" + orf_id, file=SEQOUT)
				print ('\n'.join(wrap(orf_seq, width = options.line_width)), file=SEQOUT)
				orf_sn += 1
			print("%d rows finished\r" % count, end=' ', file=sys.stderr)
		print ("\n", file=sys.stderr)

	
	SEQOUT.close()
	INFOUT.close()
	
	logging.info ("Calculate coding probability ...")
	coding_prediction(options.logit_model, options.out_file + '.ORF_info.tsv', options.out_file)	#output options.out_file + '.ORF_prob.tsv'
	
	if options.mode == 'p':
		logging.info ("Select ORF with the highest coding probability ...")
		col_index = 9
	elif options.mode == 'l':
		logging.info ("Select the longest ORF ...")
		col_index = 6
		
	BEST = open((options.out_file + '.ORF_prob.best.tsv'), 'w')
	best_candidates = {}
	for l in open((options.out_file + '.ORF_prob.tsv'), 'r'):
		l = l.strip()
		if l.startswith('ID'):
			print ("seq_ID\t" + l, file=BEST)
			continue
		f = l.split('\t')
		seq_id = f[0].split('_ORF_')[0]
		prob = float(f[col_index])
		if seq_id not in best_candidates:
			best_candidates[seq_id] = f
		else:
			if prob > float(best_candidates[seq_id][col_index]):
				best_candidates[seq_id] = f
	
	for k,v in best_candidates.items():
		print (k + '\t' + '\t'.join(v), file=BEST)
	
	BEST.close()		
	logging.info ("Done!")
	
	finish_up(options.out_file, options.n_top_orf, options.min_orf_len)
	
		
if __name__=='__main__':
	main()	
