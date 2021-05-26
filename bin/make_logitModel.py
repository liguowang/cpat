#!/usr/bin/env python
'''
-----------------------------------------
build logit model from training datasets
-----------------------------------------
'''
import sys
if sys.version_info[0] < 3 or sys.version_info[1] < 5:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " Needs python3.5 or newer!\n", file=sys.stderr)
	sys.exit()

import logging
import pysam
from optparse import OptionParser
from cpmodule import fickett
from cpmodule  import find_orfs
from cpmodule  import FrameKmer
from cpmodule  import ireader
from cpmodule.utils import bed_or_fasta,index_fasta,signal_handler,make_logit

__author__ = "Liguo Wang"
__contributor__="Liguo Wang, Hyun Jung Park, Wei Li"
__copyright__ = "Copyright 2012, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="3.0.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu; wangliguo78@gmail.com"
__status__ = "Production"



def extract_feature_from_bed(inbed,refgenome,stt,stp,c_tab,g_tab,min_orf):
	'''extract features of sequence from bed line'''
	transtab = str.maketrans("ACGTNX","TGCANX")
	mRNA_seq = ''
	mRNA_size = 0
	if inbed.strip():
		try:
			fields = inbed.split()
			chrom = fields[0]
			tx_start = int( fields[1] )
			#tx_end = int( fields[2] )
			geneName = fields[3]
			strand = fields[5].replace(" ","_")
			exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
			exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
			exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
			exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
			exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));
		except:
			print("Wrong format!" + inbed, file=sys.stderr)
			return None
		mRNA_size = sum(exon_sizes)
		for st,end in zip(exon_starts, exon_ends):
			exon_coord = chrom + ':' + str(st +1) + '-' + str(end)
			tmp1 = pysam.faidx(refgenome,exon_coord)
			mRNA_seq += ''.join([i.rstrip('\n\r') for i in tmp1.split()[1:]])
		if strand =='-':
			mRNA_seq = mRNA_seq.upper().translate(transtab)[::-1]
		tmp2 = find_orfs.ORFFinder(mRNA_seq, min_orf=min_orf)
		ORFs = tmp2.orf_candidates(start_coden=stt, stop_coden=stp, antisense = False, n_candidate = 1)
		if len(ORFs) == 0:
			return None
		(direction, frame, ORF_start, ORF_end, CDS_size, CDS_seq) = ORFs[0]
		#print (ORFs)
		fickett_score = fickett.fickett_value(CDS_seq)
		hexamer = FrameKmer.kmer_ratio(CDS_seq,6,3,c_tab,g_tab)
		return (geneName, mRNA_size, CDS_size, fickett_score, hexamer)

def extract_CDS_from_bed(inbed,refgenome,stt,stp,c_tab,g_tab,min_orf):
	'''extract CDS sequence from bed line'''
	transtab = str.maketrans("ACGTNX","TGCANX")
	CDS_seq = ''
	mRNA_size = 0
	if inbed.strip():
		try:
			fields = inbed.split()
			chrom = fields[0]
			tx_start = int( fields[1] )
			#tx_end = int( fields[2] )
			geneName = fields[3]
			strand = fields[5].replace(" ","_")
			cdsStart = int(fields[6])
			cdsEnd = int(fields[7])
			exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
			exon_starts = list(map(int, fields[11].rstrip( ',\n' ).split( ',' ) ))
			exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
			exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
			exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));
		except:
			print("Wrong format!" + inbed, file=sys.stderr)
			return None
		mRNA_size = sum(exon_sizes)

		for base,offset in zip( exon_starts, exon_sizes):
			if (base + offset) < cdsStart: continue
			if base > cdsEnd: continue
			cds_exon_start = max(base, cdsStart )
			cds_exon_end = min(base+offset, cdsEnd )
			exon_coord = chrom + ':' + str(cds_exon_start +1) + '-' + str(cds_exon_end)
			tmp1 = pysam.faidx(refgenome,exon_coord)
			CDS_seq += ''.join([i.rstrip('\n\r') for i in tmp1.split()[1:]])
		if strand =='-':
			CDS_seq = CDS_seq.upper().translate(transtab)[::-1]
		CDS_size = len(CDS_seq)
		fickett_score = fickett.fickett_value(CDS_seq)
		hexamer = FrameKmer.kmer_ratio(CDS_seq,6,3,c_tab,g_tab)
		return (geneName, mRNA_size, CDS_size, fickett_score, hexamer)

def extract_feature_from_seq(seq,stt,stp,c_tab,g_tab,min_orf):
	'''extract features of sequence from fasta entry'''

	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = find_orfs.ORFFinder(mRNA_seq, min_orf = min_orf)
	ORFs = tmp.orf_candidates(start_coden=stt, stop_coden=stp, antisense = False, n_candidate = 3)
	if len(ORFs) == 0:
		return None
	(direction, frame, ORF_start, ORF_end, CDS_size, CDS_seq) = ORFs[0]
	fickett_score1 = fickett.fickett_value(CDS_seq)
	hexamer = FrameKmer.kmer_ratio(CDS_seq,6,3,c_tab,g_tab)
	return (mRNA_size, CDS_size, fickett_score1, hexamer)


def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-c","--cgene",action="store",dest="coding_file",help="Genomic sequnences of protein-coding RNAs in FASTA (https://en.wikipedia.org/wiki/FASTA_format) or standard 12-column BED (https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. It is recommended to use *short* and *unique* sequence identifiers (such as Ensembl transcript id) in FASTA and BED file. The input FASTA or BED file could be a regular text file or compressed file (*.gz, *.bz2) or accessible URL (http://, https://, ftp://). When BED file is provided, use the ORF defined in the BED file (the 7th and 8th columns in BED file define the positions of 'start codon, and 'stop codon', respectively). When FASTA file is provided, searching for the longet ORF. For well annotated genome, we recommend using BED file as input because the longest ORF predicted from RNA sequence might not be the real ORF. If this is a BED file, reference genome ('-r/--ref') should be specified.")
	parser.add_option("-n","--ngene",action="store",dest="noncoding_file",help="Genomic sequences of non-coding RNAs in FASTA (https://en.wikipedia.org/wiki/FASTA_format) or standard 12-column BED (https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. It is recommended to use *short* and *unique* sequence identifiers (such as Ensembl transcript id) in FASTA and BED file. The input FASTA or BED file could be a regular text file or compressed file (*.gz, *.bz2) or accessible URL (http://, https://, ftp://). If this is a BED file, reference genome ('-r/--ref') should be specified.")
	parser.add_option("-o","--outfile",action="store",dest="out_file",help="The prefix of output files.")
	parser.add_option("-x","--hex",action="store",dest="hexamer_dat",help="Hexamer frequency table. CPAT has prebuilt hexamer frequency tables for Human, Mouse, Fly, Zebrafish. Run 'make_hexamer_tab.py' to generate this table.")
	parser.add_option("-r","--ref",action="store",dest="ref_genome",help="Reference genome sequences in FASTA format. Ignore this option if mRNA sequences file was provided to '-g'. Reference genome file will be indexed automatically if the index file  *.fai) does not exist.")
	parser.add_option("-s","--start",action="store",dest="start_codons",default='ATG',help="Start codon (use 'T' instead of 'U') used to define the start of open reading frame (ORF). default=%default")
	parser.add_option("-t","--stop",action="store",dest="stop_codons",default='TAG,TAA,TGA',help="Stop codon (use 'T' instead of 'U') used to define the end of open reading frame (ORF). Multiple stop codons are separated by ','. default=%default")
	parser.add_option("--min-orf",action="store",type="int", dest="min_orf_len",default=30,help="Minimum ORF length in nucleotides.  default=%default")
	parser.add_option("--log-file",action="store",type="string", dest="log_file",default='make_logitModel_run_info.log',help="Name of log file. default=\"%default\"")
	parser.add_option("--verbose",action="store_true",dest="debug",default=False,help="Logical to determine if detailed running information is printed to screen.")
	(options,args)=parser.parse_args()

	#check input and output files
	for file in ([options.coding_file,options.noncoding_file,options.out_file,options.hexamer_dat]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	if options.debug:
		logging.basicConfig(filename='%s' % options.log_file, filemode='w',format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	else:
		logging.basicConfig(filename='%s' % options.log_file, filemode='w',format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
	#logging to console
	logFormat = logging.Formatter("%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S')
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormat)
	logging.getLogger().addHandler(consoleHandler)

	start_codons = options.start_codons.strip().split(',')
	stop_codons = options.stop_codons.strip().split(',')
	logging.info ("Start codons used: [%s]" % ','.join(start_codons))
	logging.info ("Stop codons used: [%s]" % ','.join(stop_codons))


	#data used to build logit model
	train_data=[]
	coding_label = 1
	noncoding_label = 0
	header = ['ID','mRNA','ORF','Fickett','Hexamer','Label']

	#build hexamer table from hexamer frequency file
	logging.info ("Reading hexamer frequency table file: \"%s\"" % options.hexamer_dat)
	coding={}
	noncoding={}
	for line in open(options.hexamer_dat, 'r'):
		line = line.strip()
		fields = line.split()
		if fields[0] == 'hexamer':continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] =  float(fields[2])

	######################################################################################
	# process protein coding transcripts
	######################################################################################
	count = 0
	logging.info("Process protein-coding RNA file: \"%s\"" % options.coding_file)


	file_format = bed_or_fasta(options.coding_file)
	if file_format == 'UNKNOWN':
		logging.error("Error: unknown file format \"%s\"" % options.coding_file)
		parser.print_help()
		sys.exit(0)
	elif file_format == 'BED':
		logging.info("Protein-coding RNA file \"%s\" is in BED format" % options.coding_file)
		if not options.ref_genome:
			logging.error("Error: Reference genome file must be provided")
			parser.print_help()
			sys.exit(0)
		index_fasta(options.ref_genome)

		for line in ireader.reader(options.coding_file):
			count +=1
			if count % 10 == 0:
				print("%d rows finished\r" % count, end=' ', file=sys.stderr)
			if line.startswith('track'):continue
			if line.startswith('#'):continue
			if line.startswith('browser'):continue

			#option-1: extract mRNA seq from BED and then find the longest ORF as CDS
			#features = extract_feature_from_bed(line, refgenome = options.ref_genome, stt = start_codons, stp = stop_codons, c_tab=coding, g_tab=noncoding, min_orf = options.min_orf_len)
			#if features is None:
			#	logging.warning("No ORF found for: %s" % '\t'.join(line.split()[0:6]))
			#	continue
			#(gene_id, mRNA_size, CDS_size, fickett_score, hexamer) = features

			#option-2: extract CDS directly from BED
			(gene_id, mRNA_size, CDS_size, fickett_score, hexamer) = extract_CDS_from_bed(line, refgenome = options.ref_genome, stt = start_codons, stp = stop_codons, c_tab=coding, g_tab=noncoding, min_orf = options.min_orf_len)

			train_data.append([gene_id, mRNA_size, CDS_size, fickett_score, hexamer, coding_label])
		logging.info("Total %d coding rows finished." % count)
	elif file_format == 'FASTA':
		if options.ref_genome:
			logging.warning("Reference genome sequence [-r] will be ignored when input file is in FASTA format.")
		logging.info("Protein-coding RNA file \"%s\" is in FASTA format" % options.coding_file)
		for sname,seq in FrameKmer.seq_generator(options.coding_file):
			count +=1
			if count % 10 == 0:
				print("%d sequences finished\r" % count, end=' ', file=sys.stderr)
			features = extract_feature_from_seq(seq = seq, stt = start_codons, stp = stop_codons, c_tab=coding, g_tab=noncoding, min_orf = options.min_orf_len)
			if features is None:
				continue
			(mRNA_size, CDS_size, fickett_score,hexamer) = features
			train_data.append([sname, mRNA_size, CDS_size, fickett_score, hexamer, coding_label])
		logging.info("Total %d coding sequences finished." % count)

	######################################################################################
	# process Non-protein coding transcripts
	######################################################################################
	count=0
	logging.info("Process non-coding RNA file: \"%s\"" % options.noncoding_file)

	file_format = bed_or_fasta(options.noncoding_file)
	if file_format == 'UNKNOWN':
		logging.error("Error: unknown file format \"%s\"" % options.noncoding_file)
		parser.print_help()
		sys.exit(0)
	elif file_format == 'BED':
		logging.info("Non-coding RNA file \"%s\" is in BED format" % options.noncoding_file)
		if not options.ref_genome:
			logging.error("Error: Reference genome file must be provided")
			parser.print_help()
			sys.exit(0)
		index_fasta(options.ref_genome)

		for line in ireader.reader(options.noncoding_file):
			count +=1
			if count % 10 == 0:
				print("%d genes finished\r" % count, end=' ', file=sys.stderr)
			if line.startswith('track'):continue
			if line.startswith('#'):continue
			if line.startswith('browser'):continue
			fields = line.split()
			if int(fields[1]) != int(fields[6]):
				logging.warning("This seems to be protein-coding:%s" % '\t'.join(fields[0:6]) )

			#if not line.strip(): continue
			features = extract_feature_from_bed(line, refgenome = options.ref_genome, stt = start_codons, stp = stop_codons, c_tab=coding, g_tab=noncoding, min_orf = options.min_orf_len)
			if features is None:
				logging.warning("No ORF found for: %s" % '\t'.join(line.split()[0:6]))
				continue
			(gene_id, mRNA_size, CDS_size, fickett_score,hexamer) = features
			train_data.append([gene_id, mRNA_size, CDS_size, fickett_score,hexamer,noncoding_label])
		logging.info("Total %d non-coding rows finished." % count)
	elif file_format == 'FASTA':
		if options.ref_genome:
			logging.warning("Reference genome sequence [-r] will be ignored when input file is in FASTA format.")
		logging.info("Non-coding RNA file \"%s\" is in FASTA format" % options.noncoding_file)
		for sname,seq in FrameKmer.seq_generator(options.noncoding_file):
			count +=1
			if count % 10 == 0:
				print("%d sequences finished\r" % count, end=' ', file=sys.stderr)
			#geneSeq = fa.getSeq(seqID = geneID)
			features = extract_feature_from_seq(seq = seq, stt = start_codons, stp = stop_codons, c_tab=coding, g_tab=noncoding, min_orf = options.min_orf_len)
			if features is None:
				continue
			(mRNA_size, CDS_size, fickett_score,hexamer) = features
			train_data.append([sname, mRNA_size, CDS_size, fickett_score,hexamer,noncoding_label])
		logging.info("Total %d non-coding sequences finished." % count)
	######################################################################################
	# writing data
	######################################################################################
	logging.info("Wrting to \"%s\"" % (options.out_file + '.feature.xls'))
	TMP = open(options.out_file + '.feature.xls', 'w')
	print('\t'.join(header), file=TMP)
	for i in train_data:
		print('\t'.join([str(j) for j in i]), file=TMP)
	TMP.close()

	#print("build logi model ...", file=sys.stderr)
	make_logit(options.out_file + '.feature.xls',options.out_file + '.make_logitModel.r', options.out_file + '.logit.RData')


if __name__ == '__main__':
	main()
