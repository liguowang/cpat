#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 18:47:37 2020

"""
import sys
import os
import signal
import subprocess
import pysam
import logging
from cpmodule import ireader
from cpmodule import find_orfs
from cpmodule import fickett
from cpmodule import FrameKmer


def make_logit(infile, rscript, outfile):
    '''make logit model'''
    logging.info(
        "Making logistic regression model from \"%s\" ..." % infile)
    RCMD = open(rscript, 'w')
    print(
        'data <- read.table(file=\"%s\",sep="\\t",header=T)' % (infile),
        file=RCMD)
    print('attach(data)', file=RCMD)
    print(
        'mylogit <- glm(%s ~ %s + %s + %s + %s, family=binomial(link=\"logit\"), na.action=na.pass)'
        % ("Label", "mRNA", "ORF", "Fickett", "Hexamer"), file=RCMD)
    print(
        'save.image(\"%s\")' % (outfile), file=RCMD)
    RCMD.close()
    try:
        subprocess.call("Rscript " + rscript, shell=True)
    except:
        logging.info("Cannot run R script: %s ...", rscript)
        sys.exit()
    logging.info("The logistic regression model is saved as %s." % outfile)


def signal_handler(sig, frame):
    logging.info('\nYou pressed Ctrl+C. Exit!')
    sys.exit(0)


signal.signal(signal.SIGINT, signal_handler)


def reverse_comp(s):
    swap = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "X": "X"}
    return "".join(swap[b] for b in s)[::-1]


def coding_prediction(rdata, idata, outfile):
    '''
    rdata stored the linear regression model, idata is data matrix
    containing features.
    '''
    RCMD = open(outfile + '.r', 'w')
    print('load(\"%s\")' % (rdata), file=RCMD)
    print(
        'test <- read.table(file=\"%s\",sep="\\t",header=T)'
        % (idata), file=RCMD)
    print(
        'test$Coding_prob <- predict(mylogit,newdata=test,type="response")',
        file=RCMD)
    print(
        'write.table(test, file=\"%s\", quote=F, sep="\\t",row.names=FALSE, col.names=TRUE)'
        % (outfile + '.ORF_prob.tsv'), file=RCMD)
    RCMD.close()
    try:
        subprocess.call("Rscript " + outfile + '.r', shell=True)
    except:
        pass
    logging.info("Removing file %s", idata)
    os.remove(idata)


def bed_or_fasta(infile):
    '''
    Determine if the input file is bed or fasta format.
    '''
    format = "UNKNOWN"
    for line in ireader.reader(infile):
        # line = line.strip()
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            format = "FASTA"
            return format
        elif len(line.split()) >= 12:
            format = 'BED'
            return format
    return format


def index_fasta(infile):
    '''
    Index fasta file using samTools.
    '''
    try:
        if os.path.getsize(infile + '.fai'):
            logging.debug("\"%s\" exists. Skip indexing!", (infile + '.fai'))
            pass
    except OSError:
        logging.warning(
            "Can not find the index file: \"%s\"", (infile + '.fai'))
        logging.info(
            "Indexing \"%s\" using the \"pysam\" module...", infile)
        pysam.faidx(infile)
        logging.info("Done!")


def seq_from_bed(inbed, refgenome):
    '''
    Extract seq of sequence from bed line.
    '''

    transtab = str.maketrans("ACGTNX", "TGCANX")
    mRNA_seq = ''
    if inbed.strip():
        try:
            fields = inbed.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            # tx_end = int(fields[2])
            geneName = fields[3]
            strand = fields[5].replace(" ", "_")
            # exon_num = int(fields[9])
            # exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
            exon_starts = list(map(int, fields[11].rstrip(',\n').split(',')))
            exon_starts = list(map((lambda x: x + tx_start), exon_starts))
            exon_ends = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
            # intron_starts = exon_ends[:-1]
            # intron_ends = exon_starts[1:]
        except:
            logging.error("Wrong format!" + inbed)
            return None
        # mRNA_size = sum(exon_sizes)
        for st, end in zip(exon_starts, exon_ends):
            exon_coord = chrom + ':' + str(st + 1) + '-' + str(end)
            tmp = pysam.faidx(refgenome, exon_coord)
            mRNA_seq += ''.join([i.rstrip('\n\r') for i in tmp.split()[1:]])
        if strand == '-':
            mRNA_seq = mRNA_seq.upper().translate(transtab)[::-1]
        return(geneName, mRNA_seq)


def finish_up(prefix, n_top, min_len):
    print("  Output files:", file=sys.stderr)
    print("\t*%s.ORF_seqs.fa: The top %d ORF sequences (at least %d nucleotides long) in FASTA format."
          % (prefix, n_top, min_len), file=sys.stderr)
    print("\t*%s.ORF_prob.tsv: ORF information (strand, frame, start, end, size, Fickett TESTCODE score, Hexamer score) and coding probability)"
          % prefix, file=sys.stderr)
    print("\t*%s.ORF_prob.best.tsv: The information of the best ORF. This file is a subset of \"%s.ORF_prob.tsv\""
          % (prefix, prefix), file=sys.stderr)
    print("\t*%s.no_ORF.txt: Sequence IDs or BED entried with no ORF found."
          % prefix, file=sys.stderr)
    print("\t*%s.r: Rscript file."
          % prefix, file=sys.stderr)


def extract_feature_from_bed(inbed, refgenome, stt, stp, c_tab, g_tab, min_orf):
    '''extract features of sequence from bed line'''
    transtab = str.maketrans("ACGTNX", "TGCANX")
    mRNA_seq = ''
    mRNA_size = 0
    if inbed.strip():
        try:
            fields = inbed.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            # tx_end = int(fields[2])
            geneName = fields[3]
            strand = fields[5].replace(" ", "_")
            exon_sizes = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_starts = list(map(int, fields[11].rstrip(',\n').split(',')))
            exon_starts = list(map((lambda x: x + tx_start), exon_starts))
            exon_ends = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
        except:
            print("Wrong format!" + inbed, file=sys.stderr)
            return None
        mRNA_size = sum(exon_sizes)
        for st, end in zip(exon_starts, exon_ends):
            exon_coord = chrom + ':' + str(st + 1) + '-' + str(end)
            tmp1 = pysam.faidx(refgenome, exon_coord)
            mRNA_seq += ''.join([i.rstrip('\n\r') for i in tmp1.split()[1:]])
        if strand == '-':
            mRNA_seq = mRNA_seq.upper().translate(transtab)[::-1]
        tmp2 = find_orfs.ORFFinder(mRNA_seq, min_orf=min_orf)
        ORFs = tmp2.orf_candidates(start_coden=stt, stop_coden=stp,
                                   antisense=False, n_candidate=1)
        if len(ORFs) == 0:
            return None
        (direction, frame, ORF_start, ORF_end, CDS_size, CDS_seq) = ORFs[0]
        # print (ORFs)
        fickett_score = fickett.fickett_value(CDS_seq)
        hexamer = FrameKmer.kmer_ratio(CDS_seq, 6, 3, c_tab, g_tab)
        return (geneName, mRNA_size, CDS_size, fickett_score, hexamer)


def extract_CDS_from_bed(inbed, refgenome, stt, stp, c_tab, g_tab, min_orf):
    '''extract CDS sequence from bed line'''
    transtab = str.maketrans("ACGTNX", "TGCANX")
    CDS_seq = ''
    mRNA_size = 0
    if inbed.strip():
        try:
            fields = inbed.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            # tx_end = int( fields[2] )
            geneName = fields[3]
            strand = fields[5].replace(" ", "_")
            cdsStart = int(fields[6])
            cdsEnd = int(fields[7])
            exon_sizes = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_starts = list(map(int, fields[11].rstrip(',\n').split(',')))
            exon_starts = list(map((lambda x: x + tx_start), exon_starts))
            exon_ends = list(map(int, fields[10].rstrip(',\n').split(',')))
            exon_ends = list(map((lambda x, y: x + y), exon_starts, exon_ends))
        except:
            print("Wrong format!" + inbed, file=sys.stderr)
            return None
        mRNA_size = sum(exon_sizes)

        for base, offset in zip(exon_starts, exon_sizes):
            if (base + offset) < cdsStart:
                continue
            if base > cdsEnd:
                continue
            cds_exon_start = max(base, cdsStart)
            cds_exon_end = min(base + offset, cdsEnd)
            exon_coord = chrom + ':' + str(cds_exon_start + 1) + '-' + str(cds_exon_end)
            tmp1 = pysam.faidx(refgenome, exon_coord)
            CDS_seq += ''.join([i.rstrip('\n\r') for i in tmp1.split()[1:]])
        if strand == '-':
            CDS_seq = CDS_seq.upper().translate(transtab)[::-1]
        CDS_size = len(CDS_seq)
        fickett_score = fickett.fickett_value(CDS_seq)
        hexamer = FrameKmer.kmer_ratio(CDS_seq, 6, 3, c_tab, g_tab)
        return (geneName, mRNA_size, CDS_size, fickett_score, hexamer)


def extract_feature_from_seq(seq, stt, stp, c_tab, g_tab, min_orf):
    '''extract features of sequence from fasta entry'''

    mRNA_seq = seq.upper()
    mRNA_size = len(seq)
    tmp = find_orfs.ORFFinder(mRNA_seq, min_orf=min_orf)
    ORFs = tmp.orf_candidates(
        start_coden=stt, stop_coden=stp, antisense=False, n_candidate=3)
    if len(ORFs) == 0:
        return None
    (direction, frame, ORF_start, ORF_end, CDS_size, CDS_seq) = ORFs[0]
    fickett_score1 = fickett.fickett_value(CDS_seq)
    hexamer = FrameKmer.kmer_ratio(CDS_seq, 6, 3, c_tab, g_tab)
    return (mRNA_size, CDS_size, fickett_score1, hexamer)
