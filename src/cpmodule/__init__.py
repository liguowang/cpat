import sys
import logging
from textwrap import wrap
from optparse import OptionParser

from cpmodule import fickett
from cpmodule import FrameKmer
from cpmodule import ireader
from cpmodule import find_orfs

from cpmodule.utils import index_fasta, make_logit
from cpmodule.utils import extract_feature_from_seq
from cpmodule.utils import extract_feature_from_bed
from cpmodule.utils import extract_CDS_from_bed
from cpmodule.utils import seq_from_bed
from cpmodule.utils import bed_or_fasta
from cpmodule.utils import coding_prediction
from cpmodule.utils import finish_up
from cpmodule.FrameKmer import kmer_freq_file

__author__ = "Liguo Wang"
__contributor__ = "Liguo Wang"
__copyright__ = "Copyright 2024, Mayo Clinic"
__credits__ = []
__license__ = "GPLv2"
__version__ = "3.0.5"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"


def cpat():
    usage = "\n%prog  [options]"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-g", "--gene",
                      action="store", type="string",
                      dest="gene_file",
                      help="Genomic sequnence(s) of RNA in FASTA \
                      (https://en.wikipedia.org/wiki/FASTA_format) or \
                      standard 12-column  BED \
                      (https://genome.ucsc.edu/FAQ/FAQformat.html#format1) \
                       format. It is recommended to use *short* and *unique* \
                       sequence identifiers (such as Ensembl transcript id) \
                       in FASTA and BED file. If this is a BED file, \
                       reference genome ('-r/--ref') should be specified. \
                       The input FASTA or BED file could be a regular text \
                       file or compressed file (*.gz, *.bz2) or accessible \
                       URL (http://, https://, ftp://). URL file cannot be a \
                       compressed file.")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="out_file", help="The prefix of output files.")
    parser.add_option("-d", "--logitModel", action="store", dest="logit_model",
                      help="Logistic regression model. The prebuilt models \
                          for Human, Mouse, Fly, Zebrafish are availablel. \
                          Run 'make_logitModel.py' to build logistic \
                          regression model for your own training datset.")
    parser.add_option("-x", "--hex",
                      action="store",
                      dest="hexamer_dat", help="The hexamer frequency table. \
                      The prebuilt tables for Human, Mouse, Fly, Zebrafish \
                      are availablel. Run 'make_hexamer_tab.py' to make this \
                      table for your own training dataset.")
    parser.add_option("-r", "--ref",
                      action="store",
                      dest="ref_genome",
                      help="Reference genome sequences in FASTA format. \
                      Reference genome file will be indexed automatically \
                      if the index file ( *.fai) does not exist. Will be \
                      ignored if FASTA file was provided to '-g/--gene'.")
    parser.add_option("--antisense",
                      action="store_true",
                      dest="antisense",
                      default=False,
                      help="Logical to determine whether to search for ORFs \
                      from the anti-sense strand. *Sense strand* (or coding \
                      strand) is DNA strand that carries the translatable \
                      code in the 5′ to 3′ direction. default=False (i.e. \
                      only search for ORFs from the sense strand)")
    parser.add_option("--start",
                      action="store",
                      type="string",
                      dest="start_codons",
                      default='ATG',
                      help="Start codon (use 'T' instead of 'U') used to \
                      define the start of open reading frame (ORF). \
                      default=%default")
    parser.add_option("--stop",
                      action="store",
                      type="string",
                      dest="stop_codons",
                      default='TAG, TAA, TGA',
                      help="Stop codon (use 'T' instead of 'U') used to \
                      define the end of open reading frame (ORF). Multiple \
                      stop codons are separated by ','. default=%default")
    parser.add_option("--min-orf",
                      action="store",
                      type="int",
                      dest="min_orf_len",
                      default=75,
                      help="Minimum ORF length in nucleotides. \
                      default=%default")
    parser.add_option("--top-orf",
                      action="store",
                      type="int",
                      dest="n_top_orf",
                      default=5,
                      help="Number of ORF candidates reported. RNAs may \
                      have dozens of putative ORFs, in most cases, the real \
                      ORF is ranked (by size) in the top several. It is not \
                      necessary to calculate \"Fickett score\", \
                      \"Hexamer score\" and \"coding probability\" for every \
                      ORF. default=%default")
    parser.add_option("--width",
                      action="store",
                      type="int",
                      dest="line_width",
                      default=100,
                      help="Line width of output ORFs in FASTA format. \
                      default=%default")
    parser.add_option("--log-file",
                      action="store",
                      type="string",
                      dest="log_file",
                      default='CPAT_run_info.log',
                      help="Name of log file. default=\"%default\"")
    parser.add_option("--best-orf",
                      action="store",
                      type="string",
                      dest="mode",
                      default='p',
                      help="Criteria to select the best ORF: \"l\"=length, \
                      selection according to the \"ORF length\"; \
                      \"p\"=probability, selection according to the \
                      \"coding probability\". default=\"%default\"")
    parser.add_option("--verbose",
                      action="store_true",
                      dest="debug",
                      default=False,
                      help="Logical to determine if detailed running \
                      information is printed to screen.")
    (options, args) = parser.parse_args()

    for file in ([options.gene_file, options.hexamer_dat,
                  options.logit_model, options.out_file]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    if options.line_width < 1:
        sys.exit(0)

    if options.mode not in ["p", "l"]:
        print("Please specifiy either \"p\" or \"l\" to --best-orf.",
              file=sys.stderr)
        sys.exit(0)

    # logging to file
    if options.debug:
        logging.basicConfig(filename='%s' % options.log_file, filemode='w',
                            format="%(asctime)s [%(levelname)s] %(message)s",
                            datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
    else:
        logging.basicConfig(filename='%s' % options.log_file, filemode='w',
                            format="%(asctime)s [%(levelname)s]  %(message)s",
                            datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)

    # logging to console
    logFormat = logging.Formatter("%(asctime)s [%(levelname)s]  %(message)s",
                                  datefmt='%Y-%m-%d %I:%M:%S')
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormat)
    logging.getLogger().addHandler(consoleHandler)

    logging.info("Running CPAT version %s..." % (__version__))
    start_codons = options.start_codons.replace(' ', '').split(',')
    stop_codons = options.stop_codons.replace(' ', '').split(',')

    SEQOUT = open(options.out_file + '.ORF_seqs.fa', 'w')
    INFOUT = open(options.out_file + '.ORF_info.tsv', 'w')
    NOORF = open(options.out_file + '.no_ORF.txt', 'w')

    logging.info("Start codons used: [%s]" % ','.join(start_codons))
    logging.info("Stop codons used: [%s]" % ','.join(stop_codons))

    # build hexamer table from hexamer frequency file
    logging.info("Reading %s" % options.hexamer_dat)
    coding = {}
    noncoding = {}
    for line in open(options.hexamer_dat):
        line = line.strip()
        fields = line.split()
        if fields[0] == 'hexamer':
            continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])

    count = 0
    logging.info("Checking format of \"%s\"" % options.gene_file)
    file_format = bed_or_fasta(options.gene_file)
    if file_format == 'UNKNOWN':
        logging.error("Unknown file format:%s" % options.gene_file)
        sys.exit(0)

    elif file_format == 'FASTA':
        logging.info("Input gene file is in FASTA format")
        if options.ref_genome:
            logging.warning(
                "\"%s\" is a sequence file. The reference genome file \"%s\" \
                will be ignored." % (options.gene_file, options.ref_genome))

        logging.info("Searching for ORFs ...")
        # do NOT change these labels, they are R variable names in the model.
        print("\t".join(["ID", "mRNA", "ORF_strand", "ORF_frame", "ORF_start",
                         "ORF_end", "ORF", "Fickett", "Hexamer"]),
              file=INFOUT)
        for name, seq in FrameKmer.seq_generator(options.gene_file):
            count += 1
            RNA_len = len(seq)
            # ORF serial number, starting from 1
            orf_sn = 1
            tmp1 = find_orfs.ORFFinder(
                seq=seq, min_orf=options.min_orf_len)
            ORFs = tmp1.orf_candidates(antisense=options.antisense,
                                       n_candidate=options.n_top_orf,
                                       start_coden=start_codons,
                                       stop_coden=stop_codons)
            if len(ORFs) == 0:
                logging.warning("No ORFs found for %s" % name)
                print(name, file=NOORF)
                continue
            for orf in ORFs:
                # (direction, frame_number+1, orf_start, orf_end, L, sequence)
                orf_seq = orf[-1]
                if orf[0] == '+':
                    # change 0-based into 1-based to be consistent with NCBI
                    # ORFfinder output
                    # (https://www.ncbi.nlm.nih.gov/orffinder/)
                    orf[2] = orf[2] + 1
                elif orf[0] == '-':
                    orf[2] = RNA_len - (orf[2])
                    orf[3] = RNA_len - orf[3] + 1

                orf_id = name + '_ORF_' + str(orf_sn) + '\t' + str(RNA_len) + \
                    '\t' + '\t'.join([str(i) for i in orf[:-1]])

                fickett_score = fickett.fickett_value(orf_seq)
                hexamer_score = FrameKmer.kmer_ratio(
                    orf_seq, 6, 3, coding, noncoding)
                print(orf_id + '\t' + str(fickett_score) + '\t' +
                      str(hexamer_score), file=INFOUT)

                print(">" + orf_id, file=SEQOUT)
                print('\n'.join(wrap(orf_seq, width=options.line_width)),
                      file=SEQOUT)
                orf_sn += 1
            print("%d sequences finished\r" % count, end=' ', file=sys.stderr)
        print("\n", file=sys.stderr)

    elif file_format == 'BED':
        logging.info("Input gene file is in BED format")
        if not options.ref_genome:
            logging.error("Reference genome file (-r/--ref) must be provided.")
            parser.print_help()
            sys.exit(0)

        # do NOT change these labels, they are R variable names in the model.
        logging.info("Searching for ORFs ...")
        print("\t".join(["ID", "mRNA", "ORF_strand", "ORF_frame", "ORF_start",
                         "ORF_end", "ORF", "Fickett", "Hexamer"]), file=INFOUT)

        index_fasta(options.ref_genome)

        for line in ireader.reader(options.gene_file):
            count += 1
            if line.startswith('track'):
                continue
            if line.startswith('#'):
                continue
            if line.startswith('browser'):
                continue
            name, seq = seq_from_bed(line, options.ref_genome)

            RNA_len = len(seq)
            # ORF serial number, starting from 1
            orf_sn = 1
            tmp1 = find_orfs.ORFFinder(seq=seq,
                                       min_orf=options.min_orf_len)
            ORFs = tmp1.orf_candidates(antisense=options.antisense,
                                       n_candidate=options.n_top_orf,
                                       start_coden=start_codons,
                                       stop_coden=stop_codons)
            if len(ORFs) == 0:
                logging.warning("No ORFs found for %s" % name)
                print(line, file=NOORF)
                continue
            for orf in ORFs:
                # (direction, frame_number+1, orf_start, orf_end, L, sequence)
                orf_seq = orf[-1]
                if orf[0] == '+':
                    # change 0-based into 1-based to be consistent with NCBI
                    # ORFfinder output
                    # (https://www.ncbi.nlm.nih.gov/orffinder/)
                    orf[2] = orf[2] + 1
                elif orf[0] == '-':
                    orf[2] = RNA_len - (orf[2])
                    orf[3] = RNA_len - orf[3] + 1

                orf_id = name + '_ORF_' + str(orf_sn) + '\t' + str(RNA_len) + \
                    '\t' + '\t'.join([str(i) for i in orf[:-1]])

                fickett_score = fickett.fickett_value(orf_seq)
                hexamer_score = FrameKmer.kmer_ratio(
                    orf_seq, 6, 3, coding, noncoding)
                print(orf_id + '\t' + str(fickett_score) + '\t' +
                      str(hexamer_score), file=INFOUT)

                print(">" + orf_id, file=SEQOUT)
                print('\n'.join(wrap(orf_seq, width=options.line_width)),
                      file=SEQOUT)
                orf_sn += 1
            print("%d rows finished\r" % count, end=' ', file=sys.stderr)
        print("\n", file=sys.stderr)

    SEQOUT.close()
    INFOUT.close()

    logging.info("Calculate coding probability ...")
    coding_prediction(
        options.logit_model,
        options.out_file + '.ORF_info.tsv',
        options.out_file)

    if options.mode == 'p':
        logging.info("Select ORF with the highest coding probability ...")
        col_index = 9
    elif options.mode == 'l':
        logging.info("Select the longest ORF ...")
        col_index = 6

    BEST = open((options.out_file + '.ORF_prob.best.tsv'), 'w')
    best_candidates = {}
    for line in open((options.out_file + '.ORF_prob.tsv'), 'r'):
        line = line.strip()
        if line.startswith('ID'):
            print("seq_ID\t" + line, file=BEST)
            continue
        f = line.split('\t')
        seq_id = f[0].split('_ORF_')[0]
        prob = float(f[col_index])
        if seq_id not in best_candidates:
            best_candidates[seq_id] = f
        else:
            if prob > float(best_candidates[seq_id][col_index]):
                best_candidates[seq_id] = f

    for k, v in best_candidates.items():
        print(k + '\t' + '\t'.join(v), file=BEST)

    BEST.close()
    logging.info("Done!")
    finish_up(options.out_file, options.n_top_orf, options.min_orf_len)


def make_hexamer_tab():
    usage = "\n%prog  [options]"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-c", "--cod",
                      action="store",
                      dest="coding_file",
                      help="Coding sequence (must be CDS without UTR, i.e. \
                      from start coden to stop coden) in fasta format")
    parser.add_option("-n",
                      "--noncod",
                      action="store",
                      dest="noncoding_file",
                      help="Noncoding sequences in fasta format")
    (options, args) = parser.parse_args()

    if not options.coding_file and not options.noncoding_file:
        parser.print_help()
        sys.exit(0)
    cod = kmer_freq_file(
        fastafile=options.coding_file,
        word_size=6,
        step_size=3,
        frame=0)
    noncod = kmer_freq_file(
        fastafile=options.noncoding_file,
        word_size=6,
        step_size=1,
        frame=0)

    cod_sum = 0.0
    cod_sum += sum(cod.values())
    noncod_sum = 0.0
    noncod_sum += sum(noncod.values())

    print('hexamer' + '\t' + 'coding' + '\t' + 'noncoding')
    for kmer in cod:
        if 'N' in kmer:
            continue
        print(kmer + '\t' + str(float(cod[kmer]/cod_sum)) +
              '\t' + str(float(noncod[kmer]/noncod_sum)))


def make_logitModel():
    usage = "\n%prog  [options]"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-c",
                      "--cgene",
                      action="store",
                      dest="coding_file",
                      help="Genomic sequnences of protein-coding RNAs in \
                          FASTA or standard 12-column BED format. It is \
                          recommended to use *short* and *unique* sequence \
                          identifiers (such as Ensembl transcript id) in \
                          FASTA and BED file. The input FASTA or BED file \
                          could be a regular text file or compressed file \
                          (*.gz, *.bz2) or accessible URL (http://, https://, \
                          ftp://). When BED file is provided, use the ORF \
                          defined in the BED file (the 7th and 8th columns in \
                          BED file define the positions of 'start codon, \
                          and 'stop codon', respectively). When FASTA file \
                          is provided, searching for the longet ORF. For \
                          well annotated genome, we recommend using BED file \
                          as input because the longest ORF predicted from \
                          RNA sequence might not be the real ORF. If this is \
                          a BED file, reference genome ('-r/--ref') should be \
                          specified.")
    parser.add_option("-n", "--ngene",
                      action="store",
                      dest="noncoding_file",
                      help="Genomic sequences of non-coding RNAs in FASTA or \
                      standard 12-column BED format. It is recommended to \
                      use *short* and *unique* sequence identifiers (such as \
                      Ensembl transcript id) in FASTA and BED file. The \
                      input FASTA or BED file could be a regular text file \
                      or compressed file (*.gz, *.bz2) or accessible URL \
                      (http://, https://, ftp://). If this is a BED file, \
                      reference genome ('-r/--ref') should be specified.")
    parser.add_option("-o", "--outfile",
                      action="store",
                      dest="out_file",
                      help="The prefix of output files.")
    parser.add_option("-x", "--hex",
                      action="store",
                      dest="hexamer_dat",
                      help="Hexamer frequency table. CPAT has prebuilt hexamer\
                      frequency tables for Human, Mouse, Fly, Zebrafish. Run \
                      'make_hexamer_tab.py' to generate this table.")
    parser.add_option("-r", "--ref",
                      action="store",
                      dest="ref_genome",
                      help="Reference genome sequences in FASTA format.\
                      Ignore this option if mRNA sequences file was provided \
                      to '-g'. Reference genome file will be indexed \
                      automatically if the index file  *.fai) does not exist.")
    parser.add_option("-s", "--start",
                      action="store", dest="start_codons",
                      default='ATG',
                      help="Start codon (use 'T' instead of 'U') used to \
                      define the start of open reading frame (ORF). \
                      default=%default")
    parser.add_option("-t", "--stop",
                      action="store",
                      dest="stop_codons",
                      default='TAG, TAA, TGA',
                      help="Stop codon (use 'T' instead of 'U') used to \
                      define the end of open reading frame (ORF). \
                      Multiple stop codons are separated by ','.\
                      default=%default")
    parser.add_option("--min-orf", action="store",
                      type="int",
                      dest="min_orf_len",
                      default=30,
                      help="Minimum ORF length in nucleotides. \
                      default=%default")
    parser.add_option("--log-file",
                      action="store",
                      type="string",
                      dest="log_file",
                      default='make_logitModel_run_info.log',
                      help="Name of log file. default=\"%default\"")
    parser.add_option("--verbose",
                      action="store_true",
                      dest="debug",
                      default=False,
                      help="Logical to determine if detailed running \
                      information is printed to screen.")
    (options, args) = parser.parse_args()

    # check input and output files
    for file in ([options.coding_file, options.noncoding_file,
                 options.out_file, options.hexamer_dat]):
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
    logging.info("Start codons used: \"%s\"", ','.join(start_codons))
    logging.info("Stop codons used: \"%s\"", ','.join(stop_codons))

    # data used to build logit model
    train_data = []
    coding_label = 1
    noncoding_label = 0
    header = ['ID', 'mRNA', 'ORF', 'Fickett', 'Hexamer', 'Label']

    # build hexamer table from hexamer frequency file
    logging.info("Reading hexamer frequency table file: \"%s\"",
                 options.hexamer_dat)
    coding = {}
    noncoding = {}
    for line in open(options.hexamer_dat, 'r'):
        line = line.strip()
        fields = line.split()
        if fields[0] == 'hexamer':
            continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])

    # process protein coding transcripts
    count = 0
    logging.info(
        "Process protein-coding RNA file: \"%s\"", options.coding_file)

    file_format = bed_or_fasta(options.coding_file)
    if file_format == 'UNKNOWN':
        logging.error(
            "Error: unknown file format \"%s\"", options.coding_file)
        parser.print_help()
        sys.exit(0)
    elif file_format == 'BED':
        logging.info(
            "Protein-coding RNA file \"%s\" is in BED format",
            options.coding_file)
        if not options.ref_genome:
            logging.error(
                "Error: Reference genome file must be provided")
            parser.print_help()
            sys.exit(0)
        index_fasta(options.ref_genome)

        for line in ireader.reader(options.coding_file):
            count += 1
            if count % 10 == 0:
                print("%d rows finished\r" % count, end=' ', file=sys.stderr)
            if line.startswith('track'):
                continue
            if line.startswith('#'):
                continue
            if line.startswith('browser'):
                continue

            # option-2: extract CDS directly from BED
            (gene_id, mRNA_size, CDS_size, fickett_score, hexamer) = \
                extract_CDS_from_bed(
                    line,
                    refgenome=options.ref_genome,
                    stt=start_codons,
                    stp=stop_codons,
                    c_tab=coding,
                    g_tab=noncoding,
                    min_orf=options.min_orf_len)

            train_data.append([gene_id,
                               mRNA_size,
                               CDS_size,
                               fickett_score,
                               hexamer,
                               coding_label])
        logging.info("Total %d coding rows finished.", count)
    elif file_format == 'FASTA':
        if options.ref_genome:
            logging.warning(
                "Reference genome sequence [-r] will be ignored when input \
                file is in FASTA format.")
        logging.info("Protein-coding RNA file \"%s\" is in FASTA format"
                     % options.coding_file)
        for sname, seq in FrameKmer.seq_generator(options.coding_file):
            count += 1
            if count % 10 == 0:
                print("%d sequences finished\r" % count, end=' ',
                      file=sys.stderr)
            features = extract_feature_from_seq(
                seq=seq,
                stt=start_codons,
                stp=stop_codons,
                c_tab=coding,
                g_tab=noncoding,
                min_orf=options.min_orf_len)
            if features is None:
                continue
            (mRNA_size, CDS_size, fickett_score, hexamer) = features
            train_data.append([
                sname,
                mRNA_size,
                CDS_size,
                fickett_score,
                hexamer,
                coding_label])
        logging.info("Total %d coding sequences finished.", count)

    # process Non-protein coding transcripts
    count = 0
    logging.info(
        "Process non-coding RNA file: \"%s\"", options.noncoding_file)

    file_format = bed_or_fasta(options.noncoding_file)
    if file_format == 'UNKNOWN':
        logging.error(
            "Error: unknown file format \"%s\"", options.noncoding_file)
        parser.print_help()
        sys.exit(0)
    elif file_format == 'BED':
        logging.info(
            "Non-coding RNA file \"%s\" is in BED format",
            options.noncoding_file)
        if not options.ref_genome:
            logging.error("Error: Reference genome file must be provided")
            parser.print_help()
            sys.exit(0)
        index_fasta(options.ref_genome)

        for line in ireader.reader(options.noncoding_file):
            count += 1
            if count % 10 == 0:
                print("%d genes finished\r" % count, end=' ', file=sys.stderr)
            if line.startswith('track'):
                continue
            if line.startswith('#'):
                continue
            if line.startswith('browser'):
                continue
            fields = line.split()
            if int(fields[1]) != int(fields[6]):
                logging.warning(
                    "This seems to be protein-coding:%s",
                    '\t'.join(fields[0:6]))

            # if not line.strip(): continue
            features = extract_feature_from_bed(
                line,
                refgenome=options.ref_genome,
                stt=start_codons,
                stp=stop_codons,
                c_tab=coding,
                g_tab=noncoding,
                min_orf=options.min_orf_len)
            if features is None:
                logging.warning(
                    "No ORF found for: %s", '\t'.join(line.split()[0:6]))
                continue
            (gene_id, mRNA_size, CDS_size, fickett_score, hexamer) = features
            train_data.append([
                gene_id,
                mRNA_size,
                CDS_size,
                fickett_score,
                hexamer,
                noncoding_label])
        logging.info("Total %d non-coding rows finished.", count)
    elif file_format == 'FASTA':
        if options.ref_genome:
            logging.warning(
                "Reference genome sequence [-r] will be ignored when input \
                file is in FASTA format.")
        logging.info(
            "Non-coding RNA file \"%s\" is in FASTA format",
            options.noncoding_file)
        for sname, seq in FrameKmer.seq_generator(options.noncoding_file):
            count += 1
            if count % 10 == 0:
                print("%d sequences finished\r" % count, end=' ',
                      file=sys.stderr)
            # geneSeq = fa.getSeq(seqID = geneID)
            features = extract_feature_from_seq(
                seq=seq,
                stt=start_codons,
                stp=stop_codons,
                c_tab=coding,
                g_tab=noncoding,
                min_orf=options.min_orf_len)
            if features is None:
                continue
            (mRNA_size, CDS_size, fickett_score, hexamer) = features
            train_data.append([
                sname,
                mRNA_size,
                CDS_size,
                fickett_score,
                hexamer,
                noncoding_label])
        logging.info(f"Total {count} non-coding sequences finished.")
    # writing data
    logging.info("Wrting to \"%s\"", (options.out_file + '.feature.xls'))
    TMP = open(options.out_file + '.feature.xls', 'w')
    print('\t'.join(header), file=TMP)
    for i in train_data:
        print('\t'.join([str(j) for j in i]), file=TMP)
    TMP.close()

    # print("build logi model ...", file=sys.stderr)
    make_logit(options.out_file + '.feature.xls', options.out_file + '.make_logitModel.r', options.out_file + '.logit.RData')
