# add path of this file (e.g., a scripts directory) to sys.path
import copy, io, os, re, subprocess, sys, tempfile, urllib.parse
from pathlib import Path
sys.path.append(str(Path(__file__).resolve(strict=True).parent))

import numpy as np
import pandas as pd

import utils_files, utils

# region ------ General bioinformatics tools

# Source: https://pythonforbiologists.com/dictionaries
# Codons are sorted alphabetically with respect to each amino acid
STANDARD_CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

def rcindex(index, length):
    '''
    Get the position of `index` for a reverse-complemented sequence of length `length`.

    Example: str(Bio.Seq.Seq(seq[index]).complement()) == seq.reverse_complement()[RCIndex(index, len(seq))]
    '''
    return length - index - 1

def reverse_translate(seq, codon_table=None, mode='best_codon', seed=''):
    '''
    Reverse translate amino acid sequence to DNA/RNA using a codon table.

    Args
    - seq: str
        Amino acid sequence to reverse translate
    - codon_table: dict (str -> dict (str -> float)). default=None
        Codon table mapping amino acid to dictionary of codon frequencies.
        Structure: {amino_acid: {codon1: freq1, codon2: freq2, ...}, ...}
        If None, uses STANDARD_CODON_TABLE with evenly distributed frequencies.
    - mode: str. default='best_codon'
        'best_codon': Always choose the most frequent codon. If multiple codons have
          the same frequency, the first codon in the dictionary is chosen.
        'harmonized': Sample codons based on their frequency in the codon table.
    - seed: int or 1-d array_like. default=''
        If of type str, ignored. Otherwise, passed to numpy.random.seed without type checking.
        Seed is set at the start of this function.

    Returns: str
      Reverse-translated sequence.
    '''
    if not isinstance(seed, str):
        np.random.seed(seed)
    assert mode in ('best_codon', 'harmonized')

    if codon_table is None:
        codon_table = {}
        if mode == 'best_codon':
            for codon, aa in STANDARD_CODON_TABLE.items():
                codon_table.setdefault(aa, codon)
        else:
            for codon, aa in STANDARD_CODON_TABLE.items():
                if aa in codon_table:
                    codon_table[aa].update({codon: 1})
                else:
                    codon_table[aa] = {codon: 1}
            for aa, freqs in codon_table.items():
                s = 1/sum(freqs.values())
                codon_table[aa] = {codon: s*freq for codon, freq in freqs.items()}
    assert len(set([type(codon_table[aa]) for aa in codon_table])) == 1

    if mode == 'best_codon':
        if not isinstance(codon_table[list(codon_table.keys())[0]], str):
            codon_table = copy.deepcopy(codon_table)
            for aa, freqs in codon_table.items():
                codon_table[aa] = max(list(codon_table[aa].keys()), key=lambda codon: codon_table[aa][codon])
        return ''.join([codon_table[aa] for aa in seq])

    return ''.join([np.random.choice(list(codon_table[aa].keys()), p=list(codon_table[aa].values())) for aa in seq])

# endregion --- General bioinformatics tools

# region ------ File type columns

BED_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd',
                'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

GFF3_COLNAMES = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

# see https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
DTYPES_GFF3 = {
    'seqid': 'category',
    'source': 'category',
    'type': 'category',
    'start': int,
    'end': int,
    'score': float,
    'strand': pd.CategoricalDtype(categories=['+', '-', '.', '?']),
    'phase': pd.CategoricalDtype(categories=['0', '1', '2', '.']),
    'attributes': str
}

DTYPES_GFF3_NO_SCORE = DTYPES_GFF3 | {'score': pd.CategoricalDtype(categories=['.'])}

# endregion --- File type columns

# region ------ FASTA tools

## FASTA Header Parsers ##

FASTA_HEADER_REGEX_DEFAULT = re.compile(
    r'(?P<name>\S+)\s*' +
    r'(?P<desc>.*)'
)

FASTA_HEADER_REGEX_ENSEMBL_PEPTIDE = re.compile(
    r'(?P<id>\S+)\s+' +
    r'(?P<seqtype>[^:\s]+)' +
    r'(:(?P<status>\S+))?\s+' + 
    r'(?P<coord_system>[^:]+):'
    r'(?P<version>[^:]+):' +
    r'(?P<name>[^:]+):' +
    r'(?P<start>[^:]+):' +
    r'(?P<end>[^:]+):' +
    r'(?P<strand>\S+)\s+' +
    r'gene:(?P<gene>\S+)\s+' +
    r'transcript:(?P<transcript>\S+)\s*' +
    r'(gene_biotype:(?P<gene_biotype>\S+)\s*)*' +
    r'(transcript_biotype:(?P<transcript_biotype>\S+)\s*)*' +
    r'(gene_symbol:(?P<gene_symbol>\S+)\s*)*' +
    r'(description:(?P<description>[^[]+)\s*)*' +
    r'(\[Source:(?P<source>[^;]+);)*' +
    r'(Acc:(?P<accession>[^]]+)\])*'
)

FASTA_HEADER_REGEX_UCSC = re.compile(
    r'(?P<name>\S+)\s+' +
    r'range=(?P<chrom>[^:\s]+):' +
    r'(?P<chromStart>[^-]+)-' + 
    r'(?P<chromEnd>\S+)\s+' + 
    r'5\'pad=(?P<pad5>\S+)\s+' +
    r'3\'pad=(?P<pad3>\S+)\s+' + 
    r'strand=(?P<strand>\S+)\s+' + 
    r'repeatMasking=(?P<repeatMasking>\S+).*'
)

FASTA_HEADER_REGEX_UNIPROT = re.compile(
    r'(?P<db>sp|tr)\|(?P<id>[^|]+)\|(?P<uniprotName>\S+)\s+' +
    r'(?P<proteinName>.*(?=\sOS=)) '
    r'OS=(?P<os>.*(?=\sOX=))\s+' +
    r'OX=(?P<ox>.*(?=\sPE=))\s+' +
    r'PE=(?P<pe>.*(?=\sSV=))\s+' +
    r'SV=(?P<sv>.*)'
)

FASTA_HEADER_REGEX_SGD_PROTEIN = re.compile(
    r'(?P<systematic_name>\S+)\s+' +
    r'(?P<standard_name>\S+)\s+'
    r'SGDID:(?P<SGDID>[^,]+),\s+' +
    r'(Chr )?(?P<chromosome>.+?(?= from))\s+from\s+' +
    r'(?P<start>\d+)-(\d+,\d+[-])*' +
    r'(?P<end>\d+),\s+' +
    r'(Genome Release (?P<release>[^,]+),\s+)?' +
    r'(?P<strand>(reverse complement)?)(,\s+)?' +
    r'(?P<status>[^,]+)' +
    r'(,\s+"(?P<description>[^"]+)")?'
)

def parseDefaultHeader(header, header_prefix='>', pattern=FASTA_HEADER_REGEX_DEFAULT):
    '''
    Parse FASTA header.

    Args
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    - pattern: str or re.Pattern. default=FASTA_HEADER_REGEX_DEFAULT
        Regex pattern with grouped names

    Returns: dict
    - Map of metadata of protein sequence.
    - Keys: depends on `pattern`
    '''

    # process and validate `pattern` argument
    if isinstance(pattern, str):
        pattern = re.compile(pattern)
    assert isinstance(pattern, re.Pattern)

    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]

    # extract key, value pairs from regex match to dict
    try:
        m = pattern.match(header)
        data = m.groupdict()

        # strip whitespace from dict values
        data = {key: value.strip() for key, value in data.items() if value not in ('', None)}
    except AttributeError:
        print(f"Failed to match header: {header}")
        data = {}
    except Exception as err:
        print(f"Failed to parse header: {header}")
        print(f"Error: {err}")
        data = {}

    return data

def parseEnsemblPepHeader(header, **kwargs):
    '''
    Parse Ensembl Peptide FASTA header. Wrapper for parseDefaultHeader().

    Returns: dict
    - Map of metadata of protein sequence.
    - Keys: not all may be present, e.g., the GRCh37.63 release only contains up to 'transcript'
      - id: ENSP ID
      - seqtype: sequence type (pep)
      - status: status of model
          known: can be mapped to species-specific entries UniProt or RefSeq
          novel: cannot be mapped (e.g., predicted based on homology)
      - coord_system: coordinate system (chromosome, supercontig)
      - version: coordinate system version (e.g., GRCh38)
      - name: location name (e.g., 12 - twelfth chromosome)
      - start: start coordinate
      - end: end coordinate
      - strand: strandedness (1, -1)
      - gene: ENSG ID
      - transcript: ENST ID
      - gene_biotype: see https://ensembl.org/info/genome/genebuild/biotypes.html
      - transcript_biotype: see https://ensembl.org/info/genome/genebuild/biotypes.html
      - gene_symbol: HGNC gene symbol
      - description: gene description
      - source: 
      - accession: accession id of sequence from source

    References
    - See the README file in the FTP directory (e.g.,
      http://ftp.ensembl.org/pub/release-63/fasta/homo_sapiens/pep/README) for a
      detailed description of the header line format and the file naming conventions.
    - See https://uswest.ensembl.org/info/data/ftp/index.html for a description of
      the LOCATION (e.g., 'chromosome:NCBI35:1:904515:910768:1') attribute.
    '''

    return parseDefaultHeader(header, pattern=FASTA_HEADER_REGEX_ENSEMBL_PEPTIDE, **kwargs)

def parseSGDProteinHeader(header, **kwargs):
    '''
    Parse SGD Protein FASTA header. Wrapper for parseDefaultHeader().

    Returns: dict
    - Map of metadata of protein sequence.
    - Keys: not all may be present
      - systematic_name: systematic name
      - standard_name: standard name
      - SGDID: SGD unique identifier
      - chromosome: chromsome
      - start: start coordinate
      - end: end coordinate
      - release: genome version
      - strand: strandedness (reverse complement, NA)
      - status: "Verified", "Uncharacterized", or "Dubious"
      - description: description

    References
    - Yeast ORF naming conventions: https://sites.google.com/view/yeastgenome-help/sgd-general-help/glossary
    - Yeast ORF statuses: https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_protein.README
    '''

    return parseDefaultHeader(header, pattern=FASTA_HEADER_REGEX_SGD_PROTEIN, **kwargs)

def parseUniProtHeader(header, header_prefix='>'):
    '''
    Parse UniProt FASTA header. See https://www.uniprot.org/help/fasta-headers.

    Args
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix

    Returns: dict
    - Map of metadata of protein sequence.
    - Keys:
      - db: 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL
      - id: primary accession number of the UniProtKB entry (e.g., P01308)
      - uniprotName: [EntryName] entry name of the UniProtKB entry (e.g., INS_HUMAN)
      - proteinName: [ProteinName] recommended name of the UniProtKB entry as annotated in the RecName field (e.g., Insulin)
      - os: [OrganismName] scientific name of the organism of the UniProtKB entry
      - ox: [OrganismIdentifier] unique identifier of the source organism, assigned by the NCBI
      - gn: [GeneName] first gene name of the UniProtKB entry (may be empty
      - pe: [ProteinExistence] numerical value describing the evidence for the existence of the protein
      - sv: [SequenceVersion] version number of the sequence
    '''

    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]

    # extract gene name if present
    split = re.split(r'GN=(?P<gn>.*)(?=\sPE=)\s+', header)
    m_gn = ''
    if len(split) not in [1,3]:
        raise ValueError
    if len(split) > 1:
        m_gn = split[1]
        header = split[0] + split[2]

    # extract key, value pairs from regex match to dict
    pattern = FASTA_HEADER_REGEX_UNIPROT
    m = pattern.match(header)
    data = m.groupdict()

    # add gene name if present
    data['gn'] = m_gn

    # remove leading/trailing whitespace from each value in dict
    data = {key: value.strip() for key, value in data.items()}
    return data

def parseUCSCHeader(header, header_prefix='>', retainKeys=True, toInt=True):
    '''
    Parse UCSC Table Browser FASTA header.

    Example: hg38_knownGene_ENST00000376838.5_0 range=chr1:11130526-11131568 5'pad=10 3'pad=3 strand=- repeatMasking=none

    Args
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    - retainKeys: bool. default=True
        Retain original FASTA header keys, e.g, 5'pad and 3'pad, instead of the converted 
        valid Python identifiers, e.g., pad5 and pad3.
    - toInt: bool. default=True
        Where possible, convert string values to int. May impact performance.

    Returns: dict
    - Map of metadata of protein sequence.
    - Keys
      - name: sequence name
      - chrom: chromosome (chr#)
      - chromStart: start coordinate (browser format: 1-based start and end)
      - chromEnd: end coordinate (browser format: 1-based start and end)
      - 5'pad: extra bases at the 5' end of the feature
      - 3'pad: extra bases at the 3' end of the feature
      - strand: +/-
      - repeatMasking: mask repeats
        - none: no repeat masking
        - N: repeats are masked to N's
        - lower: repeats are masked to lower case
      - See https://genomebrowser.wustl.edu/goldenPath/help/hgTextHelp.html#FASTA for an older description
        of the FASTA header.
      - In the Table Browser, these options are specified after clicking 'get output'
    '''

    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]

    # extract key, value pairs from regex match to dict
    pattern = FASTA_HEADER_REGEX_UCSC
    m = pattern.match(header)
    data = m.groupdict()

    if retainKeys:
        data['5\'pad'] = data.pop('pad5')
        data['3\'pad'] = data.pop('pad3')

    # remove leading/trailing whitespace from each value in dict
    if toInt:
        # convert str to int if applicable
        for key, value in data.items():
            value = value.strip()
            if utils.isint(value):
                value = int(value)
            data[key] = value
    else:
        data = {key: value.strip() for key, value in data.items()}
    return data

def fastaToDf(file, header_prefix='>', headerParser=parseDefaultHeader, keepRawHeader=None, **kwargs):
    '''
    Parse FASTA file into pandas DataFrame.

    Args
    - file: str or file object
        Path to FASTA file, or file object. Standard compression extensions accepted.
    - header_prefix: str. default='>'
        FASTA header line prefix
    - headerParser: function. default = parseUniProtHeader
        Function to parse header line into dict
    - keepRawHeader: str. default=None
        Column name in which to store original raw header including the header prefix.
    - **kwargs
        Additional keyword arguments to pass to headerParser().

    Returns: pandas.DataFrame
      Rows: protein / DNA entries
      Columns: data about entries. Always includes 'seq' (sequence) column.
    '''

    if isinstance(file, str):
        f = utils_files.createFileObject(file, 'rt')
    else:
        f = file
        if not isinstance(file, (io.IOBase, tempfile._TemporaryFileWrapper)):
            print(f'Could not verify that file {file} is a file object. Continuing anyway...', file=sys.stderr)

    entries = []
    entry = {'seq': ''}
    while True:
        line = f.readline()

        if line == '':
            entries.append(entry)
            break

        if line.startswith(header_prefix):
            # add previous entry to running list of entries
            if entry['seq'] != '':
                # add sequence to entry
                entries.append(entry)

            # parse new entry
            entry = headerParser(line, header_prefix=header_prefix, **kwargs)
            entry['seq'] = ''
            if isinstance(keepRawHeader, str):
                entry[keepRawHeader] = line.strip()
        else:
            entry['seq'] += line.strip()
    f.close()

    # construct pandas DataFrame from list of dicts
    df = pd.DataFrame(entries)
    return df

def dfToFasta(df, file=None, close=True, name_col='name', seq_col='seq', header_prefix='>', spacer='\n'):
    '''
    Write FASTA file from dataframe of sequences.

    Args
    - df: pandas.DataFrame
        Table of names and sequences
    - file: str or file object. default=None
        Path to save FASTA file, or file object. If None, prints to stdout.
    - close: bool. default=True
        Close file object before returning
    - name_col: str. default='name'
        Column to use as FASTA header
    - seq_col: str. default='seq'
        Column with sequences
    - header_prefix: str. default='>'
        FASTA header line prefix
    - spacer: str. default='\n'
        Space between sequences.

    Returns: None
    '''

    try:
        if isinstance(file, str):
            f = utils_files.createFileObject(file, 'w')
        else:
            f = file
            if not isinstance(file, (io.IOBase, tempfile._TemporaryFileWrapper)):
                print(f'Could not verify that file {file} is a file object. Continuing anyway...', file=sys.stderr)

        for _, row in df.iterrows():
            print(header_prefix + row[name_col], row[seq_col], sep='\n', end=spacer, file=f)
    finally:
        if close:
            try:
                f.close()
            except:
                pass

def fastqToDf(file, quality=True, reads_keep=None, reads_ignore=None, close=True, truncate_name=False):
    '''
    Args
    - file: str or file object
        Path to FASTQ file, or file object. Standard compression extensions accepted.
        Assumes file uses exactly 4 lines per read; sequence and quality strings are not split over multiple lines.
    - quality: bool. default=True
        Read quality values
    - reads_keep: iterable of str. default=None
        Only keep these read names
    - reads_ignore: iterable of str. default=None
        Ignore these read names
    - close: bool. default=True
        Only applicable if `file` is a file object. Close file object before returning.
    - truncate_name: bool. default=False
        Truncate name of read at first whitespace.

    Returns: pandas.DataFrame or None
      Rows: FASTQ entries
      Columns: 'name', 'seq'[, 'quality']
      Returns None if file could not be parsed correctly
    '''
    df = None
    if reads_keep is not None:
        reads_keep = set(reads_keep)
    if reads_ignore is not None:
        reads_ignore = set(reads_ignore)
    if isinstance(file, str):
        assert Path(file).exists()
        if file.endswith('.gz'):
            cat = 'zcat'
        elif file.endswith('.bz2'):
            cat = 'bzcat'
        elif file.endswith('.xz'):
            cat = 'xzcat'
        else:
            cat = 'cat'
        if quality:
            awk_prog = r'BEGIN {RS="@"; FS="\n"; OFS="\t"} NR>1 {print $1, $2, $4}'
        else:
            awk_prog = r'BEGIN {RS="@"; FS="\n"; OFS="\t"} NR>1 {print $1, $2}'
        cmd = f'{cat} {file} | awk \'{awk_prog}\''
        if reads_keep is not None:
            cmd += f' | grep -E -e \'{"|".join(map(re.escape, reads_keep))}\''
        if reads_ignore is not None:
            cmd += f' | grep -E -v -e \'{"|".join(map(re.escape, reads_ignore))}\''
        df = pd.read_csv(
            io.StringIO(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout),
            sep='\t', header=None, names=['name', 'seq'] + (['quality'] if quality else []))
        if truncate_name and len(df) > 0:
            df['name'] = df['name'].str.split('\s', expand=True)[0]
    elif isinstance(file, (io.IOBase, tempfile._TemporaryFileWrapper)):
        try:
            entries = []
            entry = {}
            NR = 1
            for line in file:
                if NR == 1:
                    if line[0] != '@':
                        print('FASTQ file not formatted correctly. Entry name must begin with \'@\'.', file=sys.stderr)
                        return None
                    entry['name'] = line[1:].strip()
                    if truncate_name:
                        entry['name'] = entry['name'].split()[0]
                elif NR == 2:
                    entry['seq'] = line.strip()
                elif NR == 4:
                    if quality:
                        entry['quality'] = line.strip()
                    if reads_keep is None or entry['name'] in reads_keep:
                        if reads_ignore is None or not entry['name'] in reads_ignore:
                            entries.append(entry)
                    entry = {}
                    NR = 0
                NR += 1
            if NR != 1:
                print('FASTQ file appears truncated.', file=sys.stderr)
            df = pd.DataFrame(entries, columns=['name', 'seq'] + (['quality'] if quality else []))
        finally:
            if close:
                try:
                    file.close()
                except:
                    pass
    return df


# endregion --- FASTA tools

# region ------ BED tools

def bedSort(bed, tmpColName='chrom_mod', inplace=False):
    '''
    Sort a BED file by chrom, chromStart, chromEnd numerically then lexicographically

    Args
    - bed: pandas.DataFrame
        BED table, column order assumed to be chrom, chromStart, and chromEnd, ...
        Actual names of columns do not matter.
        First (chrom) column assumed to follow UCSC chromosome format: chr1, chr2, ...
    - tmpColName: str. default='chrom_mod'
        Unique temporary column name to use for sorting.
    - inplace: bool. default=False
        True: Perform sorting inplace (edits input DataFrame, which is passed by reference) and return None.
        False: Make a deep copy of the input DataFrame before sorting. Returns sorted DataFrame.

    Returns: pandas.DataFrame or None
      Sorted BED table

    TODO: Handle other chromosome naming formats.
    '''
    assert(tmpColName not in bed.columns)

    def pad_fn(x):
        '''
        Pad 0 to single-digit chromosomes: chr# --> chr0#
        '''
        if x[3:].isdigit() and int(x[3:]) < 10:
            return 'chr0' + x[3:]
        else:
            return x

    if inplace:
        bed[tmpColName] = bed.iloc[:,0].map(pad_fn)
        bed.sort_values(by=[tmpColName, bed.columns[1], bed.columns[2]], inplace=True)
        bed.drop(tmpColName, axis=1, inplace=True)
        return None
    else:
        bed = bed.copy()
        bed[tmpColName] = bed.iloc[:,0].map(pad_fn)
        bed = bed.sort_values(by=[tmpColName, bed.columns[1], bed.columns[2]]).drop(tmpColName, axis=1)
    return bed

def bedOverlapIndices(bed):
    '''
    Find overlapping intervals in a BED file.

    Args
    - bed: pandas.DataFrame
        BED format, must contain first 3 columns

    Returns: list of 2-tuples
      Each tuple represents a pair of overlapping intervals.
      The tuple elements are the DataFrame indices of the intervals.

    Algorithm based on https://www.geeksforgeeks.org/check-if-any-two-intervals-overlap-among-a-given-set-of-intervals/.
    '''

    # make sure all start positions are less than end positions; flip if necessary
    for i in range(bed.shape[0]):
        start = bed.iloc[i,1]
        if start > bed.iloc[i,2]:
            bed.iloc[i,1] = bed.iloc[i,2]
            bed.iloc[i,2] = start

    overlaps = []

    # sort intervals by start
    bed = bed.sort_values(by=[bed.columns[1], bed.columns[2]])

    # if start of an interval < end of previous interval, then there is an overlap
    for i in range(1, bed.shape[0]):
        j = 1
        while bed.iloc[i,1] < bed.iloc[i-j,2] and i >= j:
            overlaps.append(tuple(sorted([bed.index[i-j], bed.index[i]])))
            j += 1

    return overlaps

# endregion --- BED tools

# region ------ SAM tools

# References
# - SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
# - SAM flags tool: https://www.samformat.info/sam-format-flag

# See section 1.4 (The alignment section: mandatory fields) of the SAM specification
SAM_COLUMN_TYPES = {
    'QNAME': str,
    'FLAG': np.uint16,
    'RNAME': str,
    'POS': np.uint32,
    'MAPQ': np.uint8,
    'CIGAR': str,
    'RNEXT': str,
    'PNEXT': np.uint32,
    'TLEN': np.int64,
    'SEQ': str,
    'QUAL': str,
    'OTHER': str
}

def fixSamColtypes(df, copy=True):
    '''
    Change column types of SAM DataFrame to those specified by SAM_COLUMN_TYPES.
    
    Args
    - df: pandas.DataFrame
        DataFrame of SAM data
    - copy: bool. default=True
        Make a deep copy of the Dataframe before modifying column types
    
    Returns: pandas.DataFrame
      DataFrame of SAM data with updated column types.
    '''
    if copy:
        df = df.copy()
    for column, dtype in SAM_COLUMN_TYPES.items():
        if column not in df.columns:
            continue
        try:
            df[column] = df[column].astype(dtype)
        except Exception as e:
            print(e, file=sys.stderr)
            print(f'Failed to convert column {column} to dtype {dtype}.', file=sys.stderr)
    return df

def samToDf(file, exclude=None, use_pandas=True, **kwargs):
    '''
    Parse SAM file into dataframe.
    
    For large SAM files, use the following arguments to return a TextFileReader object for iteration:
    - use_pandas=True
    - iterator=True
    - [optional] chunksize=<int>

    Args
    - file: str or file object
        Path to SAM file, or file object. Standard compression extensions accepted.
    - exclude: dict from int -> set of str. default=None
        Map from column number (0-indexed) to values to exclude.
        Example: To exclude read names (QNAME) of 'abc' or 'def', pass exclude={0: set(['abc', 'def'])}
        Only supported if `use_pandas=False`.
    - use_pandas: bool. default=True
        True: Parse SAM file using pandas.read_csv(). Only the first optional field will be kept.
          Empirically ~3x faster.
        False: All optional fields are kept but concatenated with a single space delimiter.
    - **kwargs
        Keyword arguments to pass to pandas.read_csv().
        If `use_pandas=False`, only supports 'sep', 'comment', and 'nrows' arguments.
        Default values:
        - sep='\t',
        - header=None,
        - index_col=False,
        - names=list(SAM_COLUMN_TYPES.keys()),
        - usecols=list(SAM_COLUMN_TYPES.keys()),
        - dtype=SAM_COLUMN_TYPES,
        - nrows=None,
        - comment='@'

    Returns: pandas.DataFrame
      Columns and types are based on SAM_COLUMN_TYPES. If a column is unable to be
      cast to int type, it remains as a str/object type.
      - 'QNAME': str
      - 'FLAG': int
      - 'RNAME': str
      - 'POS': int
      - 'MAPQ': int
      - 'CIGAR': str
      - 'RNEXT': str
      - 'PNEXT': int
      - 'TLEN': int
      - 'SEQ': str
      - 'QUAL': str
      - 'OTHER': str
    '''
    args = dict(
        sep='\t',
        header=None,
        index_col=False,
        names=list(SAM_COLUMN_TYPES.keys()),
        usecols=list(SAM_COLUMN_TYPES.keys()),
        dtype=SAM_COLUMN_TYPES,
        nrows=None,
        comment='@'
    )
    args.update(kwargs)

    if use_pandas:
        df = pd.read_csv(file, **args)
    else:
        if isinstance(file, str):
            f = utils_files.createFileObject(file, 'r')
        else:
            f = file
            if not isinstance(file, (io.IOBase, tempfile._TemporaryFileWrapper)):
                print(f'Could not verify that file {file} is a file object. Continuing anyway...', file=sys.stderr)
        with f:
            entries = []
            for i, line in enumerate(f):
                if args['nrows'] is not None and i > args['nrows']:
                    break
                if line.startswith(args['comment']):
                    continue
                data = line.strip().split(args['sep'])
                if len(data) > 11:
                    data = data[:11] + [' '.join(data[11:])]
                if len(data) == 11:
                    data.append('')
                if exclude is not None:
                    for idx, values in exclude.items():
                        if data[idx] in values:
                            continue
                entries.append(data)
        df = pd.DataFrame(data=entries, columns=SAM_COLUMN_TYPES.keys())
    fixSamColtypes(df, copy=False)
    
    return df

def bamToDf(file, samtools_path=None, use_tempfile=True, parser=samToDf, verbose=True, **kwargs):
    '''
    Parse BAM file into dataframe. Requires samtools.

    Args
    - file: str
        Path to BAM file.
    - samtools_path: str. default=None
        Path to samtools executable. If None, assumes samtools is in PATH.
    - use_tempfile: bool. default=True
        True: converts BAM file to a temporary SAM file first; uses less memory
        False: converts BAM file to SAM file in memory
    - parser: function. default=samToDf
        Function to parse file once the BAM file has been converted to SAM format.
    - verbose: bool. default=True
        Print path to temporary SAM file if use_tempfile is True.
    - **kwargs
        Additional keyword arguments to pass to parser.

    Returns: depends on `parser`. The default samToDf() parser returns a pandas.DataFrame.
    '''
    if samtools_path is None:
        result = subprocess.run(['where' if os.name == 'nt' else 'which', 'samtools'], capture_output=True)
        if result.returncode != 0:
            raise NotImplementedError('Parsing BAM files relies on samtools, which was not found.')
        samtools_path = result.stdout.splitlines()[0]
    assert Path(file).exists()

    if use_tempfile:
        with tempfile.NamedTemporaryFile(mode='r') as sam_file:
            subprocess.run([samtools_path, 'view', file, '-o', sam_file.name])
            if verbose:
                print(f'Using temporary decompressed SAM file at {sam_file.name}')
            return samToDf(sam_file.name)

    sam_file = io.StringIO(subprocess.run([samtools_path, 'view', file], capture_output=True, text=True).stdout)
    return parser(sam_file)

def dfToSam(df, file, close=True):
    '''
    Write dataframe to SAM file with minimal required header.

    Args
    - df: pandas.DataFrame
        DataFrame of SAM data
    - file: str or file object
        Path to save SAM file, or file object. Standard compression extensions accepted.
    - close: bool. default=True
        Close file object before returning

    Returns: None
    '''
    if isinstance(file, str):
        f = utils_files.createFileObject(file, 'w')
    else:
        f = file
        if not isinstance(file, (io.IOBase, tempfile._TemporaryFileWrapper)):
            print(f'Could not verify that file {file} is a file object. Continuing anyway...', file=sys.stderr)

    try:
        sam_columns_strict = list(SAM_COLUMN_TYPES.keys())
        sam_columns_strict.remove('OTHER')
        other_columns = [col for col in df.columns if col not in sam_columns_strict]
        print('@HD', 'VN:1.6', file=f, sep='\t')
        df_header = df[['RNAME', 'TLEN']].drop_duplicates('RNAME')
        df_header = df_header[df_header['RNAME'] != '*']
        df_header['TLEN'] = abs(df_header['TLEN']).astype(int)
        for _, row in df_header.iterrows():
            print('@SQ', f'SN:{row["RNAME"]}', f'LN:{row["TLEN"]}', file=f, sep='\t')
        df = fixSamColtypes(df, copy=True)
        df[sam_columns_strict + other_columns].to_csv(f, sep='\t', header=False, index=False)
    finally:
        if close:
            try:
                f.close()
            except Exception as e:
                 print(e)

def dfToBam(df, file, samtools_path=None, use_tempfile=True):
    '''
    Write dataframe to BAM file with minimal required header.

    Args
    - df: pandas.DataFrame
        DataFrame of SAM data
    - file: str
        Path to BAM file.
    - samtools_path: str. default=None
        Path to samtools executable. If None, assumes samtools is in PATH.
    - use_tempfile: bool. default=True
        True: writes table to a temporary SAM file first before converting to BAM file; uses less memory
        False: writes table directly to BAM file

    Returns: None
    '''
    if samtools_path is None:
        result = subprocess.run(['where' if os.name == 'nt' else 'which', 'samtools'], capture_output=True)
        if result.returncode != 0:
            raise NotImplementedError('Parsing BAM files relies on samtools, which was not found.')
        samtools_path = result.stdout.splitlines()[0]

    if use_tempfile:
        with tempfile.NamedTemporaryFile(mode='w') as sam_file:
            dfToSam(df, sam_file.name)
            subprocess.run([samtools_path, 'view', '-h', '-b', '-o', file, sam_file.name])
    else:
        with io.StringIO() as buf:
            dfToSam(df, buf, close=False)
            subprocess.run([samtools_path, 'view', '-h', '-b', '-o', file], input=buf.getvalue().encode())

# endregion --- SAM tools

# region ------ Sequence alignment tools

## BLAST
# 
# Hit: entry (sequence) in target database
# High-scoring pair (HSP): regions of significant alignment between query and hit sequences
# 
# BioPython Bio.Blast.NCBIXML parser
# - Each query corresponds to a record
# - Each record contains zero or more alignments (hits)
# - Each alignment contains one or more HSPs

def blastToDf(records):
    '''
    Convert BLAST records to DataFrame.

    Args
    - records: iterator of Bio.Blast.Record.Blast
        BLAST records (each record corresponds to 1 query)

    Returns: pandas.DataFrame or None (no HSPs in records)
    - Rows: HSPs
    - Columns
      - align_length: int
          alignment length
      - bits: float
          score (in bits)
      - expect: float
          E-value
      - frame: 2-tuple of int
          (translation frame of query, translation frame of subject)
      - gaps: int
          number of gaps
      - hit_id: str
          SeqId of subject
          Ex: gi|930125|emb|CAA35112.1|
      - identities: int
          number of identities
      - match: str
          formating middle line
          Ex: +ACD CRKKK+KC
      - num_alignments: ??
      - positives: int
          number of positives
      - query: str
          alignment string for the query (with gaps)
          Ex: RACDQCRKKKIKC
      - query_def: str.
          Definition line of query.
          Ex: Rgt1p [Saccharomyces cerevisiae S288C]
      - query_end: int.
      - query_id: str
          SeqId of query
          Ex: NP_015333.1
      - query_ start: int
      - sbjct: str
          alignment string for subject (with gaps)
          Ex: QACDACRKKKLKC
      - sbjct_end: int
      - sbjct_start: int
      - score: float
      - strand: 2-tuple of ??

    Reference
    - See https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd for XML DTD schema

    TODO: add support for iterator of Bio.SearchIO._model.query.QueryResult
    '''
    hsps = [] # high-scoring pairs
    for record in records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                d = vars(hsp)
                d.update({'hit_id': alignment.hit_id, 'query_def': record.query, 'query_id': record.query_id})
                hsps.append(d)
    if len(hsps) == 0:
        return None
    return pd.DataFrame(hsps)

# endregion --- Sequence alignment tools

# region --- GTF/GFF tools

# see https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
DTYPES_GFF3 = {
    'seqid': 'category',
    'source': 'category',
    'type': 'category',
    'start': int,
    'end': int,
    'score': float,
    'strand': pd.CategoricalDtype(categories=['+', '-', '.', '?']),
    'phase': pd.CategoricalDtype(categories=['0', '1', '2', '.']),
    'attributes': str
}

DTYPES_GFF3_NO_SCORE = DTYPES_GFF3 | {'score': pd.CategoricalDtype(categories=['.'])}

def extract_gtf_attribute(
    s_attributes: pd.Series,
    attribute: str,
    table_type: str = 'GTF'
) -> pd.Series:
    '''
    Extract attribute from GTF/GFF attributes column.

    Args
    - s_attributes: GTF/GFF attributes column
    - attribute: name of attribute to extract
    - table_type: 'GTF' or 'GFF'

    Returns: extracted attribute values
    '''
    if table_type == 'GTF':
        regex = r'(?:^| )' + attribute + r' "([^"]+)"'
        return s_attributes.str.extract(regex, expand=False)
    elif table_type == 'GFF':
        regex = r'(?:^|;)' + attribute + r'=([^,=;]+)'
        # In GFF attributes, tags or values containing ",=;" characters are encoded using URL escape rules
        # see https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
        extracted = s_attributes.str.extract(regex, expand=False)
        extracted.loc[~extracted.isna()] = extracted.loc[~extracted.isna()].map(urllib.parse.unquote)
        return extracted
    else:
        raise ValueError(f"table_type {table_type} not understood; must be either 'GTF' or 'GFF'")

def add_columns_from_attributes(
    df: pd.DataFrame,
    attributes: None | list[str] = None,
    rename_attributes: None | dict[str, str] = None,
    col_attributes: str = 'attributes',
    table_type: str = 'GTF'
) -> pd.DataFrame:
    '''
    Add attributes from the attributes column of a GTF table to its own column.

    Args
    - df: GTF/GFF table
    - attributes: attributes to extract from attributes column
    - rename_attributes: mapping of attributes to new column names
    - col_attributes: name of attributes column
    - table_type: whether attributes are given in GFF format (tag1=value1;tag2=value2)
        or GTF format (tag1 "value1"; tag2 "value2")

    Returns: table with the extracted attributes as new columns
    '''
    if attributes is None:
        return df
    for attribute in attributes:
        if rename_attributes:
            new_col_name = rename_attributes[attribute]
        else:
            new_col_name = attribute
        df[new_col_name] = extract_gtf_attribute(df[col_attributes], attribute, table_type)
    return df

# endregion --- GTF/GFF tools