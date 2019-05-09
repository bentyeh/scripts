# add path of this file (e.g., a scripts directory) to sys.path
import sys, os
sys.path.append(os.path.dirname(__file__))

import io, re, warnings
from collections.abc import Iterable
from multiprocessing import Pool

import pandas as pd
import Bio.Entrez
import requests

import utils_files, utils

# region --- General bioinformatics tools

def rcindex(index, length):
    '''
    Get the position of `index` for a reverse-complemented sequence of length `length`.

    Example: str(Bio.Seq.Seq(seq[index]).complement()) == seq.reverse_complement()[RCIndex(index, len(seq))]
    '''
    return length - index - 1

# endregion --- General bioinformatics tools

# region --- File type columns

BED_COLNAMES = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd',
                'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

GFF3_COLNAMES = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

# endregion --- File type columns

# region --- FASTA tools

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

def fastaToDF(file, save='', header_prefix='>', headerParser=parseDefaultHeader, keepRawHeader=None, **kwargs):
    '''
    Parse FASTA file into pandas DataFrame.

    Args
    - file: str or io.IOBase
        Path to FASTA file, or file object. Gzip-compressed files with extension '.gz'
        are accepted.
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
        f = utils_files.createFileObject(file)
    elif isinstance(file, io.IOBase):
        f = file
    else:
        raise ValueError('`file` must be a string or file object')

    entries = []
    entry = {'seq': ''}
    while True:
        line = f.readline()

        if line == '':
            entries.append(entry)
            break

        if line.startswith(header_prefix):
            # add previous entry to running list of entries
            if entry['seq'] is not '':
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

def dfToFasta(df, file=None, name_col='name', seq_col='seq', header_prefix='>'):
    '''
    Args
    - df: pandas.DataFrame
        Table of names and sequences
    - file: str. default=None
        Path to save FASTA file
    - name_col: str. default='name'
        Column to use as FASTA header
    - seq_col: str. default='seq'
        Column with sequences
    - header_prefix: str. default='>'
        FASTA header line prefix
    '''

    if isinstance(file, str):
        f = utils_files.createFileObject(file, 'wt')
    elif isinstance(file, io.IOBase):
        f = file
    else:
        raise ValueError('`file` must be a string or file object')

    for _, row in df.iterrows():
        print(header_prefix + row[name_col], row[seq_col], sep='\n', end='\n\n', file=f)
    f.close()

# endregion --- FASTA tools

# region --- BED tools

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
        Perform sorting inplace (edits input DataFrame). Otherwise, makes a copy of the input DataFrame.

    Returns: pandas.DataFrame
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

    if not inplace:
        bed = bed.copy()
    bed[tmpColName] = bed.iloc[:,0].map(pad_fn)
    bed = bed.sort_values(by=[tmpColName, bed.columns[1], bed.columns[2]], inplace=True)
    bed = bed.drop(tmpColName, axis=1, inplace=True)
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

# region --- Database identifier conversion tools

## NCBI

def search_NCBIGene(term, email=None, useDefaults=True, useSingleIndirectMatch=True, verbose=True, **kwargs):
    '''
    Search official human gene names and aliases in NCBI Gene database for a match to term, returning official IDs,
    names, and aliases.

    Dependencies: Biopython

    Args
    - term: iterable of str
        Gene name / alias
    - email: str. default=None
        Email registered with NCBI. If not None, sets Bio.Entrez.email. If Bio.Entrez.email is None, defaults to
        'A.N.Other@example.com'.
    - useDefaults: bool. default=True
        Apply the following fields to the query:
        - orgn: 'Homo sapiens'
        - prop: 'alive'
    - useSingleIndirectMatch: bool. default=True
        If the Entrez Gene Database query returns only 1 NCBI Gene ID, use the match even if the
        term does not exactly match the gene symbol or an alias.
    - verbose: bool. default=True
        Print submitted query.
    - kwargs
        Additional fields to add to the query. If multiple values are given, an OR relationship is assumed.
        Use keyword 'verbatim' to directly add text to append to the query.
        Common fields: chr (Chromosome), go (Gene Ontology), pmid (PubMed ID), prop (Properties), taxid (Taxonomy ID)
        References:
        - https://www.ncbi.nlm.nih.gov/books/NBK3841/#_EntrezGene_Query_Tips_How_to_submit_deta_
            Official NCBI Help page.
        - https://www.ncbi.nlm.nih.gov/books/NBK49540/
            Fields (and abbreviated fields) available for sequence databases (Nucleotide, Protein, EST, GSS), but most
            of them apply to the NCBI Gene database as well.

    Returns: dict: int -> dict: str -> str or list of str
      {<NCBI Gene ID>: {'symbol': <gene symbol>, 'name': <gene name>, 'aliases': [aliases]}}
    '''

    # Check email registered in Biopython
    if email is not None:
        Bio.Entrez.email = email
    if Bio.Entrez.email is None:
        Bio.Entrez.email = 'A.N.Other@example.com'

    # construct query parameters
    params = {}
    defaults = {
        'orgn': 'Homo sapiens',
        'prop': 'alive'
    }
    if useDefaults:
        params.update(defaults)
    params.update(kwargs)

    # build query
    query = '{}[gene]'.format(term)
    verbatim = ''
    if 'verbatim' in params:
        verbatim = params.pop('verbatim')
    for key, values in params.items():
        if isinstance(values, str) or not isinstance(values, Iterable):
            values = set([values])
        query += ' AND ({})'.format(' OR '.join(['{}[{}]'.format(value, key) for value in values]))
    query += verbatim
    if verbose:
        print(query)

    # initialize return variable
    matches = {}

    # submit query
    handle = Bio.Entrez.esearch(db="gene", term=query)
    idList = Bio.Entrez.read(handle)['IdList']
    for id in idList:
        handle = Bio.Entrez.esummary(db='gene', id=id)
        record = Bio.Entrez.read(handle)
        symbol = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        name = record['DocumentSummarySet']['DocumentSummary'][0]['Description']
        aliases = record['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases'].split(', ')
        terms = [match.lower() for match in [name] + aliases]
        if (term.lower() in terms):
            matches[id] = int(id)
            matches[id] = {'symbol': symbol, 'name': name, 'aliases': aliases}
    if useSingleIndirectMatch and len(matches) == 0 and len(idList) == 1:
        handle = Bio.Entrez.esummary(db='gene', id=id)
        record = Bio.Entrez.read(handle)
        symbol = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        name = record['DocumentSummarySet']['DocumentSummary'][0]['Description']
        aliases = record['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases'].split(', ')
        matches[id] = {'symbol': symbol, 'name': name, 'aliases': aliases}
        if verbose:
            print(f'{term} failed to match any name or alias, returning the single unique result provided by Entrez.')
    return matches

def multisearch_NCBIGene(terms, nProc=None, **kwargs):
    '''
    Search official human gene names and aliases in NCBI Gene database for terms.

    Args
    - terms: iterable of str
        Terms to lookup
    - email: str
        email registered with NCBI
    - nProc: int. default=None
        Number of processes to use. If None, uses as many processes as available CPUs.

    Returns: dict: str -> dict: int -> dict: str -> str or list of str
      {<term>: <NCBI Gene ID>: {'symbol': <gene symbol>, 'name': <gene name>, 'aliases': [aliases]}}
    '''

    if nProc is None:
        try:
            nProc = len(os.sched_getaffinity(0))
        except AttributeError:
            warnings.warn('Unable to get accurate number of available CPUs. Assuming all CPUs are available.')
            nProc = os.cpu_count()
    assert isinstance(nProc, int) and nProc > 0

    results = []
    with Pool(processes=nProc) as pool:
        for i in range(len(terms)):
            results.append(pool.apply_async(search_NCBIGene, (terms[i],), kwargs))
        return {terms[i]: results[i].get() for i in range(len(terms))}

## UniProt

# See https://www.uniprot.org/help/accession_numbers
UNIPROT_ACCESSION_REGEX = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

def convertUniProt(source, to, query, toDf=True):
    '''
    Convert database identifiers using UniProt's Retrieve/ID mapping service.

    Args
    - source: str
        Query identifier type.
    - to: str
        Target identifier type.
    - query: str
        Identifiers separated by spaces or new lines.
    - toDf: bool. default=True
        Return results as pandas DataFrame.

    Returns: str or pandas.DataFrame
      UniProt's mapping service returns a tab-delimited table with column headers "From" and "To".
      If `toDf` is True, the raw text is returned. Otherwise, it is parsed into a DataFrame.

    Reference: https://www.uniprot.org/help/api_idmapping
    '''

    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': source,
        'to': to,
        'format': 'tab',
        'query': query
    }

    r = requests.post(url=url, data=params)
    if not r.ok:
        r.raise_for_status()

    if toDf:
        return pd.read_csv(io.StringIO(r.text), sep='\t', header=0, index_col=False)
    return r.text

# endregion --- Database identifier conversion tools

# region --- Sequence alignment tools

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
    - records: list of Bio.Blast.Record
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