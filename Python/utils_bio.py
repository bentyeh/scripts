# add path of this file (e.g., a scripts directory) to sys.path
import sys, os
sys.path.append(os.path.dirname(__file__))

import io, re, warnings
from collections.abc import Iterable
from multiprocessing import Pool

import pandas as pd
import Bio.Entrez

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
    r'(?P<name>[^\s]+)\s*' +
    r'(?P<desc>.*)'
)

FASTA_HEADER_REGEX_ENSEMBL_PEPTIDE = re.compile(
    r'(?P<id>[^\s]+)\s+' +
    r'(?P<seqtype>[^:\s]+)' +
    r'(?P<status>[^\s]*)\s+' + 
    r'(?P<coord_system>[^:]+):'
    r'(?P<version>[^:]+):' +
    r'(?P<name>[^:]+):' +
    r'(?P<start>[^:]+):' +
    r'(?P<end>[^:]+):' +
    r'(?P<strand>[^\s]+)\s+' +
    r'gene:(?P<gene>[^\s]+)\s+' +
    r'transcript:(?P<transcript>[^\s]+)\s*' +
    r'(gene_biotype:(?P<gene_biotype>[^\s]+)\s*)*' +
    r'(transcript_biotype:(?P<transcript_biotype>[^\s]+)\s*)*' +
    r'(gene_symbol:(?P<gene_symbol>[^\s]+)\s*)*' +
    r'(description:(?P<description>[^[]+)\s*)*' +
    r'(\[Source:(?P<source>[^;]+);)*' +
    r'(Acc:(?P<accession>[^]]+)\])*'
)

FASTA_HEADER_REGEX_UCSC = re.compile(
    r'(?P<name>[^\s]+)\s+' +
    r'range=(?P<chrom>[^:\s]+):' +
    r'(?P<chromStart>[^-]+)-' + 
    r'(?P<chromEnd>[^\s]+)\s+' + 
    r'5\'pad=(?P<pad5>[^\s]+)\s+' +
    r'3\'pad=(?P<pad3>[^\s]+)\s+' + 
    r'strand=(?P<strand>[^\s]+)\s+' + 
    r'repeatMasking=(?P<repeatMasking>[^\s]+).*'
)

FASTA_HEADER_REGEX_UNIPROT = re.compile(
    r'(?P<db>sp|tr)\|(?P<id>[^|]+)\|(?P<uniprotName>\S+)\s+' +
    r'(?P<proteinName>.*(?=\sOS=)) '
    r'OS=(?P<os>.*(?=\sOX=))\s+' +
    r'OX=(?P<ox>.*(?=\sPE=))\s+' +
    r'PE=(?P<pe>.*(?=\sSV=))\s+' +
    r'SV=(?P<sv>.*)'
)

def parseDefaultHeader(header, header_prefix='>'):
    '''
    Parse FASTA header.
    
    Args
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    
    Returns: dict
    - Map of metadata of protein sequence.
    - Keys: name, desc
    '''

    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]
    
    # extract key, value pairs from regex match to dict
    p = FASTA_HEADER_REGEX_DEFAULT
    m = p.match(header)
    data = m.groupdict()
    
    # strip whitespace from dict values
    data = {key: value.strip() for key, value in data.items() if value is not ''}
    return data

def parseEnsemblPepHeader(header, header_prefix='>'):
    '''
    Parse Ensembl Peptide FASTA header.
    
    Args
    - header: str
        FASTA header line
    - header_prefix: str. default='>'
        FASTA header line prefix
    
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
    
    # strip whitespace and prefix
    header = header.strip()
    if header.startswith(header_prefix):
        header = header[len(header_prefix):]
    
    # extract key, value pairs from regex match to dict
    p = FASTA_HEADER_REGEX_ENSEMBL_PEPTIDE
    m = p.match(header)
    data = m.groupdict()
    
    # strip whitespace from dict values
    data = {key: value.strip() for key, value in data.items() if value is not None}
    return data

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
      Keys: db, id, uniprotName, proteinName, os, ox, gn (may be empty), pe, sv
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
    p = FASTA_HEADER_REGEX_UNIPROT
    m = p.match(header)
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
    
    Args:
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
    p = FASTA_HEADER_REGEX_UCSC
    m = p.match(header)
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

def fastaToDF(file, save='', header_prefix='>', headerParser=parseUniProtHeader, **kwargs):
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
        else:
            entry['seq'] += line.strip()
    f.close()
    
    # construct pandas DataFrame from list of dicts
    df = pd.DataFrame(entries)
    return df

# endregion --- FASTA tools

# region --- BED tools

def bedSort(bed, tmpColName='chrom_mod'):
    '''
    Sort a BED file by chrom, chromStart, chromEnd numerically then lexicographically

    Args
    - bed: pandas.DataFrame
        BED table, must contain columns chrom, chromStart, and chromEnd.
        Assumes UCSC chromosome format: chr1, chr2, ...
    - tmpColName: str. default='chrom_mod'
        Temporary column name to use for sorting.

    Returns: pandas.DataFrame
      Sorted BED table
    
    TODO: Handle other chromosome naming formats.
    '''

    def pad_fn(x):
        '''
        Pad 0 to single-digit chromosomes: chr# --> chr0#
        '''
        if x[3:].isdigit() and int(x[3:]) < 10:
            return 'chr0' + x[3:]
        else:
            return x

    bed['chrom_mod'] = bed['chrom'].map(pad_fn)
    bed = bed.sort_values(by=['chrom_mod', 'chromStart', 'chromEnd']).drop('chrom_mod', axis=1)
    return bed

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

# endregion --- Database identifier conversion tools