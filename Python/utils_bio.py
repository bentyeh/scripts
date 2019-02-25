import io, re
import pandas as pd
import utils_files

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
    data = {key: value.strip() for key, value in data.items() if value is not None}
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
        raise
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
            if isint(value):
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

def searchNCBIGeneNames(term, email, useSingleIndirectMatch = True):
    '''
    Search official human gene names and aliases in NCBI Gene database for a match to term, returning offical names and IDs.
    Dependencies: Biopython
    
    Args:
    - term: str
        gene name / alias
    - email: str
        email registered with NCBI
    - useSingleIndirectMatch
        If the Entrez Gene Database query returns only 1 NCBI Gene ID, even if the term does not exactly match the gene symbol
        or an alias, use the match.
    
    Returns: dict: str -> list
        "names": list of matched official gene name(s)
        "ids": list of NCBI Gene IDs corresponding to matched official gene name(s)
    '''
    from Bio import Entrez

    names, ids = [], []
    Entrez.email = email
    handle = Entrez.esearch(db="gene", term='(' + term + '[gene]) AND (Homo sapiens[orgn]) AND alive[prop] NOT newentry[gene]')
    idList = Entrez.read(handle)['IdList']
    for id in idList:
        handle = Entrez.esummary(db='gene', id=id)
        record = Entrez.read(handle)
        name = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        aliases = record['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases'].split(', ')
        if (term in [name] + aliases):
            names.append(name)
            ids.append(id)
    if useSingleIndirectMatch:
        if len(names) == 0 and len(idList) == 1:
            ids.append(idList[0])
            handle = Entrez.esummary(db='gene', id=id)
            record = Entrez.read(handle)
            names.append(record['DocumentSummarySet']['DocumentSummary'][0]['Name'])
    return({"names": names, "ids": ids})

def mtSearchNCBIGeneNames(terms, email, nThreads = None):
    '''
    Search official human gene names and aliases in NCBI Gene database for terms.
    
    Args
    - terms: list of str
        terms to lookup
    - email: str
        email registered with NCBI
    - nThreads: int
        None: Uses a ThreadPool of a default number of threads as returned by os.cpu_count()
        1+: Uses a ThreadPool of nThreads
    
    Returns: list of str
      Where no matches found, returns the empty string. Otherwise, returns the first matched official gene symbol.
    '''
    from multiprocessing.pool import ThreadPool
    
    dict_results = []
    gene_symbols = []
    pool = ThreadPool(nThreads)
    print("Using {:d} threads...".format(pool._processes))
    for i in range(len(terms)):
        dict_results.append(pool.apply_async(searchNCBIGeneNames, (terms[i], email)))
    pool.close()
    pool.join()
    for i in range(len(terms)):
        match = dict_results[i].get()
        if len(match["names"]) == 0:
            gene_symbols.append("")
        elif terms[i] in match["names"]:
            gene_symbols.append(terms[i])
        else:
            gene_symbols.append(match["names"][0])
    return(gene_symbols)
