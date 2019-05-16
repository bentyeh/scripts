# add path of this file (e.g., a scripts directory) to sys.path
import sys, os
sys.path.append(os.path.dirname(__file__))

import io, multiprocessing, re, time
from collections.abc import Iterable

import pandas as pd
import requests
import Bio.Blast.NCBIWWW
import Bio.SearchIO

import utils_files
import utils_bio

# region ------ UniProt

PROTEINS_HEADERS = {
    'xml':   'application/xml',
    'json':  'application/json',
    'fasta': 'text/x-fasta',
    'flat':  'text/x-flatfile',
    'gff':   'text/x-gff',
    'peff':  'text/x-peff'
}

PROTEINS_BASEURL = 'https://www.ebi.ac.uk/proteins/api'

def get_UniProt(uniprot_id, export, toDf=False):
    '''
    Get UniProt entry.

    Args
    - uniprot_id: str
        UniProt Accession ID (e.g., P12345).
    - export: str.
        Format: 'txt', 'fasta', 'xml', 'rdf', 'gff'
    - toDf: bool. default=False
        Parse data into a pandas.DataFrame. Currently only 'fasta' and 'gff' export types are supported.

    Returns: str or pandas.DataFrame
    '''

    assert export in ('txt', 'fasta', 'xml', 'rdf', 'gff')
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.{export}'
    headers = {'accept': 'text/html'}
    r = requests.get(url=url, headers=headers)
    if not r.ok:
        r.raise_for_status()

    if toDf:
        if export == 'fasta':
            return utils_bio.fastaToDF(io.StringIO(r.text), headerParser=utils_bio.parseUniProtHeader)
        if export == 'gff':
            return pd.read_csv(io.StringIO(r.text), sep='\t', header=None, comment='#',
                               names=utils_bio.GFF3_COLNAMES, index_col=False)
    else:
        return r.text

def get_Proteins(service='/proteins', entry=None, export='json', verbose=False, **kwargs):
    '''
    Retrieve UniProt protein entries using the Proteins API.
    See https://www.ebi.ac.uk/proteins/api/doc/.

    Args
    - service: str. default='/proteins'
        Service requested through REST interface. See examples below.
    - entry: dict. default=None
        Desired entry for services that do not take key/value pairs. See examples below.
    - export: str. default='json'
        Response content type: 'fasta', 'json', 'xml', 'flat' (UniProt text format), 'gff', or 'peff'.
        All services can return XML or JSON formatted results. Other return types are service-specific.
    - verbose: bool. default=False
        Print request URL.
    - **kwargs:
        GET parameters (given as key/value pairs in the URL after a question mark).
        Common parameters listed below:
        - accession: comma-separated UniProt accessions, up to 100
        - offset: page starting point, default=0
        - size: page size, default=100. If -1, returns all records and offset will be ignored

    Returns: str

    Examples
    - Request single protein entry, FASTA format
        Desired GET request: https://www.ebi.ac.uk/proteins/api/proteins/P24928
        --> get_Proteins(service='/proteins/{accession}', entry={'accession': 'P24928'}, export='fasta')
    - Request multiple protein entries, JSON format
        Desired GET request: https://www.ebi.ac.uk/proteins/api/proteins?accession=P24928,P30876
        --> get_Proteins(service='/proteins', entry=None, export='json', accession='P24928,P30876')
    '''

    if re.search(r'\{.*\}', service) and entry is not None:
        service = service.format(**entry)
    url = PROTEINS_BASEURL + service
    headers = {'Accept': PROTEINS_HEADERS[export]}
    r = requests.get(url, params=kwargs, headers=headers)

    if verbose:
        print(r.url)

    if not r.ok:
        r.raise_for_status()
    return r.text

# endregion --- UniProt

# region ------ InterPro

URL_InterPro_protein = 'https://www.ebi.ac.uk/interpro/protein/'

regex_InterPro_GOmap = re.compile(
    r'InterPro:(?P<InterPro_ID>[^\s]+)\s+' +
    r'(?P<InterPro_desc>[^>]+)\s+[>]\s+' + 
    r'GO:(?P<GO_desc>[^;]+)\s+[;]\s+' +
    r'(?P<GO_ID>GO[:]\d+)'
)

def get_InterPro_protein(uniprot_id, export='tsv', toDf=True):
    '''
    Get InterPro protein sequence or entry annotation table.

    Args
    - uniprot_id: str
    - export: str. default='tsv'
        'tsv' or 'fasta'
    - toDF: bool. default=True
        Parse data into a pandas.DataFrame. Otherwise returns the raw text.
        - export=='tsv': parse TSV format file directly into data frame
        - export=='fasta': parse FASTA to pandas.DataFrame using utils_bio.fastaToDF()

    Returns: str or pandas.DataFrame
    '''

    assert export in ('tsv', 'fasta')

    url = URL_InterPro_protein + uniprot_id
    headers = {'accept': 'text/html'}
    params = {'export': export}
    r = requests.get(url=url, params=params, headers=headers)
    if not r.ok:
        r.raise_for_status()

    # strip whitespace from each line
    # - lines without data in every column often have too many tab ('\t') separators
    #   such that parsing into a pandas.DataFrame will throw an error
    text = '\n'.join([line.strip() for line in r.text.strip().splitlines()])

    if toDf:
        if export == 'tsv':
            return pd.read_csv(io.StringIO(text), sep='\t', index_col=False)
        else:
            return utils_bio.fastaToDF(io.StringIO(r.text), headerParser=utils_bio.parseDefaultHeader)
    return text

def process_InterPro_GOmap(path):
    '''
    Process InterPro entry to GO term mapping to a data frame.

    Args
    - path: str
        Path to interpro2go file

    Returns: pandas.DataFrame
      Columns: InterPro_ID, InterPro_desc, GO_desc, GO_ID
    '''

    lines = utils_files.readFile(path)
    result = []
    for line in lines:
        m = regex_InterPro_GOmap.match(line)
        if m is not None:
            data = m.groupdict()
            data = {key: value.strip() for key, value in data.items()}
            result.append(data)
    return pd.DataFrame(result)

# endregion --- InterPro

# region ------ QuickGO

GO_RELATIONS = ['is_a', 'part_of', 'occurs_in', 'regulates']

def parse_QuickGO_JSON(results):
    '''
    Parse JSON body returned by QuickGO 'search' API into pandas.DataFrame with same columns
    as returned by the 'download' API

    Arg: list of dict

    Returns: pandas.DataFrame
      Columns will have the same names and order as returned by the QuickGO 'download' API
    '''

    goAspect_map = {'cellular_component': 'C', 'biological_process': 'P', 'molecular_function': 'F'}

    colNames_map = OrderedDict([
        # entries in order by 'search' API method (see getQuickGO())
        ('geneProductDb', 'GENE PRODUCT DB'),
        ('geneProductId', 'GENE PRODUCT ID'),
        ('symbol', 'SYMBOL'),
        ('qualifier', 'QUALIFIER'),
        ('goId', 'GO TERM'),
        ('goAspect', 'GO ASPECT'),
        ('evidenceCode', 'ECO ID'),
        ('goEvidence', 'GO EVIDENCE CODE'),
        ('reference', 'REFERENCE'),
        ('withFrom', 'WITH/FROM'),
        ('taxonId', 'TAXON ID'),
        ('assignedBy', 'ASSIGNED BY'),
        ('extensions', 'ANNOTATION EXTENSION'),
        ('date', 'DATE'),

        # additional entries retrieveable only via 'download' method
        ('goName', 'GO NAME'),
        ('synonyms', 'SYNONYMS')
    ])

    for i in range(len(results)):
        # extract 'withFrom' dictionaries into strings
        withFrom_list = []
        if results[i]['withFrom'] is not None:
            for j in range(len(results[i]['withFrom'])):
                for k in range(len(results[i]['withFrom'][j]['connectedXrefs'])):
                    withFrom_list.append(':'.join([results[i]['withFrom'][j]['connectedXrefs'][0]['db'],
                                                   results[i]['withFrom'][j]['connectedXrefs'][0]['id']]))
        results[i]['withFrom'] = '|'.join(withFrom_list)

        # extract 'geneProductId' into 'GENE PRODUCT DB' and 'GENE PRODUCT ID'
        results[i]['geneProductDb'], results[i]['geneProductId'] = results[i]['geneProductId'].split(':')

        # convert 'goAspect' to single-character symbol
        results[i]['goAspect'] = goAspect_map[results[i]['goAspect']]

    # rename and reorder columns
    df = pd.DataFrame(results).rename(mapper=colNames_map, axis='columns')
    df = df[list(colNames_map.values())]

    # convert 'DATE' to date format; convert None values to np.nan
    df['DATE'] = pd.to_datetime(df['DATE'], yearfirst=True, format='%Y%m%d')
    df.loc[df['ANNOTATION EXTENSION'].isnull(), 'ANNOTATION EXTENSION'] = np.nan
    return df

def get_QuickGO_annotations(goIds, useDefaults=True, method='opt', parseResponse=True,
                            attempts=5, sleep=0.5, nProc=1, **kwargs):
    '''
    Download annotations from QuickGO. See https://www.ebi.ac.uk/QuickGO/api/index.html.

    Args
    - goIds: str
        Comma-separated string of GO IDs (e.g., 'GO:0016592' for mediator complex)
    - useDefaults: bool. default=True
        Apply the following 'default' filters:
        - taxonId: 9606 (Homo sapiens)
        - geneProductType: 'protein'
        - geneProductSubset: 'Swiss-Prot' (only reviewed UniProtKB entries)
        - proteome: 'gcrpCan' (Gene Centric Reference Proteome Canonical)
        - goUsage: 'descendants'
        - goUsageRelationships: 'is_a,part_of,occurs_in' (excludes 'regulates')
        - limit: 100 (maximum limit per page for 'search' API)
        - downloadLimit: 50000 (maximum limit for 'download' API)
    - method: str. default='opt'
        'search': Use QuickGO 'search' API, which returns a JSON response body.
          Required for large (>50000 returned entries) query
        'download': Use QuickGO 'download' API, which returns a text (TSV) response body.
          Maximum 50000 returned entries.
        'opt': Try 'search'. If the number of hits is < 50000, switch to 'download'.
    - parseResponse: bool. default=True
        True: parse response into pandas.DataFrame
        False: return respose bodies
    - attempts: int. default=5
        Number of attempts to retry a request upon error
    - sleep: float. default=0.5
        Seconds to sleep for in between each attempt
    - nProc: int. default=1
        Number of processes to use. Error handling is not implemented if nProc > 1.
    - **kwargs
        Additional parameters to pass to send in the body of the HTTP request

    Returns: pandas.DataFrame, list of requests.Response, or None
      `parseResponse` == True --> pandas.DataFrame
      `parseResponse` == False --> list of requests.Response
      Error in request after `attempt` attempts --> None

    Notes
    - The following parameters do not seem to have any effect using the download API:
      includeFields, limit
    '''

    # validate arguments
    assert(type(attempts) is int and attempts > 0)
    assert(not(nProc > 1 and method == 'download'))
    assert(method in ['search', 'download', 'opt'])

    if method in ['search', 'opt']:
        url = 'https://www.ebi.ac.uk/QuickGO/services/annotation/search'
        headers = {'Accept': 'application/json'}
    else:
        url = 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch'
        headers = {'Accept': 'text/tsv'}

    # set HTTP request body parameters
    defaults = {
        'taxonId': 9606,
        'geneProductType': 'protein',
        'geneProductSubset': 'Swiss-Prot',
        'proteome': 'gcrpCan',
        'goUsage': 'descendants',
        'goUsageRelationships': 'is_a,part_of,occurs_in',
        'includeFields': 'goName,synonyms',
        'limit': 100,
        'downloadLimit': 50000
    }
    params = {'goId': goIds}
    if useDefaults:
        params.update(defaults)
    params.update(kwargs)

    (r, attempts_remaining) = _get_QuickGO_helper(url, params, headers, method, attempts, sleep)
    if attempts_remaining > 0:
        allResponses = [r]
    else:
        return None

    if method == 'download':
        if parseResponse:
            df = pd.read_csv(io.StringIO(r.text), sep='\t')
            if 'DATE' in df.columns:
                df['DATE'] = pd.to_datetime(df['DATE'], yearfirst=True, format='%Y%m%d')
            return df
        else:
            return allResponses

    response_body = r.json()
    totalPages = response_body['pageInfo']['total']
    params.update({'page': response_body['pageInfo']['current']})
    if totalPages > 1 and params['page'] < totalPages:
        if method == 'opt' and response_body['numberOfHits'] < 50000:
            print('Switching to \'download\' API...')
            return get_QuickGO_annotations(goIds, useDefaults, 'download', parseResponse, **kwargs)

        params['page'] += 1
        if nProc > 1:
            allResponses.extend(_get_QuickGO_mt(url, params, headers, attempts, sleep, nProc,
                                                pages=range(params['page'],totalPages+1)))
        else:
            while True:
                (r, attempts_remaining) = _get_QuickGO_helper(url, params, headers, attempts, sleep)
                if attempts_remaining > 0:
                    allResponses.append(r)
                else:
                    print('Skipping page {}'.format(params['page']))
                if params['page'] >= totalPages: 
                    break
                params.update({'page': params['page'] + 1})

    if parseResponse:
        results = list(itertools.chain.from_iterable([r.json()['results'] for r in allResponses]))
        return parseQuickGOJSON(results)
    else:
        return allResponses

def _get_QuickGO_helper(url, params, headers, method, attempts=5, sleep=0.5):
    '''
    Returns: (request.Response, int)
      HTTP response and unused attempts
    '''

    while attempts > 0:
        try:
            r = requests.get(url, params=params, headers=headers)
            if not r.ok:
                r.raise_for_status()
            if method is not 'download':
                response_body = r.json()
                print('numberOfHits: {}'.format(response_body['numberOfHits']),
                      response_body['pageInfo'], sep=', ')
            break
        except Exception as err:
            print(err)
            if (attempts > 1):
                print('Attempts remaining: {}'.format(attempts-1))
                time.sleep(sleep)
            else:
                print('All attempts exhausted')
        attempts -= 1

    return (r, attempts)

def _get_QuickGO_mt(url, params, headers, attempts, sleep, nProc, pages):
    '''
    Multiprocess implementation of getQuickGO(method='search') without error handling
    '''

    with multiprocessing.Pool(nProc) as pool:
        print("Using {:d} processes...".format(pool._processes))
        sleep = max(sleep, 0.1*nProc)

        allResponses = []
        for i in pages:
            params.update({'page': i})
            allResponses.append(pool.apply_async(_get_QuickGO_helper,
                                                 (url, params, headers, 'search', attempts, sleep)))
        pool.close()
        pool.join()

    for i in range(len(allResponses)):
        allResponses[i] = allResponses[i].get()

    return allResponses

def get_QuickGO_descendants(goIds, relations=GO_RELATIONS, depth=-1, dset=None):
    '''
    Get descendant terms.

    Args
    - goIds: iterable of str
    - relations: iterable. default=GO_RELATIONS
        By default, includes all relations: is_a, part_of, occurs_in, regulates
    - depth: int. default=-1
        Number of layers of descendants to get (e.g., set depth=1 to retrieve only direct descendants).
        Set depth < 0 to return all descendants.
    - dset: set. default=set()
        Set of descendants already counted. Primarily used for recursive calls.

    Returns: set of str
    '''

    assert not isinstance(goIds, str) and isinstance(goIds, Iterable)

    if dset is None:
        dset = set()

    if goIds <= dset or depth == 0:
        return dset | goIds
    dset.update(goIds)

    # comma-separated GO IDs
    goIds_cs = ','.join(goIds)

    if depth < 0:
        requestURL = f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goIds_cs}/descendants'
        r = requests.get(requestURL, params={'relations': ','.join(relations)}, headers={'Accept': 'application/json'})
    else:
        requestURL = f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goIds_cs}/complete'
        r = requests.get(requestURL, headers={'Accept': 'application/json'})

    if not r.ok:
        r.raise_for_status()

    responseBody = r.json()

    if depth > 0:
        depth -= 1
        for i in range(responseBody['numberOfHits']):
            if 'children' not in responseBody['results'][i]:
                continue
            children = responseBody['results'][i]['children']
            descendants = set([children[j]['id'] for j in range(len(children)) if children[j]['relation'] in relations])
            new_descendants = descendants - dset
            dset.update(get_QuickGO_descendants(new_descendants, relations=relations, depth=depth, dset=dset))
    else:
        for i in range(responseBody['numberOfHits']):
            dset.update(responseBody['results'][i]['descendants'])

    return dset

def get_QuickGO_secondaryIds(goIds):
    '''
    Get secondary IDs.

    Args
    - goIds: iterable of str

    Returns: set of str
    '''

    goIds_cs = ','.join(goIds)
    requestURL = f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goIds_cs}/complete'
    r = requests.get(requestURL, headers={'Accept': 'application/json'})

    if not r.ok:
        r.raise_for_status()

    responseBody = r.json()

    secondaryIds = []
    for i in range(responseBody['numberOfHits']):
        if 'secondaryIds' in responseBody['results'][i]:
            secondaryIds.extend(responseBody['results'][i]['secondaryIds'])
    return set(secondaryIds)

# endregion --- QuickGO

# region ------ NCBI

def blastAndParse(program, database, sequence, format_type='XML', blast_fun=Bio.Blast.NCBIWWW.qblast,
                  log=sys.stdout, parse=None, save=None, compress=True, **kwargs):
    '''
    BLAST sequences and save results to file.

    Args
    - program: str.
        blastn, megablast, blastp, blastx, tblastn, or tblastx
    - database: str
        For DNA
          nt: Nucleotide collection
          refseq_rna: NCBI Transcript Reference Sequences
          pdbnt: PDB nucleotide database
        For Protein
          nr: Non-redundant
          refseq_protein: NCBI Protein Reference Sequences
          swissprot: Non-redundant UniProtKB/SwissProt sequences
          pdbaa: PDB protein database
    - sequence: str.
        FASTA, bare sequence, or newline-delimited identifiers
    - format_type: str. default='XML'
        HTML, Text, XML, XML2, JSON2, or Tabular. default=XML
    - blast_fun: function. default=Bio.Blast.NCBIWWW.qblast
        BLAST function. Should accept similar arguments as Bio.Blast.NCBIWWW.qblast and return a file object
        (e.g., io.StringIO).
    - log: str. default=sys.stdout
        Path to save log file.
    - save: str. default=None
        Path to save raw BLAST results
    - parse: str. default=None
        Directory in which to save parsed BLAST results (e.g., results from individual queries).
        Filenames are constructed from query IDs.
    - compress: bool. default=True
        Gzip-compress parsed files.
        Does not affect compression of `save` file, which is determined automatically from the supplied filename.
    - **kwargs: additional arguments to pass to blast_fun()

    Returns: io.StringIO or None
      Returns result of blast_fun if successful, otherwise None.

    Notes
    - program, database, sequence, and format_type are arguments directly passed to blast_fun()
    - Any (sub)directories required for the paths provided in save and parse will be automatically created.
    - Parsing of 'HTML', 'Text', 'XML2', and 'JSON2' format_type is not supported
    '''

    # Formats supported by BioPython's SearchIO module
    # Note that BLAST Text format_type is supported via 'blast-text', but writing is not
    parse_types = {'XML': 'blast-xml', 'Tabular': 'blast-tab'}
    extension = {'XML': '.xml', 'Text': '.txt', 'XML2': '.xml', 'JSON2': '.json', 'Tabular': '.tsv'}
    if compress:
        extension = {k: v + '.gz' for k, v in extension.items()}

    try:
        result = blast_fun(program, database, sequence, format_type=format_type, **kwargs)
    except Exception as err:
        print(err, file=log)
        return None

    if save is not None:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        f = utils_files.createFileObject(save, 'wt')
        f.write(result.read())
        f.close()
        result.seek(0)

    if parse is not None:
        if format_type in parse_types:
            os.makedirs(os.path.dirname(parse), exist_ok=True)
            for qresult in Bio.SearchIO.parse(result, format=parse_types[format_type]):
                filename = utils_files.createValidPath(qresult.id) + extension[format_type]
                with utils_files.createFileObject(os.path.join(parse, filename), 'wt') as f:
                    Bio.SearchIO.write(qresult, f, parse_types[format_type])
        else:
            print(f'Parsing of format_type {format_type} is not supported.', file=log)

    return result

def stop_blast(rids, url=Bio.Blast.NCBIWWW.NCBI_BLAST_URL):
    for rid in rids:
        requests.delete(url=url, params={'CMD': 'DELETE', 'RID': rid})
        requests.delete(url=url, data={'CMD': 'DELETE', 'RID': rid})
        requests.get(url=url, params={'CMD': 'DELETE', 'RID': rid})
        requests.get(url=url, data={'CMD': 'DELETE', 'RID': rid})
        requests.put(url=url, params={'CMD': 'DELETE', 'RID': rid})
        requests.put(url=url, data={'CMD': 'DELETE', 'RID': rid})

def my_qblast(program, database, sequence, url_base=Bio.Blast.NCBIWWW.NCBI_BLAST_URL,
              composition_based_statistics=None, entrez_query='(none)', matrix_name=None,
              expect=None, filter=None, gapcosts=None, genetic_code=None,
              hitlist_size=None, nucl_penalty=None, nucl_reward=None,
              threshold=None, word_size=None, alignments=None, descriptions=None,
              format_object=None, format_type=None, ncbi_gi=None,
              lock=None, queue=None, semaphore=None, initial_delay=10, server_delay=10, response_delay=60):
    '''
    Submit queries to NCBI BLAST server using the Common URL API. Supports multiprocessing
    synchronization while respecting server limits.

    BLAST Common URL API: https://ncbi.github.io/blast-cloud/dev/api.html

    Args: Modified from Bio.Blast.NCBIWWW.qblast with the following changes
    - Removed arguments not officially supported in the BLAST URL API (except entrez_query):
        auto_format, db_genetic_code, endpoints, genetic_code, i_thresh, other_advanced,
        perc_ident, phi_pattern, query_file, query_believe_defline, query_from,
        query_to, searchsp_eff, service, ungapped_alignment, alignment_view,
        entrez_links_new_window, expect_low, expect_high, format_entrez_query,
        template_type, template_length
    - Reset arguments to default values
        expect=None, hitlist_size=None, alignments=None, descriptions=None, format_type=None
    - Remap matrix_name to MATRIX
    - New arguments
      - Multiprocessing synchronization
        - lock: multiprocessing.Lock. default=None
            Once acquired, grants exclusive contact with server.
        - queue: multiprocessing.Queue. default=None
            Queue into which to pass (sequence, RID) tuple to communicate with parent process.
        - semaphore: multiprocessing.Semaphore. default=None
            Semaphore to limit number of parallel requests to the server.
      - Server limits based on https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        - server_delay: int. default=10
            Minimum number of seconds between any contact with the server
        - response_delay: int. default=60
            Minimum number of seconds between polling for result of the current RID

    Returns: io.StringIO, or requests.Response
      If BLAST is successful, returns io.StringIO of response text.
      If initial BLAST query submission is unsuccessful, returns the Response object.
    '''
    if lock is None:
        lock = multiprocessing.Lock()
    if semaphore is None:
        semaphore = multiprocessing.Semaphore()

    data = {
        'CMD': 'PUT',
        'QUERY': sequence,
        'DATABASE': database,
        'PROGRAM': program,
        'FILTER': filter,
        'EXPECT': expect,
        'NUCL_REWARD': nucl_reward,
        'NUCL_PENALTY': nucl_penalty,
        'GAPCOSTS': gapcosts,
        'MATRIX': matrix_name,
        'HITLIST_SIZE': hitlist_size,
        'THRESHOLD': threshold,
        'WORD_SIZE': word_size,
        'COMPOSITION_BASED_STATISTICS': composition_based_statistics,
        'ENTREZ_QUERY': entrez_query
    }
    data = {k: v for k, v in data.items() if v is not None}

    semaphore.acquire()
    lock.acquire()
    time.sleep(server_delay)
    r = requests.put(url=url_base, data=data)
    lock.release()
    if not r.ok:
        print('BLAST query submission error.', file=sys.stderr)
        semaphore.release()
        return r

    rid, _ = Bio.Blast.NCBIWWW._parse_qblast_ref_page(io.StringIO(r.text))
    if queue is not None:
        queue.put((sequence, rid))

    params = {
        'CMD': 'GET',
        'ALIGNMENTS': alignments,
        'DESCRIPTIONS': descriptions,
        'FORMAT_OBJECT': format_object,
        'FORMAT_TYPE': format_type,
        'NCBI_GI': ncbi_gi,
        'RID': rid
    }
    params = {k: v for k, v in params.items() if v is not None}

    previous = time.time()

    while True:
        # respect response_delay
        current = time.time()
        wait = previous + response_delay - current - server_delay
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        # respect server_delay
        lock.acquire()
        time.sleep(server_delay)
        r = requests.get(url=url_base, params=params)
        lock.release()

        results = r.text

        # Can see an "\n\n" page while results are in progress,
        if results == "\n\n":
            continue

        # XML results don't have the Status tag when finished
        if "Status=" not in results:
            break

        status = re.search(r'Status=(\S+)', results).group(1)
        if status.upper() in ['READY', 'FAILED', 'UNKNOWN', None]:
            if status.upper() == 'UNKNOWN':
                print(f'Search {rid} expired', file=sys.stderr)
            if status.upper() == 'FAILED':
                print(f'', file=sys.stdout)
            break
    semaphore.release()
    return io.StringIO(results)

def blastp_post(**kwargs):
    '''
    Submit query to NCBI BLAST server using a POST request.

    A POST request is what the web form (https://blast.ncbi.nlm.nih.gov/Blast.cgi) submits.
    An alternative is to use PUT-based requests, as implemented by Bio.Blast.NCBIWWW.qblast()
    in the BioPython package.

    References
    - https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
        Web form documentation

    General
    - QUERY: FASTA, bare sequence, or newline-delimited identifiers
    - db: nucleotide (blastn, tblastn), protein (blastp, blastx)
    - QUERY_FROM: Specify segment of the query sequence to BLAST
    - QUERY_TO: Specify segment of the query sequence to BLAST
    - QUERYFILE: Path to FASTA file
    - JOB_TITLE: Title of BLAST search
    - SUBJECTS: FASTA, bare sequence, or newline-delimited identifiers
    - stype: nucleotide (blastn, tblastn), protein (blastp, blastx)
    - SUBJECTS_FROM: Specify segment of the subject sequence to BLAST
    - SUBJECTS_TO: Specify segment of the subject sequence to BLAST
    - SUBJECTFILE: Path to FASTA file
    - EQ_MENU: Restrict search to specified organisms. Ex: 'Fungi (taxid:4751)'
    - ORG_EXCLUDE: Set to 'on' to exclude the organisms specified in EQ_MENU

    blastp
    - DATABASE: default=nr
        nr (Non-redundant protein sequences), refseq_protein (RefSeq), landmark (Model Organisms),
        swissprot (Swiss-Prot), pataa (Patented protein sequences), pdb (Protein Data Bank),
        env_nr (Metagenomic proteins), or tsa_nr (Transcriptome Shotgun Assembly proteins)
    - BLAST_PROGRAMS: default='blastp'
        kmerBlastp (Quick BLASTP), blastp, psiBlast (PSI-BLAST), phiBlast (PHI-BLAST), or
        deltaBlast (DELTA-BLAST)
    - MAX_NUM_SEQ: default=[blastp, kmerBlastp] 100; [psiBlast, phiBlast, deltaBlast] 500
        Maximum number of aligned sequences to display (web interface ranges from 100 to 20000)
    - SHORT_QUERY_ADJUST: default='on'
        Automatically adjust parameters for short input sequences
    - EXPECT: E-value threshold. default=10
    - WORD_SIZE: default=[blastp, kmerBlastp] 6; [psiBlast, phiBlast, deltaBlast] 3
        Length of the seed that initiates an alignment.
    - HSP_RANGE_MAX: default=0
        Limit the number of matches to a query range.
    - MATRIX_NAME: default=BLOSUM62
        PAM30, PAM70, PAM250, BLOSUM80, BLOSUM62, BLOSUM45, BLOSUM50, or BLOSUM90
    - GAPCOSTS: default='11 1'
        Cost to create and extend a gap in an alignment.
    - COMPOSITION_BASED_STATISTICS: default=2
        Matrix adjustment method to compensate for amino acid composition of sequences.
        0 (No adjustment), 1 (Composition-based statistics), 2 (Conditional compositional
        score matrix adjustment), 3 (Universal compositional score matrix adjustment)

    Unset
    - PHI_PATTERN: [phiBlast only] PHI pattern to start the search
    - PSSM: [psiBlast, phiBlast only] Path to PSSM file
    - I_THRESH: [psiBlast, phiBlast, deltaBlast only] default=0.005
        Statistical significance threshold to include a sequence in the model used by PSI-BLAST to create
        the PSSM on the next iteration. 
    - DI_THRESH: [deltaBlast only] default=0.05
        Statistical significance threshold to include a domain in the model used by DELTA-BLAST to create the PSSM 
    - PSI_PSEUDOCOUNT: [psiBlast, phiBlast, deltaBlast only] default=0
        Pseduocount parameter. If zero is specified, then the parameter is automatically determined through
        a minimum length description principle (PMID 19088134). A value of 30 is suggested in order to obtain
        the approximate behavior before the minimum length principle was implemented.
    - FILTER: include this argument once for each filter to set
        L (low complexity regions), M (mask for lookup table only)
    - LCASE_MASK: set to 'on' to enable mask
        Mask any letters that were lower-case in the FASTA input.

    Hidden fields for default formatting parameters
    - SHOW_OVERVIEW: on
    - SHOW_LINKOUT: on
    - GET_SEQUENCE: on
    - FORMAT_OBJECT: Alignment
    - FORMAT_TYPE: HTML
    - ALIGNMENT_VIEW: Pairwise
    - MASK_CHAR: 2
    - MASK_COLOR: 1
    - DESCRIPTIONS: 100
    - ALIGNMENTS: 100
    - LINE_LENGTH: 60
    - NEW_VIEW: on
    - OLD_VIEW: false
    - NCBI_GI: ''
    - SHOW_CDS_FEATURE: ''
    - NUM_OVERVIEW: 100
    - FORMAT_EQ_TEXT: ''
    - FORMAT_ORGANISM: ''
    - EXPECT_LOW: ''
    - EXPECT_HIGH: ''
    - PERC_IDENT_LOW: ''
    - PERC_IDENT_HIGH: ''
    - QUERY_INDEX: 0
    - FORMAT_NUM_ORG: 1
    - CONFIG_DESCR: 2,3,4,5,6,7,8

    Hidden fields
    - CLIENT: web
    - SERVICE: [blastp, kmerBlastp, psiBlast, phiBlast] plain; [deltaBlast] delta_blast
    - CMD: request
    - PAGE: Proteins
    - PROGRAM: blastp
    - MEGABLAST: ''
    - RUN_PSIBLAST: [blastp, kmerBlastp] ''; [psiBlast, phiBlast, deltaBlast] on
    - WWW_BLAST_TYPE: ''
    - CDD_SEARCH: on
    - ID_FOR_PSSM: ''
    - DB_DISPLAY_NAME: same as DATABASE
    - ORG_DBS: ''
    - SAVED_PSSM: ''
    - SELECTED_PROG_TYPE: same as DATABASE
    - SAVED_SEARCH: ''
    - BLAST_SPEC: ''
    - MIXED_DATABASE: ''
    - QUERY_BELIEVE_DEFLINE: ''
    - DB_DIR_PREFIX: ''
    - USER_DATABASE: ''
    - USER_WORD_SIZE: ''
    - USER_MATCH_SCORES: ''
    - USER_FORMAT_DEFAULTS: ''
    - NO_COMMON: ''
    - NUM_DIFFS: 
    - NUM_OPTS_DIFFS
    - UNIQ_DEFAULTS_NAME
    - PAGE_TYPE: BlastSearch
    '''

    default_kwargs = {
        'DATABASE': 'nr',
        'db': 'protein',
        'BLAST_PROGRAMS': 'blastp',
        'MAX_NUM_SEQ': 100,
        'WORD_SIZE': 3,
        'SHORT_QUERY_ADJUST': 'on',
        'EXPECT': 10,
        'HSP_RANGE_MAX': 0,
        'MATRIX_NAME': 'BLOSUM62',
        'GAPCOSTS': '11 1',
        'COMPOSITION_BASED_STATISTICS': 2
    }

    if 'BLAST_PROGRAMS' in kwargs and kwargs['BLAST_PROGRAMS'] in ['psiBlast', 'phiBlast', 'deltaBlast']:
        default_kwargs['MAX_NUM_SEQ'] = 500
        default_kwargs['WORD_SIZE'] = 3

    FILE_ARGS = ['QUERYFILE', 'SUBJECTFILE', 'PSSM']
    for key in kwargs:
        if key in FILE_ARGS:
            kwargs[key] = (kwargs[key], open(kwargs[key], 'rb'), 'application/octet-stream')
        else:
            kwargs[key] = (None, kwargs[key])

    r = requests.post(url=Bio.Blast.NCBIWWW.NCBI_BLAST_URL, files=data, allow_redirects=False)
    return r


# endregion --- NCBI