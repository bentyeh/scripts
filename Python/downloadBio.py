# add path of this file (e.g., a scripts directory) to sys.path
import sys, os
sys.path.append(os.path.dirname(__file__))

import io, re
from collections.abc import Iterable
from multiprocessing import Pool

import pandas as pd
import requests

import utils_files
import utils_bio

# region --- UniProt

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

# region --- InterPro

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

# region --- QuickGO

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
    
    with Pool(nProc) as pool:
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