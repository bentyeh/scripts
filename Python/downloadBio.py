'''
Services
- UniProt: www.uniprot.org
- InterPro: www.ebi.ac.uk/interpro
- QuickGO: www.ebi.ac.uk/QuickGO
- NCBI
- IDT: www.idtdna.com
- Kazusa DNA Research Institute Codon Usage Database: http://www.kazusa.or.jp/codon/
- uMelt: https://www.dna-utah.org/umelt/quartz/um.php
- STRING: www.string-db.org
'''

import collections, collections.abc, io, itertools, json, multiprocessing, os, re, time, sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import Bio.Blast.NCBIWWW
import Bio.SearchIO
from tqdm import tqdm

# add path of this file (e.g., a scripts directory) to sys.path
sys.path.append(Path(__file__).resolve(strict=True).parent)
import utils_files, utils_bio

# region ------ UniProt
# 
# Accession (AC): https://www.uniprot.org/help/accession_numbers
#   Example: P12345
# Entry names (formerly ID): https://www.uniprot.org/help/entry_name
#   Example: AATM_RABIT
#   Note: Entry names are not stable identifiers and are updated for consistency
#     (for instance to ensure that related entries have similar names or when a
#     UniProtKB/TrEMBL entry is integrated into UniProtKB/Swiss-Prot)
# 
# Proteins API Notes
# - Cross-reference databases to use with /proteins/{dbtype}:{dbid} service: https://www.uniprot.org/database/

UNIPROT_ACCESSION_REGEX = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
SWISSPROT_ENTRYNAME_REGEX = re.compile(r'[A-Z0-9]{1,5}_[A-Z0-9]{1,5}')
TREMBL_ENTRYNAME_REGEX = re.compile(UNIPROT_ACCESSION_REGEX.pattern + r'_[A-Z0-9]{1,5}')

PROTEINS_HEADERS = {
    'xml':   'application/xml',
    'json':  'application/json',
    'fasta': 'text/x-fasta',
    'flat':  'text/x-flatfile',
    'gff':   'text/x-gff',
    'peff':  'text/x-peff'
}

PROTEINS_BASEURL = 'https://www.ebi.ac.uk/proteins/api'

# Available columns when querying UniProt
# - Reference: https://www.uniprot.org/uniprot/#customize-columns
# - Usage
#     https://www.uniprot.org/uniprot/?
#       query=[ accession:<uniprot_ac>[+OR+[...]] | yourlist:<list_id> [filters]]
#       [&columns=<col1>[,<col2>[,...]]]
#       [&format=[ fasta | tab | xlsx | xml | rdf | txt | gff | list ]]
# - Example
#     https://www.uniprot.org/uniprot/?query=accession:P12345&columns=entry_name,genes,organism&format=tab
UNIPROT_COLUMNS = {
    'Entry': 'id',
    'Entry name': 'entry_name',
    'Gene names': 'genes',
    'Gene names (ordered locus)': 'genes_8_OLN_9_',
    'Gene names (ORF)': 'genes_8_ORF_9_',
    'Gene names (primary)': 'genes_8_PREFERRED_9_',
    'Gene names (synonym)': 'genes_8_ALTERNATIVE_9_',
    'Organism': 'organism',
    'Organism ID': 'organism-id',
    'Protein names': 'protein_names',
    'Proteomes': 'proteome',
    'Taxonomic lineage': 'lineage_8_ALL_9_',
    'Virus hosts': 'virus_hosts',
    'Alternative products (isoforms)': 'comment_8_ALTERNATIVE_PRODUCTS_9_',
    'Alternative sequence': 'feature_8_ALTERNATIVE_SEQUENCE_9_',
    'Erroneous gene model prediction': 'comment_8_ERRONEOUS_GENE_MODEL_PREDICTION_9_',
    'Fragment': 'fragment',
    'Gene encoded by': 'encodedon',
    'Length': 'length',
    'Mass': 'mass',
    'Mass spectrometry': 'comment_8_MASS_SPECTROMETRY_9_',
    'Natural variant': 'feature_8_NATURAL_VARIANT_9_',
    'Non-adjacent residues': 'feature_8_NON_ADJACENT_RESIDUES_9_',
    'Non-standard residue': 'feature_8_NON_STANDARD_RESIDUE_9_',
    'Non-terminal residue': 'feature_8_NON_TERMINAL_RESIDUE_9_',
    'Polymorphism': 'comment_8_POLYMORPHISM_9_',
    'RNA editing': 'comment_8_RNA_EDITING_9_',
    'Sequence': 'sequence',
    'Sequence caution': 'comment_8_SEQUENCE_CAUTION_9_',
    'Sequence conflict': 'feature_8_SEQUENCE_CONFLICT_9_',
    'Sequence uncertainty': 'feature_8_SEQUENCE_UNCERTAINTY_9_',
    'Sequence version': 'version_8_sequence_9_',
    'Absorption': 'comment_8_ABSORPTION_9_',
    'Active site': 'feature_8_ACTIVE_SITE_9_',
    'Activity regulation': 'comment_8_ACTIVITY_REGULATION_9_',
    'Binding site': 'feature_8_BINDING_SITE_9_',
    'Calcium binding': 'feature_8_CALCIUM_BIND_9_',
    'Catalytic activity': 'comment_8_CATALYTIC_ACTIVITY_9_',
    'Cofactor': 'comment_8_COFACTOR_9_',
    'DNA binding': 'feature_8_DNA_BINDING_9_',
    'EC number': 'ec',
    'Function': 'comment_8_FUNCTION_9_',
    'Kinetics': 'comment_8_KINETICS_9_',
    'Metal binding': 'feature_8_METAL_BINDING_9_',
    'Nucleotide binding': 'feature_8_NP_BIND_9_',
    'Pathway': 'comment_8_PATHWAY_9_',
    'pH dependence': 'comment_8_PH_DEPENDENCE_9_',
    'Redox potential': 'comment_8_REDOX_POTENTIAL_9_',
    'Rhea Ids': 'rhea-id',
    'Site': 'feature_8_SITE_9_',
    'Temperature dependence': 'comment_8_TEMPERATURE_DEPENDENCE_9_',
    'Annotation': 'annotation_score',
    'Caution': 'comment_8_CAUTION_9_',
    'Features': 'features',
    'Keyword ID': 'keyword-id',
    'Keywords': 'keywords',
    'Matched text': 'context',
    'Miscellaneous': 'comment_8_MISCELLANEOUS_9_',
    'Protein existence': 'existence',
    'Reviewed': 'reviewed',
    'Tools': 'tools',
    'UniParc': 'uniparcid',
    'Interacts with': 'interactor',
    'Subunit structure': 'comment_8_SUBUNIT_STRUCTURE_9_',
    'Developmental stage': 'comment_8_DEVELOPMENTAL_STAGE_9_',
    'Induction': 'comment_8_INDUCTION_9_',
    'Tissue specificity': 'comment_8_TISSUE_SPECIFICITY_9_',
    'Gene ontology (biological process)': 'go_8_biological_process_9_',
    'Gene ontology (cellular component)': 'go_8_cellular_component_9_',
    'Gene ontology (GO)': 'go',
    'Gene ontology (molecular function)': 'go_8_molecular_function_9_',
    'Gene ontology IDs': 'go-id',
    'ChEBI': 'chebi',
    'ChEBI (Catalytic activity)': 'chebi_8_Catalytic_activity_9_',
    'ChEBI (Cofactor)': 'chebi_8_Cofactor_9_',
    'ChEBI IDs': 'chebi-id',
    'Allergenic properties': 'comment_8_ALLERGENIC_PROPERTIES_9_',
    'Biotechnological use': 'comment_8_BIOTECHNOLOGICAL_USE_9_',
    'Disruption phenotype': 'comment_8_DISRUPTION_PHENOTYPE_9_',
    'Involvement in disease': 'comment_8_INVOLVEMENT_IN_DISEASE_9_',
    'Mutagenesis': 'feature_8_MUTAGENESIS_9_',
    'Pharmaceutical use': 'comment_8_PHARMACEUTICAL_USE_9_',
    'Toxic dose': 'comment_8_TOXIC_DOSE_9_',
    'Intramembrane': 'feature_8_INTRAMEMBRANE_9_',
    'Subcellular location': 'comment_8_SUBCELLULAR_LOCATION_9_',
    'Topological domain': 'feature_8_TOPOLOGICAL_DOMAIN_9_',
    'Transmembrane': 'feature_8_TRANSMEMBRANE_9_',
    'Chain': 'feature_8_CHAIN_9_',
    'Cross-link': 'feature_8_CROSS_LINK_9_',
    'Disulfide bond': 'feature_8_DISULFIDE_BOND_9_',
    'Glycosylation': 'feature_8_GLYCOSYLATION_9_',
    'Initiator methionine': 'feature_8_INITIATOR_METHIONINE_9_',
    'Lipidation': 'feature_8_LIPIDATION_9_',
    'Modified residue': 'feature_8_MODIFIED_RESIDUE_9_',
    'Peptide': 'feature_8_PEPTIDE_9_',
    'Post-translational modification': 'comment_8_POST-TRANSLATIONAL_MODIFICATION_9_',
    'Propeptide': 'feature_8_PROPEPTIDE_9_',
    'Signal peptide': 'feature_8_SIGNAL_9_',
    'Transit peptide': 'feature_8_TRANSIT_9_',
    '3D': '3d',
    'Beta strand': 'feature_8_BETA_STRAND_9_',
    'Helix': 'feature_8_HELIX_9_',
    'Turn': 'feature_8_TURN_9_',
    'Mapped PubMed ID': 'citationmapping',
    'PubMed ID': 'citation',
    'Date of creation': 'created',
    'Date of last modification': 'last-modified',
    'Date of last sequence modification': 'sequence-modified',
    'Entry version': 'version_8_entry_9_',
    'Coiled coil': 'feature_8_COILED_COIL_9_',
    'Compositional bias': 'feature_8_COMPOSITIONAL_BIAS_9_',
    'Domain': 'feature_8_DOMAIN_EXTENT_9_',
    'Motif': 'feature_8_MOTIF_9_',
    'Protein families': 'families',
    'Region': 'feature_8_REGION_9_',
    'Repeat': 'feature_8_REPEAT_9_',
    'Sequence similarities': 'comment_8_SEQUENCE_SIMILARITIES_9_',
    'Zinc finger': 'feature_8_ZINC_FINGER_9_',
    'Taxonomic lineage (all)': 'lineage_8_all_9_',
    'Taxonomic lineage (CLASS)': 'lineage_8_CLASS_9_',
    'Taxonomic lineage (COHORT)': 'lineage_8_COHORT_9_',
    'Taxonomic lineage (FAMILY)': 'lineage_8_FAMILY_9_',
    'Taxonomic lineage (FORMA)': 'lineage_8_FORMA_9_',
    'Taxonomic lineage (GENUS)': 'lineage_8_GENUS_9_',
    'Taxonomic lineage (INFRACLASS)': 'lineage_8_INFRACLASS_9_',
    'Taxonomic lineage (INFRAORDER)': 'lineage_8_INFRAORDER_9_',
    'Taxonomic lineage (KINGDOM)': 'lineage_8_KINGDOM_9_',
    'Taxonomic lineage (ORDER)': 'lineage_8_ORDER_9_',
    'Taxonomic lineage (PARVORDER)': 'lineage_8_PARVORDER_9_',
    'Taxonomic lineage (PHYLUM)': 'lineage_8_PHYLUM_9_',
    'Taxonomic lineage (SPECIES)': 'lineage_8_SPECIES_9_',
    'Taxonomic lineage (SPECIES_GROUP)': 'lineage_8_SPECIES_GROUP_9_',
    'Taxonomic lineage (SPECIES_SUBGROUP)': 'lineage_8_SPECIES_SUBGROUP_9_',
    'Taxonomic lineage (SUBCLASS)': 'lineage_8_SUBCLASS_9_',
    'Taxonomic lineage (SUBCOHORT)': 'lineage_8_SUBCOHORT_9_',
    'Taxonomic lineage (SUBFAMILY)': 'lineage_8_SUBFAMILY_9_',
    'Taxonomic lineage (SUBGENUS)': 'lineage_8_SUBGENUS_9_',
    'Taxonomic lineage (SUBKINGDOM)': 'lineage_8_SUBKINGDOM_9_',
    'Taxonomic lineage (SUBORDER)': 'lineage_8_SUBORDER_9_',
    'Taxonomic lineage (SUBPHYLUM)': 'lineage_8_SUBPHYLUM_9_',
    'Taxonomic lineage (SUBSPECIES)': 'lineage_8_SUBSPECIES_9_',
    'Taxonomic lineage (SUBTRIBE)': 'lineage_8_SUBTRIBE_9_',
    'Taxonomic lineage (SUPERCLASS)': 'lineage_8_SUPERCLASS_9_',
    'Taxonomic lineage (SUPERFAMILY)': 'lineage_8_SUPERFAMILY_9_',
    'Taxonomic lineage (SUPERKINGDOM)': 'lineage_8_SUPERKINGDOM_9_',
    'Taxonomic lineage (SUPERORDER)': 'lineage_8_SUPERORDER_9_',
    'Taxonomic lineage (SUPERPHYLUM)': 'lineage_8_SUPERPHYLUM_9_',
    'Taxonomic lineage (TRIBE)': 'lineage_8_TRIBE_9_',
    'Taxonomic lineage (VARIETAS)': 'lineage_8_VARIETAS_9_',
    'Taxonomic lineage IDs': 'lineage-id'
}

# Names of cross-reference databases
# - ID_TYPE: names used in UniProt's idmapping tables
#     ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
#     Note: PIR is found in idmapping_selected.tab, but not idmapping.tab
# - Abbreviation: abbreviations used in UniProt's Retrieve/ID mapping service
UNIPROT_IDTYPE_TO_ABBREVIATION = {
    'Allergome': 'ALLERGOME_ID',
    'ArachnoServer': 'ARACHNOSERVER_ID',
    'Araport': 'ARAPORT_ID',
    'BioCyc': 'BIOCYC_ID',
    'BioGrid': 'BIOGRID_ID',
    'BioMuta': 'BIOMUTA_ID',
    'CCDS': 'CCDS_ID',
    'CGD': 'CGD',
    'ChEMBL': 'CHEMBL_ID',
    'ChiTaRS': 'CHITARS_ID',
    'CollecTF': 'COLLECTF_ID',
    'ComplexPortal': 'COMPLEXPORTAL_ID',
    'ConoServer': 'CONOSERVER_ID',
    'CRC64': 'CRC64',
    'dictyBase': 'DICTYBASE_ID',
    'DIP': 'DIP_ID',
    'DisProt': 'DISPROT_ID',
    'DMDM': 'DMDM_ID',
    'DNASU': 'DNASU_ID',
    'DrugBank': 'DRUGBANK_ID',
    'EchoBASE': 'ECHOBASE_ID',
    'EcoGene': 'ECOGENE_ID',
    'eggNOG': 'EGGNOG_ID',
    'EMBL-CDS': 'EMBL',
    'EMBL': 'EMBL_ID',
    'Ensembl': 'ENSEMBL_ID',
    'EnsemblGenome': 'ENSEMBLGENOME_ID',
    'EnsemblGenome_PRO': 'ENSEMBLGENOME_PRO_ID',
    'EnsemblGenome_TRS': 'ENSEMBLGENOME_TRS_ID',
    'Ensembl_PRO': 'ENSEMBL_PRO_ID',
    'Ensembl_TRS': 'ENSEMBL_TRS_ID',
    'ESTHER': 'ESTHER_ID',
    'euHCVdb': 'EUHCVDB_ID',
    'EuPathDB': 'EUPATHDB_ID',
    'FlyBase': 'FLYBASE_ID',
    'GeneCards': 'GENECARDS_ID',
    'GeneDB': 'GENEDB_ID',
    'GeneID': 'P_ENTREZGENEID',
    'Gene_Name': 'GENENAME',
    'Gene_OrderedLocusName': None,
    'Gene_ORFName': None,
    'GeneReviews': 'GENEREVIEWS_ID',
    'Gene_Synonym': None,
    'GeneTree': 'GENETREE_ID',
    'GeneWiki': 'GENEWIKI_ID',
    'GenomeRNAi': 'GENOMERNAI_ID',
    'GI': 'P_GI',
    'GlyConnect': 'GLYCONNECT_ID',
    'GuidetoPHARMACOLOGY': 'GUIDETOPHARMACOLOGY_ID',
    'HGNC': 'HGNC_ID',
    'H-InvDB': 'H_INVDB_ID',
    'HOGENOM': 'HOGENOM_ID',
    'HPA': 'HPA_ID',
    'KEGG': 'KEGG_ID',
    'KO': 'KO_ID',
    'LegioList': 'LEGIOLIST_ID',
    'Leproma': 'LEPROMA_ID',
    'MaizeGDB': 'MAIZEGDB_ID',
    'MEROPS': 'MEROPS_ID',
    'MGI': 'MGI_ID',
    'MIM': 'MIM_ID',
    'MINT': None,
    'mycoCLAP': 'MYCOCLAP_ID',
    'NCBI_TaxID': None,
    'neXtProt': 'NEXTPROT_ID',
    'OMA': 'OMA_ID',
    'Orphanet': 'ORPHANET_ID',
    'OrthoDB': 'ORTHODB_ID',
    'PATRIC': 'PATRIC_ID',
    'PDB': 'PDB_ID',
    'PeroxiBase': 'PEROXIBASE_ID',
    'PharmGKB': 'PHARMGKB_ID',
    'PIR': 'PIR',
    'PomBase': 'POMBASE_ID',
    'ProteomicsDB': 'PROTEOMICSDB_ID',
    'PseudoCAP': 'PSEUDOCAP_ID',
    'Reactome': 'REACTOME_ID',
    'REBASE': 'REBASE_ID',
    'RefSeq': 'P_REFSEQ_AC',
    'RefSeq_NT': 'REFSEQ_NT_ID',
    'RGD': 'RGD_ID',
    'SGD': 'SGD_ID',
    'STRING': 'STRING_ID',
    'SwissLipids': 'SWISSLIPIDS_ID',
    'TAIR': None,
    'TCDB': 'TCDB_ID',
    'TreeFam': 'TREEFAM_ID',
    'TubercuList': 'TUBERCULIST_ID',
    'UCSC': 'UCSC_ID',
    'UniParc': 'UPARC',
    'UniPathway': 'UNIPATHWAY_ID',
    'UniProtKB-ID': 'ID',
    'UniRef100': 'NF100',
    'UniRef50': 'NF50',
    'UniRef90': 'NF90',
    'VectorBase': 'VECTORBASE_ID',
    'VGNC': 'VGNC_ID',
    'WBParaSite': 'WBPARASITE_ID',
    'World-2DPAGE': 'WORLD_2DPAGE_ID',
    'WormBase': 'WORMBASE_ID',
    'WormBase_PRO': 'WORMBASE_PRO_ID',
    'WormBase_TRS': 'WORMBASE_TRS_ID',
    'Xenbase': 'XENBASE_ID',
    'ZFIN': 'ZFIN_ID'
}

UNIPROT_IDMAPPING_HEADER = ['UniProtKB-AC', 'ID_type', 'ID']

UNIPROT_DB_NAMES = ("UniProtKB", "UniProtKB-Swiss-Prot", "UniRef50", "UniRef90", "UniRef100", "UniParc", "UniProtKB_AC-ID")

def convert_via_UniProt(source, to, query, taxId=None, intermediate="UniProtKB-Swiss-Prot", poll_interval=2, max_time=600, verbose=False):
    '''
    Wrapper around convert_UniProt to convert between any 2 databases.

    The only new parameter is `intermediate`, which must be either "UniProtKB" or "UniProtKB-Swiss-Prot".
    - This specifies which intermediate database (all of UniProt, including UniProtKB/Swiss-Prot and UniProtKB/TrEMBL;
      or UniProtKB/Swiss-Prot only) to use for conversion.
    '''
    assert intermediate in ("UniProtKB", "UniProtKB-Swiss-Prot")
    if source == "UniProtKB_AC-ID" or to in ("UniProtKB", "UniProtKB-Swiss-Prot"):
        return convert_UniProt(source, to, query, taxId, poll_interval, max_time, verbose)
    tmp = convert_UniProt(source, intermediate, query, taxId, poll_interval, max_time, verbose)
    if len(tmp) == 1:
        # job did not finish before max_time
        return tmp
    else:
        results, failedIds = tmp
    failedIds = set(failedIds)
    new_query = set()
    for s in results.values():
        new_query |= s
    tmp = convert_UniProt("UniProtKB_AC-ID", to, list(new_query), taxId, poll_interval, max_time, verbose)
    if len(tmp) == 1:
        # job did not finish before max_time
        return tmp
    else:
        results2, failedIds2 = tmp
    results_final = dict()
    for id_orig, uniprot_ids in results.items():
        failed = True
        for uniprot_id in uniprot_ids:
            if uniprot_id in results2:
                failed = False
                if id_orig not in results_final:
                    results_final[id_orig] = results2[uniprot_id]
                else:
                    results_final[id_orig] |= results2[uniprot_id]
        if failed:
            failedIds.add(id_orig)
    assert len(failedIds & set(results_final.keys())) == 0
    reconstructed_queries = failedIds | set(results_final.keys())
    missing_queries = set(query) - reconstructed_queries
    if len(missing_queries) > 0:
        print('Failed to get results for', missing_queries, file=sys.stderr)
    extra_queries = reconstructed_queries - set(query)
    if len(extra_queries) > 0:
        print('Results for unsolicited queries', extra_queries, file=sys.stderr)
    # assert all(entry in reconstructed_queries for entry in query)
    return results_final, list(failedIds)

def convert_UniProt(source, to, query, taxId=None, poll_interval=2, max_time=600, verbose=False):
    '''
    Convert database identifiers using UniProt's ID mapping service.

    See https://rest.uniprot.org/configure/idmapping/fields for names of databases
    to use as arguments for "source" and "to". In general, a single server request
    can convert between any non-UniProt identifier and a UniProt identifier, but cannot
    convert between two non-UniProt identifiers. Use convert_via_UniProt() to automatically
    first convert a non-UniProt identifier to a UniProt identifier, and then convert to
    the target non-UniProt identifier.

    Args
    - source: str
        Query identifier type. (e.g., UniProtKB, Ensembl_Protein, etc.)
    - to: str
        Target identifier type. (e.g., UniProtKB, Gene_Name, etc.)
    - query: list of str, or io.TextIOBase
        List of identifiers
    - taxId: str. default=None
        NCBI taxonomy ID
    - poll_interval: int or float. default=2
        Seconds to sleep between polling the server for job status
    - max_time: int or float. default=600 (10 minutes)
        Maximum seconds to wait before returning.
    - verbose: bool. default=False
        Print messages about jobID, etc.

    Returns
    - If max_time is not exceeded:
      - suceeded: dict[str, set[str]] mapping query identifiers to target identifiers
      - failed: list[str] of query identifiers that could not be successfully mapped
    - If max_time is exceeded:
      - jobId: str

    Reference: https://www.uniprot.org/help/id_mapping
    '''
    url_mapping = 'https://rest.uniprot.org/idmapping'
    url_submit = f'{url_mapping}/run'
    url_result = f'{url_mapping}/results'

    r = requests.post(
        url=url_submit,
        data={
            'from': source,
            'to': to,
            'ids': ','.join(query),
            'taxId': taxId
        }
    )
    if not r.ok and source in UNIPROT_DB_NAMES and taxId is not None:
        print('Since source is a UniProt database, re-trying without taxonomy ID.', file=sys.stderr)
        r = requests.post(
            url=url_submit,
            data={
                'from': source,
                'to': to,
                'ids': ','.join(query)
            }
        )
    if not r.ok:
        r.raise_for_status()

    jobId = r.json()['jobId']
    if verbose:
        print('UniProt ID mapping request URL:', r.url)
        print('UniProt ID mapping job ID:', jobId)

    start = time.time()
    current_time = time.time()
    finished = False
    while current_time - start < max_time:
        status = requests.get(f"{url_mapping}/status/{jobId}")
        if not status.ok:
            status.raise_for_status()
        j = status.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                if verbose:
                    print(f"UniProt ID mapping waiting for job to complete: retrying in {poll_interval}s")
                time.sleep(poll_interval)
                current_time = time.time()
            else:
                raise Exception(j["jobStatus"])
        else:
            finished = True
            break

    if not finished:
        if verbose:
            print(f"UniProt ID mapping job did not finish within max_time of {max_time}s. Returning jobId.")
        return jobId

    if "failedIds" in j:
        failedIds = j["failedIds"]
    else:
        failedIds = []

    results = dict()
    for d in j['results']:
        k = d['from']
        v = d['to']
        # "If you map to UniProtKB, UniParc or UniRef data, the full entries will be returned to you for convenience."
        if to in ('UniProtKB', "UniProtKB-Swiss-Prot"):
            v = v['primaryAccession']
        elif to in ("UniRef50", "UniRef90", "UniRef100"):
            v = v['id']
        elif to == 'UniParc':
            v = v['uniParcId']
        if k not in results:
            results[k] = {v}
        else:
            results[k].add(v)
    while 'Link' in status.headers:
        url_next = re.match(r'<(.+)>; rel="next"', status.headers["Link"]).group(1)
        status = requests.get(url_next)
        j = status.json()
        for d in j['results']:
            k = d['from']
            v = d['to']
            if to in ('UniProtKB', "UniProtKB-Swiss-Prot"):
                v = v['primaryAccession']
            elif to in ("UniRef50", "UniRef90", "UniRef100"):
                v = v['id']
            elif to == 'UniParc':
                v = v['uniParcId']
            if k not in results:
                results[k] = {v}
            else:
                results[k].add(v)

    assert len(set(failedIds) & set(results.keys())) == 0
    reconstructed_queries = set(failedIds) | set(results.keys())
    missing_queries = set(query) - reconstructed_queries
    if len(missing_queries) > 0:
        print('Failed to get results for', missing_queries, file=sys.stderr)
    extra_queries = reconstructed_queries - set(query)
    if len(extra_queries) > 0:
        print('Results for unsolicited queries', extra_queries, file=sys.stderr)
    # assert all(entry in reconstructed_queries for entry in query)
    return results, failedIds

def get_UniProt(uniprot_ac, export, toDf=False):
    '''
    Get UniProt entry.

    Args
    - uniprot_ac: str
        UniProt Accession (e.g., P12345)
    - export: str.
        Format: 'txt', 'fasta', 'xml', 'rdf', 'gff'
    - toDf: bool. default=False
        Parse data into a pandas.DataFrame. Currently only 'fasta' and 'gff' export types are supported.

    Returns: str or pandas.DataFrame
    '''
    assert export in ('txt', 'fasta', 'xml', 'rdf', 'gff')
    url = f'https://www.uniprot.org/uniprot/{uniprot_ac}.{export}'
    headers = {'accept': 'text/html'}
    r = requests.get(url=url, headers=headers)
    if not r.ok:
        r.raise_for_status()

    if toDf:
        if export == 'fasta':
            return utils_bio.fastaToDf(io.StringIO(r.text), headerParser=utils_bio.parseUniProtHeader)
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
        Required parameter(s) specified with brackets '{param}' in `service`.
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
    - Single protein entry, FASTA format
        GET request: https://www.ebi.ac.uk/proteins/api/proteins/P12345
        --> get_Proteins(service='/proteins/{accession}', entry={'accession': 'P12345'}, export='fasta')
    - Multiple protein entries, JSON format
        GET request: https://www.ebi.ac.uk/proteins/api/proteins?accession=P12345,P67890
        --> get_Proteins(service='/proteins', entry=None, export='json', accession='P12345,P67890')
    - Single protein entry by cross-reference, JSON format
        GET request: https://www.ebi.ac.uk/proteins/api/proteins/EMBL:AAGW02052878
        --> get_Proteins(service='/proteins/{dbtype}:{dbid}', entry={'dbtype': 'EMBL', 'dbid': 'AAGW02052878'})
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

def get_Proteins_taxa(uniprot_acs, rank=None, cache=None, verbose=True, **kwargs):
    '''
    Lookup parent taxa of a protein at a given rank.

    Args
    - uniprot_ac: list of str
        UniProt Accessions (e.g., P12345)
    - rank: str. default=None
        Taxonomic rank. See NCBI_TAXONOMY_RANKS.
        If None, returns the taxonomy of the protein.
    - cache: dict of str -> int. default=None
        Map of taxa of input proteins to parent taxa at given rank.
    - **kwargs: additional arguments to pass to get_taxon_at_rank()
        Only applicable if rank is not None.

    Returns: dict of str -> int[, dict of int -> int]
      Map of protein to parent taxa id at a given rank.
      [If rank is not None] Map of taxa id of input proteins to parent taxa id at given rank.
    '''
    max_query_size = 100
    uniprot_acs = list(set(uniprot_acs)) # remove duplicates
    taxa = {}
    for i in tqdm(range(int(np.ceil(len(uniprot_acs)/max_query_size))), disable=not verbose):
        proteins = uniprot_acs[i*max_query_size:(i+1)*max_query_size]
        results = get_Proteins(service='/proteins', entry=None, export='json', accession=','.join(proteins))
        results = json.loads(results)
        for protein in results:
            taxa[protein['accession']] = protein['organism']['taxonomy']

    if rank is None:
        return taxa

    if cache is None:
        cache = {}

    taxa_at_rank = get_taxon_at_rank(list(set(taxa.values()) - set(cache.keys())), rank=rank, **kwargs)
    taxa_at_rank.update(cache)
    return ({acc: taxa_at_rank[taxa[acc]] for acc in taxa}, taxa_at_rank)

def get_taxon_at_rank(taxa, rank, check_self=True, nproc=1, verbose=True):
    '''
    Lookup parent taxa at a given rank.

    Args
    - taxa: list of int
        Taxa to lookup
    - rank: str
        Taxonomic rank. See NCBI_TAXONOMY_RANKS.
    - check_self: bool. default=True
        Before looking up full lineage, check if input taxa are already at desired rank.
    - nproc: int. default=1
        Number of processes to use. The EMBL server appears fairly robust to >10 simultaneous requests.
    - verbose: bool. default=True
        Show progress bar.

    Returns: dict of int -> int
      Map input taxa to parent taxa at a given rank. If no parent taxa at given
      rank exists (e.g., the rank provided is below the rank of the input taxa),
      the value of such input taxa will be None.
    '''
    d = {taxon: None for taxon in taxa}
    taxa = list(set(taxa)) # remove duplicates

    if verbose:
        callback = lambda *x: pbar.update()
    else:
        callback = None

    if check_self:
        with tqdm(total=len(taxa), desc='check_self', disable=not verbose) as pbar:
            max_query_size = 50
            results = []
            with multiprocessing.Pool(max(nproc, 10)) as pool:
                for i in range(int(np.ceil(len(taxa)/max_query_size))):
                    query = ','.join([str(i) for i in taxa[i*max_query_size:(i+1)*max_query_size]])
                    results.append(pool.apply_async(
                        get_Proteins,
                        args=tuple(),
                        kwds=dict(service='/taxonomy/ids/{ids}/node', entry={'ids': query}, export='json'),
                        callback=callback
                    ))
                pool.close()
                pool.join()
        for result in results:
            for taxon in json.loads(result.get())['taxonomies']:
                if taxon.get('rank') == rank:
                    d[taxon['taxonomyId']] = taxon['taxonomyId']

    results = []
    taxa = set(taxa) - set([k for k in d.keys() if d[k] is not None])
    with tqdm(total=len(taxa), desc='lookup_lineage', disable=not verbose or len(taxa) == 0) as pbar:
        with multiprocessing.Pool(nproc) as pool:
            for taxon in taxa:
                results.append(pool.apply_async(
                    get_Proteins,
                    args=tuple(),
                    kwds=dict(service='/taxonomy/lineage/{id}', entry={'id': taxon}, export='json'),
                    callback=callback
                ))
            pool.close()
            pool.join()

    lineages = [result.get() for result in results]
    for lineage in lineages:
        nodes = json.loads(lineage)['taxonomies']
        for node in nodes:
            if node['rank'] == rank:
                d[nodes[0]['taxonomyId']] = node['taxonomyId']
                break
    return d

# endregion --- UniProt

# region ------ InterPro

INTERPRO_PROTEIN_URL = 'https://www.ebi.ac.uk/interpro/protein/'

INTERPRO_GOMAP_REGEX = re.compile(
    r'InterPro:(?P<InterPro_ID>[^\s]+)\s+' +
    r'(?P<InterPro_desc>[^>]+)\s+[>]\s+' + 
    r'GO:(?P<GO_desc>[^;]+)\s+[;]\s+' +
    r'(?P<GO_ID>GO[:]\d+)'
)

def get_InterPro_protein(uniprot_ac, export='tsv', toDf=True):
    '''
    Get InterPro protein sequence or entry annotation table.

    Args
    - uniprot_ac: str
        UniProt Accession (e.g., P12345)
    - export: str. default='tsv'
        'tsv' or 'fasta'
    - toDF: bool. default=True
        Parse data into a pandas.DataFrame. Otherwise returns the raw text.
        - export=='tsv': parse TSV format file directly into data frame
        - export=='fasta': parse FASTA to pandas.DataFrame using utils_bio.fastaToDf()

    Returns: str or pandas.DataFrame
    '''

    assert export in ('tsv', 'fasta')

    url = INTERPRO_PROTEIN_URL + uniprot_ac
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
            return utils_bio.fastaToDf(io.StringIO(r.text), headerParser=utils_bio.parseDefaultHeader)
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
        m = INTERPRO_GOMAP_REGEX.match(line)
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

    colNames_map = collections.OrderedDict([
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
                    withFrom_list.append(':'.join([results[i]['withFrom'][j]['connectedXrefs'][k]['db'],
                                                   results[i]['withFrom'][j]['connectedXrefs'][k]['id']]))
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
                            attempts=5, sleep=0.5, nproc=1, **kwargs):
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
    - nproc: int. default=1
        Number of processes to use. Error handling is not implemented if nproc > 1.
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
    assert(not(nproc > 1 and method == 'download'))
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
        if nproc > 1:
            allResponses.extend(_get_QuickGO_mt(url, params, headers, attempts, sleep, nproc,
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
        return parse_QuickGO_JSON(results)
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
            if method != 'download':
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

def _get_QuickGO_mt(url, params, headers, attempts, sleep, nproc, pages):
    '''
    Multiprocess implementation of getQuickGO(method='search') without error handling
    '''

    with multiprocessing.Pool(nproc) as pool:
        print("Using {:d} processes...".format(pool._processes))
        sleep = max(sleep, 0.1*nproc)

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

    assert not isinstance(goIds, str) and isinstance(goIds, collections.abc.Iterable)

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
# 
# Accession numbers
# - Sequence ID format: https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta
#   - Also see https://blast.advbiocomp.com/doc/FAQ-Indexing.html and 
#     https://en.wikipedia.org/wiki/FASTA_format#NCBI_identifiers
# - Database-specific formats
#   - GenBank: https://www.ncbi.nlm.nih.gov/Sequin/acc.html
#   - RefSeq: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/

NCBI_IDTYPE = {
    'lcl': 'local',
    'bbs': 'GenInfo backbone seqid',
    'bbm': 'GenInfo backbone moltype',
    'gim': 'GenInfo import ID',
    'gb' : 'GenBank',
    'emb': 'EMBL',
    'pir': 'PIR',
    'sp' : 'SWISS-PROT',
    'pat': 'patent',
    'pgp': 'pre-grant patent',
    'ref': 'RefSeq',
    'gnl': 'general database reference',
    'gi' : 'GenInfo integrated database',
    'dbj': 'DDBJ',
    'prf': 'PRF',
    'pdb': 'PDB',
    'tpg': 'third-party GenBank',
    'tpe': 'third-party EMBL',
    'tpd': 'third-party DDBJ',
    'tr' : 'TrEMBL',
    'gpp': 'genome pipeline',
    'nat': 'named annotation track'
}

NCBI_IDTYPE_TO_UNIPROT_IDTYPE = {
    'gb' : 'EMBL-CDS',
    'emb': 'EMBL-CDS',
    'ref': 'RefSeq',
    'gi' : 'GI',
    'ddj': 'EMBL-CDS',
    'pdb': 'PDB',
    'pir': 'PIR',
    'tgp': 'EMBL-CDS',
    'tpe': 'EMBL-CDS',
    'tpd': 'EMBL-CDS'
}

NCBI_BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

NCBI_TAXONOMY_RANKS = [
    'no rank', 'superkingdom', 'genus', 'species', 'order', 'family', 'subspecies', 'subfamily',
    'tribe', 'phylum', 'class', 'species group', 'forma', 'suborder', 'subclass', 'varietas',
    'kingdom', 'subphylum', 'infraorder', 'superfamily', 'infraclass', 'superorder', 'subgenus',
    'superclass', 'parvorder', 'superphylum', 'species subgroup', 'subcohort', 'cohort', 'subtribe',
    'section', 'series', 'subkingdom', 'subsection']

def blast_delete(rids, url_base=NCBI_BLAST_URL, verbose=True, **kwargs):
    '''
    Delete BLAST searches by Request ID (RID).

    Args
    - rids: list of str
        Request IDs
    - url_base: str. default=NCBI_BLAST_URL
        BLAST server URL
    - verbose: bool. default=True
        Print out RID that is being deleted
    - **kwargs: arguments to pass to requests.delete()

    Returns: None

    Notes
    - This function currently does not seem to have any effect. Results are still
      retrievable from the NCBI server via blast_get() or going to
      https://blast.ncbi.nlm.nih.gov/Blast.cgi?RID=<rid>&CMD=Get
    - Not sure if using `params` or `data` argument in requests.delete() call is correct.
    '''
    for rid in rids:
        if verbose:
            print(f'Stopping {rid}')
        requests.delete(url=url_base, params={'CMD': 'DELETE', 'RID': rid}, **kwargs)
    return None

def blast_submit(program, database, query=None, sequence=None, url_base=NCBI_BLAST_URL,
                 composition_based_statistics=None, entrez_query='(none)', matrix_name=None,
                 matrix=None, expect=None, filter=None, gapcosts=None, hitlist_size=None,
                 nucl_penalty=None, nucl_reward=None, threshold=None, word_size=None,
                 ncbi_gi=None, lock=None, server_delay=10):
    '''
    Submit BLAST query and return Request ID (RID).

    BLAST Common URL API Parameters: https://ncbi.github.io/blast-cloud/dev/api.html
    - program: str
        One of blastn, blastp, blastx, tblastn, tblastx. To enable megablast, use PROGRAM=blastn&MEGABLAST=on.
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
    - query: str. default=None
        FASTA, bare sequence, or newline-delimited identifiers
        Aliased by `sequence` (used by Biopython's qblast()). `query` takes precedence if both arguments are given.
    - composition_based_statistics: int. default=None
        Matrix adjustment method to compensate for amino acid composition of sequences.
          0 (No adjustment), 1 (Composition-based statistics), 2 (Conditional compositional
          score matrix adjustment), 3 (Universal compositional score matrix adjustment)
        If None, appears to default to 2.
    - entrez_query: str. default='(none)'
        Note: This is unofficially supported by the NCBI BLAST server.
    - matrix: str. default=None
        Aliased by `matrix_name` (used by Biopython's qblast()). `matrix` takes precedence if both arguments are give.
        If None, defaults to 'BLOSUM62'.
    - expect: float. default=None.
        Expect value.
        If None, appears to default to 10.
    - filter: str. default=None.
        Low complexity filtering. 'F' to disable. 'T' or 'L' to enable. Prepend 'm' for mask at lookup (e.g., 'mL')
        If None, appears to default to 'F'.
    - gapcosts: str. default=None
        Gap existence and extension costs. Pair of positive integers separated by a space such as '11 1'.
        If None, appears to default to '11 1'.
    - hitlist_size: int. default=None
        Number of databases sequences to keep
        If None, appears to default to 100.
    - nucl_penalty: int. default=None
        Cost for mismatched bases (BLASTN and megaBLAST)
    - nucl_reward: int. default=None
        Reward for matching bases (BLASTN and megaBLAST)
    - threshold: int. default=None
        Neighboring score for initial words. Positive integer (BLASTP default is 11). Does not apply to BLASTN or MegaBLAST.
        If None, the value depends on word_size.
          blastp - (word_size, threshold): (3, 11), (6, 21)
    - word_size: int. default=None
        Size of word for initial matches.
        If None, appears to default to 3.
    - ncbi_gi: str. default=None
        'T' or 'F'. Show NCBI GIs in report.
        Does not appear to have any effect. Empirically, GIs are included only if entrez_query is not '(none)'

    Args
    - url_base: str. default=NCBI_BLAST_URL
        BLAST server URL
    - lock: multiprocessing.Lock. default=None
        Once acquired, grants exclusive contact with server. Any server delays are applied before
        releasing the lock.
    - server_delay: int. default=10
        Minimum number of seconds between any contact with the server

    Returns: str or requests.Response
      Request ID (RID). Returns the request response if unable to find RID in response (e.g., unsuccessful submission).

    Notes
    - Officially documented BLAST parameters defaults are stated as defaults. Empirically observed defaults
      are stated as 'appears to default to'.
    '''
    assert (query is not None) or (sequence is not None)
    if query is None:
        query = sequence
    if matrix is None:
        matrix = matrix_name
    if lock is None:
        lock = multiprocessing.Lock()

    data = {
        'CMD': 'PUT',
        'QUERY': query,
        'DATABASE': database,
        'PROGRAM': program,
        'FILTER': filter,
        'EXPECT': expect,
        'NUCL_REWARD': nucl_reward,
        'NUCL_PENALTY': nucl_penalty,
        'GAPCOSTS': gapcosts,
        'MATRIX': matrix,
        'HITLIST_SIZE': hitlist_size,
        'THRESHOLD': threshold,
        'WORD_SIZE': word_size,
        'COMPOSITION_BASED_STATISTICS': composition_based_statistics,
        'NCBI_GI': ncbi_gi,
        'ENTREZ_QUERY': entrez_query
    }
    data = {k: v for k, v in data.items() if v is not None}

    lock.acquire()
    r = requests.put(url=url_base, data=data)
    time.sleep(server_delay)
    lock.release()
    if not r.ok:
        return r

    try:
        rid, _ = Bio.Blast.NCBIWWW._parse_qblast_ref_page(io.StringIO(r.text))
        return rid
    except ValueError as err:
        print(err, file=sys.stderr)
        return r

def blast_get(rid, url_base=NCBI_BLAST_URL, no_wait=True,
              alignments=None, descriptions=None, format_object=None, format_type=None, ncbi_gi=None, hitlist_size=None,
              lock=None, initial_delay=0, server_delay=10, response_delay=60):
    '''
    Retrieve BLAST result by Request ID (RID).

    BLAST Common URL API Args: https://ncbi.github.io/blast-cloud/dev/api.html
    - alignments: int. default=None
        Number of alignments to print (applies to HTML and Text format_type)
    - descriptions: int. default=None
        Number of descriptions to print (applies to HTML and Text format_type)
    - format_object: str. default=None
        'SearchInfo': status check
        'Alignment': report formatting
    - format_type: str. default=None
        HTML (default if None), Text, XML, XML2, JSON2, or Tabular
    - ncbi_gi: str. default=None
        'T' or 'F'. Show NCBI GIs in report.
        Does not appear to have any effect. Empirically, GIs are included only if entrez_query was not '(none)'
        in the original BLAST query submission.
    - hitlist_size: int. default=None
        Number of databases sequences to keep
        If None, appears to default to 100.

    Args
    - rid: str
        BLAST Request ID
    - url_base: str. default=NCBI_BLAST_URL
        BLAST server URL
    - no_wait: bool. default=True
        Return request result immediately.
    - lock: multiprocessing.Lock. default=None
        Once acquired, grants exclusive contact with server. Any server delays are applied before
        releasing the lock.
    - initial_delay: int. default=0
        Minimum number of seconds before initial polling request of the current RID
    - server_delay: int. default=10
        Minimum number of seconds between any contact with the server
    - response_delay: int. default=60
        Minimum number of seconds between polling for result of the current RID

    Returns: io.StringIO
      io.StringIO of the request response
    '''
    if lock is None:
        lock = multiprocessing.Lock()

    params = {
        'CMD': 'GET',
        'ALIGNMENTS': alignments,
        'DESCRIPTIONS': descriptions,
        'HITLIST_SIZE': hitlist_size,
        'FORMAT_OBJECT': format_object,
        'FORMAT_TYPE': format_type,
        'NCBI_GI': ncbi_gi,
        'RID': rid
    }
    params = {k: v for k, v in params.items() if v is not None}

    time.sleep(initial_delay)
    while True:
        # respect server_delay
        lock.acquire()
        r = requests.get(url=url_base, params=params)
        previous = time.time()
        time.sleep(server_delay)
        lock.release()

        results = r.text

        # Can see an "\n\n" page while results are in progress
        if results == "\n\n":
            continue

        # XML results don't have the Status tag when finished
        if "Status=" not in results:
            break

        status = re.search(r'Status=(\S+)', results).group(1)
        if status.upper() in ['READY', 'FAILED', 'UNKNOWN', None] or no_wait:
            if status.upper() == 'UNKNOWN':
                print(f'Search {rid} expired', file=sys.stderr)
            if status.upper() == 'FAILED':
                print(f'Search {rid} failed', file=sys.stderr)
            break

        # respect response_delay
        current = time.time()
        wait = max(previous + response_delay - current, 0)
        time.sleep(wait)
        previous = current + wait

    return io.StringIO(results)

def my_qblast(program, database, query=None, sequence=None, semaphore=None, queue=None, job_title=None, **kwargs):
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
    - New arguments
      - job_title: str
      - Aliases to match argument names of the BLAST Common URL API
        - matrix: higher precedence alias of matrix_name
        - query: higher precedence alias of sequence
      - Multiprocessing synchronization
        - lock: multiprocessing.Lock. default=None
            Once acquired, grants exclusive contact with server. Any server delays are applied before
            releasing the lock.
        - queue: multiprocessing.Queue. default=None
            Queue into which to pass (job_title, RID) or (sequence, RID) tuple to communicate with parent process.
        - semaphore: multiprocessing.Semaphore. default=None
            Semaphore to limit number of parallel requests to the server.
      - Server limits
        - initial_delay: int. default=None
            Minimum number of seconds before initial polling request of the current RID
        - server_delay: int. default=None
            Minimum number of seconds between any contact with the server
        - response_delay: int. default=None
            Minimum number of seconds between polling for result of the current RID

    Returns: io.StringIO, or requests.Response
      If BLAST is successful, returns io.StringIO of response text.
      If initial BLAST query submission is unsuccessful, returns the Response object.
    '''
    assert (query is not None) or (sequence is not None)
    kwargs.update({
        'program': program,
        'database': database,
        'query': query,
        'sequence': sequence,
        'job_title': job_title
    })
    if semaphore is None:
        semaphore = multiprocessing.Semaphore()

    submit_args = ['program', 'database', 'query', 'sequence', 'url_base',
                   'composition_based_statistics', 'entrez_query', 'matrix_name', 'matrix',
                   'expect', 'filter', 'gapcosts', 'hitlist_size', 'nucl_penalty', 'nucl_reward',
                   'threshold', 'word_size', 'lock', 'server_delay']
    get_args = ['rid', 'url_base', 'no_wait', 'alignments', 'descriptions', 'format_object',
                'format_type', 'ncbi_gi', 'hitlist_size', 'lock', 'initial_delay',
                'server_delay', 'response_delay']

    semaphore.acquire()
    rid = blast_submit(**{arg: kwargs[arg] for arg in submit_args if arg in kwargs})
    if not isinstance(rid, str):
        semaphore.release()
        return rid

    if queue is not None:
        if job_title is None:
            job_title = sequence
        queue.put((job_title, rid))

    kwargs['rid'] = rid
    result = blast_get(**{arg: kwargs[arg] for arg in get_args if arg in kwargs})
    semaphore.release()
    return result

def blastAndParse(blast_fun=my_qblast, save=None, parse=None, parse_type=None, parse_name=None, log=sys.stdout, **kwargs):
    '''
    BLAST sequences and save results to file.

    Args
    - blast_fun: function. default=my_qblast
        BLAST function. Should accept kwargs and return a file object (subclass of io.IOBase, e.g., io.StringIO).
    - save: str. default=None
        Path to save raw BLAST results.
    - parse: str. default=None
        Directory in which to save parsed BLAST results (e.g., results from individual queries).
        Filenames are constructed from query IDs.
        If None, results are not parsed.
    - parse_type: str. default=None
        Format of results returned by blast_fun. Supported values: 'XML', 'Tabular'.
        If None, looks for `format_type` in kwargs.
    - parse_name: function. default=None
        Function that accepts an individual parsed result (Bio.SearchIO._model.query.QueryResult object,
        e.g., an element returned by Bio.SearchIO.parse()) and outputs a filename (str). If None, constructs
        a filename using the query result ID, parse_type, and applies gzip compression.
    - log: str. default=sys.stdout
        Path to save log file.
    - **kwargs: additional arguments to pass to blast_fun()

    Returns: return value of blast_fun, or None if blast_fun raises an exception.

    Notes
    - Any (sub)directories required for the paths provided in save and parse will be automatically created.
    - Gzip-compression is automatically applied if a `save` or `parse_name` path ends with '.gz'
    '''
    # Formats supported by BioPython's SearchIO module
    # Note that SearchIO supports reading BLAST Text format via 'blast-text', but writing is not supported
    parse_types = {'XML': 'blast-xml', 'Tabular': 'blast-tab'}
    extension = {'XML': '.xml.gz', 'Tabular': '.tsv.gz'}

    try:
        result = blast_fun(**kwargs)
    except Exception as err:
        print(err, file=log)
        return None

    if not issubclass(type(result), io.IOBase):
        print(('blast_fun() did not return a file object (subclass of io.IOBase). '
               'Results could not be saved.'), file=log)
        return result

    if save is not None:
        Path(save).parent.mkdir(parents=True, exist_ok=True)
        with utils_files.createFileObject(save, 'wt') as f:
            f.write(result.read())
        result.seek(0)

    if parse is not None:
        if parse_type is None:
            if 'format_type' in kwargs:
                parse_type = kwargs['format_type']
        else:
            print('No parse type provided.', file=log)
            return result

        if parse_type in parse_types:
            Path(parse).parent.mkdir(parents=True, exist_ok=True)
            for qresult in Bio.SearchIO.parse(result, format=parse_types[parse_type]):
                if parse_name is None:
                    filename = utils_files.createValidPath(qresult.id) + extension[parse_type]
                else:
                    filename = parse_name(qresult)
                with utils_files.createFileObject(Path(parse, filename), 'wt') as f:
                    Bio.SearchIO.write(qresult, f, parse_types[parse_type])
        else:
            print(f'Parsing of format_type {parse_type} is not supported.', file=log)
    return result

# def blastp_post(**kwargs):
#     '''
#     Submit query to NCBI BLAST server using a POST request.

#     A POST request is what the web form (https://blast.ncbi.nlm.nih.gov/Blast.cgi) submits.
#     An alternative is to use PUT-based requests, as implemented by Bio.Blast.NCBIWWW.qblast()
#     in the BioPython package.

#     References
#     - https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
#         Web form documentation

#     General
#     - QUERY: FASTA, bare sequence, or newline-delimited identifiers
#     - db: nucleotide (blastn, tblastn), protein (blastp, blastx)
#     - QUERY_FROM: Specify segment of the query sequence to BLAST
#     - QUERY_TO: Specify segment of the query sequence to BLAST
#     - QUERYFILE: Path to FASTA file
#     - JOB_TITLE: Title of BLAST search
#     - SUBJECTS: FASTA, bare sequence, or newline-delimited identifiers
#     - stype: nucleotide (blastn, tblastn), protein (blastp, blastx)
#     - SUBJECTS_FROM: Specify segment of the subject sequence to BLAST
#     - SUBJECTS_TO: Specify segment of the subject sequence to BLAST
#     - SUBJECTFILE: Path to FASTA file
#     - EQ_MENU: Restrict search to specified organisms. Ex: 'Fungi (taxid:4751)'
#     - ORG_EXCLUDE: Set to 'on' to exclude the organisms specified in EQ_MENU

#     blastp
#     - DATABASE: default=nr
#         nr (Non-redundant protein sequences), refseq_protein (RefSeq), landmark (Model Organisms),
#         swissprot (Swiss-Prot), pataa (Patented protein sequences), pdb (Protein Data Bank),
#         env_nr (Metagenomic proteins), or tsa_nr (Transcriptome Shotgun Assembly proteins)
#     - BLAST_PROGRAMS: default='blastp'
#         kmerBlastp (Quick BLASTP), blastp, psiBlast (PSI-BLAST), phiBlast (PHI-BLAST), or
#         deltaBlast (DELTA-BLAST)
#     - MAX_NUM_SEQ: default=[blastp, kmerBlastp] 100; [psiBlast, phiBlast, deltaBlast] 500
#         Maximum number of aligned sequences to display (web interface ranges from 100 to 20000)
#     - SHORT_QUERY_ADJUST: default='on'
#         Automatically adjust parameters for short input sequences
#     - EXPECT: E-value threshold. default=10
#     - WORD_SIZE: default=[blastp, kmerBlastp] 6; [psiBlast, phiBlast, deltaBlast] 3
#         Length of the seed that initiates an alignment.
#     - HSP_RANGE_MAX: default=0
#         Limit the number of matches to a query range.
#     - MATRIX_NAME: default=BLOSUM62
#         PAM30, PAM70, PAM250, BLOSUM80, BLOSUM62, BLOSUM45, BLOSUM50, or BLOSUM90
#     - GAPCOSTS: default='11 1'
#         Cost to create and extend a gap in an alignment.
#     - COMPOSITION_BASED_STATISTICS: default=2
#         Matrix adjustment method to compensate for amino acid composition of sequences.
#         0 (No adjustment), 1 (Composition-based statistics), 2 (Conditional compositional
#         score matrix adjustment), 3 (Universal compositional score matrix adjustment)

#     Unset
#     - PHI_PATTERN: [phiBlast only] PHI pattern to start the search
#     - PSSM: [psiBlast, phiBlast only] Path to PSSM file
#     - I_THRESH: [psiBlast, phiBlast, deltaBlast only] default=0.005
#         Statistical significance threshold to include a sequence in the model used by PSI-BLAST to create
#         the PSSM on the next iteration. 
#     - DI_THRESH: [deltaBlast only] default=0.05
#         Statistical significance threshold to include a domain in the model used by DELTA-BLAST to create the PSSM 
#     - PSI_PSEUDOCOUNT: [psiBlast, phiBlast, deltaBlast only] default=0
#         Pseduocount parameter. If zero is specified, then the parameter is automatically determined through
#         a minimum length description principle (PMID 19088134). A value of 30 is suggested in order to obtain
#         the approximate behavior before the minimum length principle was implemented.
#     - FILTER: include this argument once for each filter to set
#         L (low complexity regions), M (mask for lookup table only)
#     - LCASE_MASK: set to 'on' to enable mask
#         Mask any letters that were lower-case in the FASTA input.

#     Hidden fields for default formatting parameters
#     - SHOW_OVERVIEW: on
#     - SHOW_LINKOUT: on
#     - GET_SEQUENCE: on
#     - FORMAT_OBJECT: Alignment
#     - FORMAT_TYPE: HTML
#     - ALIGNMENT_VIEW: Pairwise
#     - MASK_CHAR: 2
#     - MASK_COLOR: 1
#     - DESCRIPTIONS: 100
#     - ALIGNMENTS: 100
#     - LINE_LENGTH: 60
#     - NEW_VIEW: on
#     - OLD_VIEW: false
#     - NCBI_GI: ''
#     - SHOW_CDS_FEATURE: ''
#     - NUM_OVERVIEW: 100
#     - FORMAT_EQ_TEXT: ''
#     - FORMAT_ORGANISM: ''
#     - EXPECT_LOW: ''
#     - EXPECT_HIGH: ''
#     - PERC_IDENT_LOW: ''
#     - PERC_IDENT_HIGH: ''
#     - QUERY_INDEX: 0
#     - FORMAT_NUM_ORG: 1
#     - CONFIG_DESCR: 2,3,4,5,6,7,8

#     Hidden fields
#     - CLIENT: web
#     - SERVICE: [blastp, kmerBlastp, psiBlast, phiBlast] plain; [deltaBlast] delta_blast
#     - CMD: request
#     - PAGE: Proteins
#     - PROGRAM: blastp
#     - MEGABLAST: ''
#     - RUN_PSIBLAST: [blastp, kmerBlastp] ''; [psiBlast, phiBlast, deltaBlast] on
#     - WWW_BLAST_TYPE: ''
#     - CDD_SEARCH: on
#     - ID_FOR_PSSM: ''
#     - DB_DISPLAY_NAME: same as DATABASE
#     - ORG_DBS: ''
#     - SAVED_PSSM: ''
#     - SELECTED_PROG_TYPE: same as DATABASE
#     - SAVED_SEARCH: ''
#     - BLAST_SPEC: ''
#     - MIXED_DATABASE: ''
#     - QUERY_BELIEVE_DEFLINE: ''
#     - DB_DIR_PREFIX: ''
#     - USER_DATABASE: ''
#     - USER_WORD_SIZE: ''
#     - USER_MATCH_SCORES: ''
#     - USER_FORMAT_DEFAULTS: ''
#     - NO_COMMON: ''
#     - NUM_DIFFS: 
#     - NUM_OPTS_DIFFS
#     - UNIQ_DEFAULTS_NAME
#     - PAGE_TYPE: BlastSearch
#     '''

#     default_kwargs = {
#         'DATABASE': 'nr',
#         'db': 'protein',
#         'BLAST_PROGRAMS': 'blastp',
#         'MAX_NUM_SEQ': 100,
#         'WORD_SIZE': 3,
#         'SHORT_QUERY_ADJUST': 'on',
#         'EXPECT': 10,
#         'HSP_RANGE_MAX': 0,
#         'MATRIX_NAME': 'BLOSUM62',
#         'GAPCOSTS': '11 1',
#         'COMPOSITION_BASED_STATISTICS': 2
#     }

#     if 'BLAST_PROGRAMS' in kwargs and kwargs['BLAST_PROGRAMS'] in ['psiBlast', 'phiBlast', 'deltaBlast']:
#         default_kwargs['MAX_NUM_SEQ'] = 500
#         default_kwargs['WORD_SIZE'] = 3

#     FILE_ARGS = ['QUERYFILE', 'SUBJECTFILE', 'PSSM']
#     for key in kwargs:
#         if key in FILE_ARGS:
#             kwargs[key] = (kwargs[key], open(kwargs[key], 'rb'), 'application/octet-stream')
#         else:
#             kwargs[key] = (None, kwargs[key])

#     r = requests.post(url=Bio.Blast.NCBIWWW.NCBI_BLAST_URL, files=data, allow_redirects=False)
#     return r

def search_NCBIGene(term, email=None, useDefaults=True, useSingleIndirectMatch=True, verbose=True, **kwargs):
    '''
    Search official human gene names and aliases in NCBI Gene database for a match to term, returning official IDs,
    names, and aliases.

    Dependencies: Biopython

    Args
    - term: iterable of str
        Gene name / alias
    - email: str. default=None
        Email registered with NCBI, used to set Bio.Entrez.email. If None, defaults to 'A.N.Other@example.com'.
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
    Bio.Entrez.email = email if email is not None else 'A.N.Other@example.com'

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
        if isinstance(values, str) or not isinstance(values, collections.abc.Iterable):
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

def multisearch_NCBIGene(terms, nproc=None, **kwargs):
    '''
    Search official human gene names and aliases in NCBI Gene database for terms.

    Args
    - terms: iterable of str
        Terms to lookup
    - email: str
        email registered with NCBI
    - nproc: int. default=None
        Number of processes to use. If None, uses as many processes as available CPUs.

    Returns: dict: str -> dict: int -> dict: str -> str or list of str
      {<term>: <NCBI Gene ID>: {'symbol': <gene symbol>, 'name': <gene name>, 'aliases': [aliases]}}
    '''

    if nproc is None:
        try:
            nproc = len(os.sched_getaffinity(0))
        except AttributeError:
            nproc = os.cpu_count()
    assert isinstance(nproc, int) and nproc > 0

    results = []
    with multiprocessing.Pool(processes=nproc) as pool:
        for i in range(len(terms)):
            results.append(pool.apply_async(search_NCBIGene, (terms[i],), kwargs))
        return {terms[i]: results[i].get() for i in range(len(terms))}

# endregion --- NCBI

# region ------ IDT

# API documentation: https://www.idtdna.com/restapi/swagger/ui/index#/
# - See also https://www.idtdna.com/pages/tools/apidoc
# - IDT apparently used to have a servers hosted at http://www.idtdna.com/SciTools/SciTools.aspx and
#   https://www.idtdna.com/AnalyzerService, the former of which was described in the publication
#   https://doi.org/10.1093/nar/gkn198

IDT_URL = 'https://www.idtdna.com'
IDT_API_URL = 'https://www.idtdna.com/restapi/v1/'

# IDT account found via http://bugmenot.com/view/idtdna.com
IDT_USERNAME = 'ev376821'
IDT_PASSWORD = 'evan1234'

IDT_DEFAULT_RENAME = {
    # OligoAnalyzer Analyze results
    "MeltTemp": "melt_temp",
    "MolecularWeight": "molecular_weight",
    "ExtCoefficient": "extinction_coefficient",
    "UgOD": "ug_od260",

    # OligoAnalyzer Hairpin results
    "deltaG": "min_hairpin_G",
}

def IDT_get_access_token(client_id, client_secret, username=IDT_USERNAME, password=IDT_PASSWORD):
    """
    Get access token for IDT SciTools Plus API. Adapted from https://www.idtdna.com/pages/tools/apidoc.

    Args
    - client_id: str
        IDT API Client ID
    - client_secret: str
        IDT API Client secret
    - username: str. default=IDT_USERNAME
        IDT account username. Apparently does not have to be the account for which the API key was generated.
    - password: str. default=IDT_PASSWORD
        IDT account password

    Returns
    - access_token: str
        The access token that can be used to authenticate requests to the IDT API.
    """
    from base64 import b64encode
    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = {
        "Content-Type": "application/x-www-form-urlencoded",
        "Authorization": "Basic " + authorization_string
    }
    request_data = {
        "grant_type": "password",
        "scope": "test",
        "username": username,
        "password": password
    }

    r = requests.post(
        "https://www.idtdna.com/Identityserver/connect/token",
        data=request_data,
        headers=request_headers
    )

    if r.status_code != 200:
        raise RuntimeError("Request failed with error code:" + str(r.status_code) + "\nBody:\n" + r.text)

    response = r.json()
    return response["access_token"]


def IDT_OligoAnalyzer_Analyze(access_token, seq, rename=False, settings=None):
    """
    Submit sequences to IDT OligoAnalyzer via IDT SciTools Plus API.
    See https://www.idtdna.com/calc/analyzer (requires login) and https://www.idtdna.com/calc/Analyzer/Home/Instructions.

    Args
    - access_token: str
        Access token for IDT SciTools Plus API
    - seq: str
        Nucleic acid sequence. Whitespace is ignored.
        - Modified bases allowed: see https://www.idtdna.com/site/catalog/Modifications/GetAllMods
        - Mixed bases allowed: see https://www.idtdna.com/calc/Analyzer/Home/Instructions and
          https://www.idtdna.com/calc/Analyzer/Home/definitions#MixedBaseDef
    - rename: bool or dict. default=False
        Rename results to a standardized nomenclature.
        If True, renames according to IDT_DEFAULT_RENAME.
        If a dict, the keys are renamed to the corresponding values.
    - settings: dict. default=None
        OligoAnalyzer settings.
        If None, defaults to nucleotide type "DNA", 50 mM Na+, 0 mM Mg2+, 0 mM dNTPs, 0.25 uM oligo.
        Acceptable range:
        - NucleotideType: "DNA" or "RNA"
        - NaConc: 5 - 1500 mM
        - MgConc: 0 - 600 mM
        - DNTPsConc: 0 - 1000 mM
        - OligoConc: 0.0001 - 100000 uM

    Returns: dict (order is not guaranteed)
    - BaseCount: int
    - Complement: str
        Complemet sequence with a space between every 3 nucleotides
    - Errors: list of dict
    - ExtCoefficient: float
        Units = L/mol/cm. Molar extinction coefficient at 260 nm wavelength.
    - GCContent: float
        Percent GC content
    - HasErrors: bool
    - HasModelErrors: bool
    - Length: int
    - MaxMeltTemp: float
        Units = Celcius. If the sequence contains mixed bases (e.g., N, R, Y, etc.), this value is the maximum melting
        temperature; otherwise, 0.
    - MeltTemp: float
        Units = Celcius. Melting temperature. If the sequence contains mixed bases, this value is the mean melting
        temperature.
    - MinMeltTemp: float
        Units = Celcius. If the sequence contains mixed bases, this value is the minimum melt temperature; otherwise, 0.
    - ModelErrors: None or list of str
    - Mods: None or list of dict
        List of modifications
    - MolecularWeight: float
        Units = g/mol
    - NmoleOD: float
        Units = nmol/mL
    - UgOD: float
        Units = ug/mL. Inverse of the molar attenuation coefficient at 260 nm wavelength and 1 cm path length
    - MgConc: float
        Units = mM. Mg2+ concentration
    - NaConc: float
        Units = mM. Na+ concentration
    - OligoConc: float
        Units = uM. Oligo concentration
    - dNTPsConc: float
        Units = mM. dNTP concentration
    - NucleotideType: str
        "DNA" or "RNA"
    - Sequence: str
        Query sequence with a space between every 3 nucleotides.
    """
    headers = {
        "Authorization": "Bearer " + access_token,
        "Accept": "application/json"
    }
    data = {
        "Sequence": seq,
        "NucleotideType": "DNA", # DNA or RNA
        "NaConc": 50, # mM
        "MgConc": 0, # mM
        "DNTPsConc": 0, # mM
        "OligoConc": 0.25 # uM
    }
    if settings is not None:
        assert "Sequence" not in settings
        data.update(settings)

    r = requests.post(
        IDT_API_URL + 'OligoAnalyzer/Analyze',
        json=data,
        headers=headers
    )
    if r.status_code != 200:
        raise RuntimeError("Request failed with error code:" + str(r.status_code) + "\nBody:\n" + r.text)

    response = r.json()
    if rename:
        name_map = IDT_DEFAULT_RENAME if rename is True else rename
        for key in name_map:
            if key in response:
                response[name_map[key]] = response.pop(key)
    return response

def IDT_OligoAnalyzer_Hairpin(access_token, seq, rename=False, settings=None):
    """
    Submit sequences to IDT OligoAnalyzer Hairpin via IDT SciTools Plus API.

    Args
    - access_token: str
        Access token for IDT SciTools Plus API
    - seq: str or list of str
        Nucleic acid sequence(s). Whitespace is ignored.
        - Modified bases allowed: see https://www.idtdna.com/site/catalog/Modifications/GetAllMods
        - Mixed bases allowed: see https://www.idtdna.com/calc/Analyzer/Home/Instructions and
          https://www.idtdna.com/calc/Analyzer/Home/definitions#MixedBaseDef
    - rename: bool or dict. default=False
        Rename results to a standardized nomenclature.
        If True, renames according to IDT_DEFAULT_RENAME.
        If a dict, the keys are renamed to the corresponding values.
    - settings: dict. default=None
        Acceptable range (defaults in parentheses):
        - NaConc: 5 - 1500 mM (50)
            The server doesn't seem to enforce this range.
        - MgConc: 0 - 600 mM (0)
        - FoldingTemp: 0 - 100 Celcius (25)
        - NucleotideType: "DNA" or "RNA" ("DNA")
            Any value other than "DNA" (case-insensitive) is apparently treated as RNA

    Returns: list of dict
    - Each entry of the list corresponds to an input sequence. It appears as if only results from the hairpin with the
      lowest Gibbs energy is returned.
      - id: str
      - sequence: str
          Query sequence (exactly as submitted)
      - thermo: float
          Units = Celcius. Melting temperature
      - deltaS: float
          Units = cal/mol/K. Entropy
      - deltaG: float
          Units = kcal/mol. Gibbs free energy
      - deltaH: float
          Units = kcal/mol. Enthalpy
      - success: bool
      - RNA: bool
          Whether the query sequence was submitted as RNA
      - errorMessages: list
    """
    headers = {
        "Authorization": "Bearer " + access_token,
        "Accept": "application/json"
    }
    data = {
        "NaConc": 50, # mM,
        "MgConc": 0, # mM,
        "FoldingTemp": 25, # Celcius,
        "NucleotideType": "DNA"
    }
    if settings is not None:
        assert "Sequence" not in settings and "Sequences" not in settings
        data.update(settings)
    if isinstance(seq, list):
        data["Sequences"] = seq
        url = IDT_API_URL + 'OligoAnalyzer/HairpinBatch'
    elif isinstance(seq, str):
        data["Sequence"] = seq
        url = IDT_API_URL + 'OligoAnalyzer/Hairpin'
    else:
        raise ValueError("seq must be a str or list of str")

    r = requests.post(url, json=data, headers=headers)
    if r.status_code != 200:
        raise RuntimeError("Request failed with error code:" + str(r.status_code) + "\nBody:\n" + r.text)

    response = r.json()
    if rename:
        for entry in response:
            name_map = IDT_DEFAULT_RENAME if rename is True else rename
            for key in name_map:
                if key in entry:
                    entry[name_map[key]] = entry.pop(key)
    return response

def IDT_OligoAnalyzer_SelfDimer(access_token, seq):
    """
    Submit sequences to IDT OligoAnalyzer SelfDimer via IDT SciTools Plus API.

    Args
    - access_token: str
        Access token for IDT SciTools Plus API
    - seq: str
        Nucleic acid sequence. Whitespace is ignored.
        - Modified bases allowed: see https://www.idtdna.com/site/catalog/Modifications/GetAllMods
        - Mixed bases allowed: see https://www.idtdna.com/calc/Analyzer/Home/Instructions and
          https://www.idtdna.com/calc/Analyzer/Home/definitions#MixedBaseDef

    Returns: list of dict
    - Each entry of the list corresponds to a potential self dimer structure.
      - StartPosition: int
          ??
      - TopLinePadding: int
      - BondLinePadding: int
      - BottomLinePadding: int
      - Bonds: list of int
          0 = no bond
          1 = additional complementary bases for that dimer structure; does not impact DeltaG value
          2 = part of the longest stretch of complementary bases; contributes to DeltaG value
      - DeltaG: float
          Units = kcal/mol. Gibbs free energy from the longest stretch of complementary bases
      - BasePairs: int
          Number of basepairs formed
      - Dimer: None
          ??
    """
    headers = {
        "Authorization": "Bearer " + access_token,
        "Accept": "application/json"
    }
    r = requests.post(
        IDT_API_URL + 'OligoAnalyzer/SelfDimer',
        params=dict(primary=seq),
        headers=headers
    )
    if r.status_code != 200:
        raise RuntimeError("Request failed with error code:" + str(r.status_code) + "\nBody:\n" + r.text)

    response = r.json()
    return response

def IDT_OligoAnalyzer_HeteroDimer(access_token, seq1, seq2):
    """
    Submit sequences to IDT OligoAnalyzer SelfDimer via IDT SciTools Plus API.

    Args
    - access_token: str
        Access token for IDT SciTools Plus API
    - seq1: str
        Nucleic acid sequence. Whitespace is ignored.
        - Modified bases allowed: see https://www.idtdna.com/site/catalog/Modifications/GetAllMods
        - Mixed bases allowed: see https://www.idtdna.com/calc/Analyzer/Home/Instructions and
          https://www.idtdna.com/calc/Analyzer/Home/definitions#MixedBaseDef
    - seq2: str
        Nucleic acid sequence. Whitespace is ignored.
        - Modified bases allowed: see https://www.idtdna.com/site/catalog/Modifications/GetAllMods
        - Mixed bases allowed: see https://www.idtdna.com/calc/Analyzer/Home/Instructions and
          https://www.idtdna.com/calc/Analyzer/Home/definitions#MixedBaseDef

    Returns: list of dict
    - Each entry of the list corresponds to a potential dimer structures.
      - StartPosition: int
          ??
      - TopLinePadding: int
      - BondLinePadding: int
      - BottomLinePadding: int
      - Bonds: list of int
          0 = no bond
          1 = additional complementary bases for that dimer structure; does not impact DeltaG value
          2 = part of the longest stretch of complementary bases; contributes to DeltaG value
      - DeltaG: float
          Units = kcal/mol. Gibbs free energy from the longest stretch of complementary bases
      - BasePairs: int
          Number of basepairs formed
      - Dimer: None
          ??
    """
    headers = {
        "Authorization": "Bearer " + access_token,
        "Accept": "application/json"
    }
    r = requests.post(
        IDT_API_URL + 'OligoAnalyzer/HeteroDimer',
        params=dict(primary=seq1, secondary=seq2),
        headers=headers
    )
    if r.status_code != 200:
        raise RuntimeError("Request failed with error code:" + str(r.status_code) + "\nBody:\n" + r.text)

    response = r.json()
    return response


def draw_dimer(seq1, seq2, entry, deltaG=True):
    """
    Visualize homodimer.

    Args
    - seq1: str
        Nucleic acid sequence (no modified bases)
    - seq2: str
        Nucleic acid sequence (no modified bases)
    - entry: dict
      - TopLinePadding: int
      - BondLinePadding: int
      - BottomLinePadding: int
      - Bonds: list of int
          0 = no bond
          1 = additional complementary bases for that dimer structure; does not impact DeltaG value
          2 = part of the longest stretch of complementary bases; contributes to DeltaG value
      - DeltaG: float
          Units = kcal/mol. Gibbs free energy from the longest stretch of complementary bases
      - BasePairs: int
          Number of basepairs formed
    - deltaG: bool. default=True
        Also print DeltaG value and number of base pairs, matching the web interface output.

    Returns: str
      Print this string to visualize homodimer structure
    """
    seq1 = "".join(seq1.split()) # remove all whitespace
    seq2 = "".join(seq2.split()) # remove all whitespace

    top = "5' " + ' ' * entry['TopLinePadding'] + seq1
    bottom = "3' " + ' ' * entry['BottomLinePadding'] + seq2[::-1]
    bonds = ' ' * (entry['BondLinePadding'] + 3)
    bond_map = {0: ' ', 1: ':', 2: '|'}
    for bond in entry['Bonds']:
        bonds += bond_map[bond]

    alignment = '\n'.join([top, bonds, bottom])
    if deltaG:
        analysis = f"Delta G: {entry['DeltaG']:.2f} kcal/mol, Base Pairs: {entry['BasePairs']}"
        return analysis + '\n' + alignment
    return alignment


def IDT_OligoAnalyzer_all(access_token, seq1, seq2=None, return_full=False, settings_oligo=None, settings_hairpin=None):
    """
    Submit sequences to IDT OligoAnalyzer using IDT SciTools Plus API.
    See https://www.idtdna.com/calc/analyzer (requires login).

    Args
    - access_token: str
        Access token for IDT SciTools Plus API
    - seq1: str
        Nucleic acid sequence (e.g., forward primer)
    - seq2: str. default=None
        Nucleic acid sequence (e.g., reverse primer)
    - return_full: bool. default=False
        Return full JSON response as dict.
    - settings_oligo: dict. default=None
        OligoAnalyzer Analyze settings.
    - settings_hairpin: dict. default=None
        OligoAnalyzer Hairpin settings.

    Returns: dict
      If `return_full`: dict representing decoded JSON results.
        - analyze: dict
        - hairpin: list of dict
        - homodimer: list of dict
        - heterodimer: list of dict
      Otherwise, a dict (str -> float):
        melt_temp: melt temperature (Celcius)
        molecular_weight: molecular weight (g/mol)
        extinction_coefficient: extinction coefficient (L/mol/cm)
        min_hairpin_G: smallest (i.e., most negative) ΔG (kcal/mol) of predicted hairpins
        min_homodimer_G': smallest (i.e., most negative) ΔG (kcal/mol) of homodimers (self-dimers)
        min_heterodimer_G': smallest (i.e., most negative) ΔG (kcal/mol) of seq1-seq2 heterodimer
          If seq2 is not provided, this value is None.
    """
    result_analyze = IDT_OligoAnalyzer_Analyze(access_token, seq1, settings=settings_oligo)
    result_hairpin = IDT_OligoAnalyzer_Hairpin(access_token, seq1, settings=settings_hairpin)[0]
    result_homodimer = IDT_OligoAnalyzer_SelfDimer(access_token, seq1)
    if seq2:
        result_heterodimer = IDT_OligoAnalyzer_HeteroDimer(access_token, seq1, seq2)
    else:
        result_heterodimer = None
    if return_full:
        result = dict(
            analyze=result_analyze,
            hairpin=result_hairpin,
            homodimer=result_homodimer,
            heterodimer=result_heterodimer
        )
    else:
        result = {
            'melt_temp': result_analyze['MeltTemp'],
            'molecular_weight': result_analyze['MolecularWeight'],
            'extinction_coefficient': result_analyze['ExtCoefficient'],
            'ug_od260': result_analyze['UgOD'],
            'min_hairpin_G': result_hairpin['deltaG'],
            'min_homodimer_G': min((x['DeltaG'] for x in result_homodimer)) if len(result_homodimer) > 0 else None,
            'min_heterodimer_G': min((x['DeltaG'] for x in result_heterodimer)) if result_heterodimer is not None and len(result_heterodimer) > 0 else None
        }
    return result

def IDT_OligoAnalyzer(seq1, seq2=None, return_full=False, username=IDT_USERNAME, password=IDT_PASSWORD):
    '''
    Submit sequences to IDT OligoAnalyzer via web interface.
    See https://www.idtdna.com/calc/analyzer (requires login).

    Args
    - seq1: str
        Nucleic acid sequence (e.g., forward primer)
    - seq2: str. default=None
        Nucleic acid sequence (e.g., reverse primer)
    - return_full: bool. default=False
        Return full JSON response as dict.

    Returns: dict
      If `return_full`: dict representing decoded JSON results
      Otherwise, a dict (str -> float):
        melt_temp: melt temperature (Celcius)
        molecular_weight: molecular weight (g/mol)
        extinction_coefficient: extinction coefficient (L/mol/cm)
        min_hairpin_G: smallest (i.e., most negative) ΔG (kcal/mol) of predicted hairpins
        min_homodimer_G': smallest (i.e., most negative) ΔG (kcal/mol) of homodimers (self-dimers)
        min_heterodimer_G': smallest (i.e., most negative) ΔG (kcal/mol) of seq1-seq2 heterodimer
          If seq2 is not provided, this value is None.
    '''
    with requests.Session() as s:
        # visit IDT
        r_idt = s.get(IDT_URL)
        if not r_idt.ok:
            r_idt.raise_for_status()

        url_login = 'https://www.idtdna.com/site/Account/Login/Gatekeeper'
        data_login = {'UserName': username, 'Password': password, 'RememberMe': 'true'}
        r_login = s.post(url_login, data=data_login)
        if not r_login.ok:
            r_login.raise_for_status()

        # defaults
        headers = {
            'origin': 'https://www.idtdna.com',
            'referer': 'https://www.idtdna.com/calc/analyzer'}

        # analyze
        url_analyze = 'https://www.idtdna.com/api/OligoAnalyzer/Analyze'
        data_analyze = {
            "Sequence": seq1,
            "NaConc": 50,
            "MgConc": 0,
            "DNTPsConc": 0,
            "OligoConc": 0.25,
            "NucleotideType": "DNA"}
        r_analyze = s.post(url_analyze, json=data_analyze, headers=headers)
        if not r_analyze.ok:
            r_analyze.raise_for_status()
        result_analyze = json.loads(r_analyze.text)

        # hairpin
        url_hairpin = 'https://www.idtdna.com/calc/analyzer/home/hairpin'
        data_hairpin = {
            'hairpinSettings': {
                "SequenceType": "Linear",
                "MaxFoldings": 20,
                "StartPos": 0,
                "StopPos": 0,
                "Temp": 25,
                "Suboptimality": 50},
            'settings': {
                "Sequence": seq1,
                "NaConc": 50,
                "MgConc": 0,
                "DNTPsConc": 0,
                "OligoConc": 0.25,
                "NucleotideType": "DNA"}}
        r_hairpin = s.post(url_hairpin, json=data_hairpin, headers=headers)
        if not r_hairpin.ok:
            r_hairpin.raise_for_status()
        result_hairpin = json.loads(r_hairpin.text)

        # self-dimer
        url_homodimer = 'https://www.idtdna.com/calc/analyzer/home/selfdimer'
        data_homodimer = {
            'settings': {
                "Sequence": seq1,
                "NaConc": 50,
                "MgConc": 0,
                "DNTPsConc": 0,
                "OligoConc": 0.25,
                "NucleotideType": "DNA"}}
        r_homodimer = s.post(url_homodimer, json=data_homodimer, headers=headers)
        if not r_homodimer.ok:
            r_homodimer.raise_for_status()
        result_homodimer = json.loads(r_homodimer.text)

        # heterodimer
        result_heterodimer = None
        if seq2 is not None:
            url_heterodimer = 'https://www.idtdna.com/calc/analyzer/home/heterodimer'
            data_heterodimer = {"primary": seq1, "secondary": seq2}
            r_heterodimer = s.post(url_heterodimer, json=data_heterodimer, headers=headers)
            if not r_heterodimer.ok:
                r_heterodimer.raise_for_status()
            result_heterodimer = json.loads(r_heterodimer.text)

        if return_full:
            return {
                'analyze': result_analyze,
                'hairpin': result_hairpin,
                'homodimer': result_homodimer,
                'heterodimer': result_heterodimer}

        return {
            'melt_temp': result_analyze['MeltTemp'],
            'molecular_weight': result_analyze['MolecularWeight'],
            'extinction_coefficient': result_analyze['ExtCoefficient'],
            'min_hairpin_G': min((x['GNode'] for x in result_hairpin['OutputObj']['Structures'])) if len(result_hairpin['OutputObj']['Structures']) > 0 else None,
            'min_homodimer_G': min((x['DeltaG'] for x in result_homodimer['Results'])) if len(result_homodimer['Results']) > 0 else None,
            'min_heterodimer_G': min((x['DeltaG'] for x in result_heterodimer['Results'])) if result_heterodimer is not None and len(result_heterodimer['Results']) > 0 else None}

def IDT_complexity_screen(seq, product_type='gblock'):
    '''
    Assess complexity of DNA sequence.

    Args
    - seq: str
        DNA sequence
    - product_type: str. default='gblock'
        gblock, gene, or megamer

    Returns: dict
      Decoded JSON response. Keys:
        ChosenIndices
        Complexities
        ComplexityScreenerResults
        FullSequence
        OptimizedSubseq
        RestrictionSites
        ComplexityColor
        ComplexityDetail
        ComplexityLimitReached
        ComplexityScore
        ComplexitySummary
    '''
    with requests.Session() as s:
        # visit IDT
        r_idt = s.get(IDT_URL)
        if not r_idt.ok:
            r_idt.raise_for_status()

        url = 'https://www.idtdna.com/CodonOpt/Home/ValidateItem'
        params = {'fullSeq': seq, 'subseq': seq, 'productType': product_type}
        r = s.post(url, json=params)
        if not r.ok:
            r.raise_for_status()
    return json.loads(r.text)

def IDT_codon_opt(seqs, sequence_type, taxon_id, product_type='gblock', return_full=False):
    '''
    Submit sequences to IDT Codon Optimization Tool.
    See https://www.idtdna.com/CodonOpt (requires login).

    Args
    - seqs: list of 2-tuple of str
        Sequences to codon optimize given as [(name1, seq1), (name2, seq2), ...]
    - sequence_type: str
        aminoAcid or dna
    - taxon_id: int or str
        Taxon ID of target organism. See taxon_id argument of downloadBio.get_codon_usage().
    - product_type: str. default='gblock'
        gblock, gene, or megamer
    - return_full: bool. default=False
        Return full JSON response as dict.

    Returns: pandas.DataFrame or dict
      If `return_full`: dict representing decoded JSON results
      Otherwise, a DataFrame:
        name: name of input sequence
        input: input sequence
        output: optimized sequence
        complexity: complexity score
          < 7: Accepted - Low Complexity
          < 10: Accepted - High Complexity
    '''
    assert sequence_type in ('aminoAcid', 'dna')

    with requests.Session() as s:
        # visit IDT
        r_idt = s.get(IDT_URL)
        if not r_idt.ok:
            r_idt.raise_for_status()

        # login to IDT - account found via http://bugmenot.com/view/idtdna.com
        url_login = 'https://www.idtdna.com/site/Account/Login/Gatekeeper'
        data_login = {'UserName': 'ev376821', 'Password': 'evan1234', 'RememberMe': 'true'}
        r_login = s.post(url_login, data=data_login)
        if not r_login.ok:
            r_login.raise_for_status()

        # build request
        url_copt = 'https://www.idtdna.com/CodonOpt/Home/Optimize'
        headers_copt = {
            'origin': 'https://www.idtdna.com',
            'referer': 'https://www.idtdna.com/CodonOpt',
        }
        data_copt = {
            'tableString': get_codon_usage(taxon_id, output='str').replace('\n', '\r\n'),
            'optimizationItems': [{
                'ItemModalError': {
                    'ShowModalMessage': 'false',
                    'ModalMessageType': 'modal-message-danger',
                    'ModalTitle': '',
                    'ModalMessage': ''
                },
                'ReadingFrameStart': -1,
                'ReadingFrameEnd': -1,
                'OptimizedSequence': '',
                'Id': idx,
                'Name': name,
                'OriginalSequence': seq.upper(),
                'SequenceType': sequence_type,
                'UpstreamSubseq': '',
                'DownstreamSubseq': '',
                'ValidatingItem': 'false',
                'Complexities': [],
                'RestrictionSiteWarning': '',
                'ManualOptimizing': 'false',
                'IsChecked': 'false',
                'AminoAcidResults': [],
                'OptSubseq': '',
                'Text': '',
                'NoComplexitiesDetected': 'false',
                'ShowComplexities': 'false',
                'ComplexitySummary': '',
                'ComplexityDetail': '',
                'ComplexityScore': 0,
                'ComplexityLimitReached': 'false',
                'ComplexityColor': ''
            } for idx, (name, seq) in enumerate(seqs)],
            'sequenceType': sequence_type,
            'productType': product_type
        }
        r_copt = s.post(url_copt, json=data_copt, headers=headers_copt)
        if not r_copt.ok:
            r_copt.raise_for_status()
        results = json.loads(r_copt.text)
        if return_full:
            return results
        return pd.DataFrame(
            [(
                result['Name'],
                result['OriginalSequence'],
                result['OptResult']['FullSequence'],
                result['OptResult']['ComplexityScore'])
                for result in results],
            columns = ['name', 'input', 'output', 'complexity'])

# endregion --- IDT

# region --- STRING

def get_interactors_string(
    identifiers: list[str],
    species: int = 10090,
    limit: int | None = None,
    required_score: int | float | None = None,
    network_type: str = 'functional',
    output_format: str = 'tsv',
    caller_identity: str = 'btyeh',
) -> str | pd.DataFrame:
    '''
    See https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set.

    Args
    - identifiers: protein identifiers
    - species: NCBI taxon ID
    - limit: number of interaction partners retrieved per protein
        If not specified, the STRING server appears to default to 10.
    - required_score: threshold of significance to include a interaction, a number between 0 and 1000
        Appears to correspond to 1000 * interaction score
    - network_type: 'functional' or 'physical'
    - output_format: 'tsv', 'tsv-no-header', 'json', 'xml', 'psi-mi', or 'psi-mi-tab'
    - caller_identity: user identity

    Returns
    - If output_format == 'tsv': returns a DataFrame
    - Otherwise, returns the server result as plain text.
    '''
    assert output_format in ['tsv', 'tsv-no-header', 'json', 'xml', 'psi-mi', 'psi-mi-tab']
    assert network_type in ['functional', 'physical']
    URL_STRING_PARTNERS = f'https://string-db.org/api/{output_format}/interaction_partners'
    params = dict(
        identifiers='\r'.join(identifiers),
        species=species,
        limit=limit,
        required_score=required_score,
        network_type=network_type,
        caller_identity=caller_identity
    )
    params = {k: v for k, v in params.items() if v is not None}
    r = requests.get(
        url=URL_STRING_PARTNERS,
        params=params
    )
    result = r.content.decode()
    if output_format == 'tsv':
        return pd.read_csv(io.StringIO(result), sep='\t')
    else:
        return result

# endregion --- STRING

# region ------ Miscellaneous

def get_codon_usage(taxon_id, codon_code=1, output='table', recalculate=False, dna=False):
    '''
    Retrieve codon usage from the Codon Usage Database hosted by the Kazusa DNA Research Institute.
    See http://www.kazusa.or.jp/codon/.

    Args
    - taxon_id: int or str
        Taxon ID of species. Append '.mitochondrion' to taxon ID for mitochondrion codon usage.
    - codon_code: int. default=1
        Codon translation code
          0: Do not translate
          1: Standard
          2: Vertebrate Mitochondrial
          3: Yeast Mitochondrial
          4: Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma
          5: Invertebrate Mitochondrial
          6: Ciliate Macronuclear and Dasycladacean
          9: Echinoderm Mitochondrial
          10: Alternative Ciliate Macronuclear
          11: Eubacterial
          12: Alternative Yeast
          13: Ascidian Mitochondrial
          14: Flatworm Mitochondrial
          15: Blepharisma Nuclear Code
    - output: str. default='table'
        table: pandas.DataFrame
          triplet: str
            Codon triplet
          amino acid: str
            [Only if codon_code > 0] Amino acid encoded by the codon triplet
          fraction: float
            [Only if codon_code > 0] For a given amino acid, the fraction coded with a particular codon triplet
          frequency: float
            Count per thousand over all CDSs of the organism.
            Sum of all frequencies in the table should equal 1000.
          number: int
            Codon count over all CDSs of the organism
        dict: dict (str -> dict (str -> float))
          Codon table mapping amino acid to dictionary of codon frequencies, such as returned by
          python_codon_tables.get_codons_table().
          Structure: {amino_acid: {codon1: freq1, codon2: freq2, ...}, ...}
        str: str
          Raw codon usage string
        html: str
          Webpage HTML
    - recalcuate: bool. default=False
        Recalculate frequency and fraction (if codon_code > 0) based on number
    - dna: bool. default=False
        Convert RNA codons to DNA (replace 'U' with 'T' in all codons)

    Returns: pandas.DataFrame or str
      See `output` argument.
    '''
    assert codon_code in (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15)
    assert output in ('dict', 'html', 'str', 'table')
    if dna and output in ('str', 'html'):
        print(f'Conversion of codons to DNA not supported for \'{output}\' output format.', file=sys.stdout)

    url = 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi'
    params = {'species': taxon_id}
    if codon_code != 0:
        params.update({'aa': codon_code, 'style': 'N'})
    r = requests.get(url, params=params)
    if not r.ok:
        r.raise_for_status()

    if output == 'html':
        return r.text

    # re.DOTALL flag necessary for the dot to match newlines
    codon_usage_str = re.search(r'<pre>(.*)</pre>', r.text, re.DOTALL | re.IGNORECASE).groups()[0].strip()
    if output == 'str':
        return codon_usage_str

    codon_usage_str = re.sub(pattern=r'\)\s+', repl=')\n', string=codon_usage_str)
    codon_usage_str = re.sub(pattern=r'\(|\)', repl=' ', string=codon_usage_str)
    codon_usage_str = re.sub(pattern=r'\n\s+', repl='\n', string=codon_usage_str)
    if codon_code == 0:
        columns = ('triplet', 'frequency', 'number')
    else:
        columns = ('triplet', 'amino acid', 'fraction', 'frequency', 'number')
    df = pd.read_csv(io.StringIO(codon_usage_str), sep='\s+', header=None, names=columns)
    if recalculate:
        df['frequency'] = df['number'] / (df['number'].sum() / 1000)
        if 'fraction' in df.columns:
            df['fraction'] = df['number'] / df.groupby('amino acid')['number'].transform('sum')
    if dna:
        df['triplet'] = df['triplet'].str.replace('U', 'T')
    if output == 'table':
        return df
    else: # output == 'dict'
        return convert_codon_usage_format(df, 'dict')

def convert_codon_usage_format(d, to):
    '''
    Convert codon table from pandas.DataFrame to dict.

    Args: see 'output' argument of get_codon_usage()
    - d: pandas.DataFrame or dict
    - to: str
        'table' or 'dict'

    Returns: see 'output' argument of get_codon_usage()
    '''
    if to == 'table':
        assert isinstance(d, dict)
        return pd.DataFrame(d).reset_index().rename(columns={'index': 'triplet'}) \
                              .melt(id_vars='triplet', var_name='amino acid', value_name='fraction') \
                              .dropna().sort_values(['amino acid', 'triplet']).reset_index(drop=True)
    elif to == 'dict':
        assert isinstance(d, pd.DataFrame)
        codon_usages = d.groupby('amino acid') \
          .apply(lambda group: group.pivot_table(index='amino acid', columns='triplet', values='fraction') \
                                    .to_dict('index')) \
          .tolist()
        return {aa: freqs for d in codon_usages for aa, freqs in d.items()}
    else:
        raise ValueError('"to" must be either "table" or "dict".')

def get_melt_curve(seq, resolution=0.5, monovalent_cations=20, mg=3, dmso=0):
    '''
    Predict fluorescent DNA melting curves of PCR products. See https://dna-utah.org/umelt/quartz/.

    Args
    - seq: str
        DNA sequence (i.e., PCR product; melt curve is calculated for homodimer)
    - resolution: float. default=0.5
        Resolution of melt curve prediction (in degrees Celcius)
    - monovalent_cations: float. default=20
        Concentration of monovalent cations in millimolar.
    - mg: float. default=3
        Concentration of Mg2+ in millimolar.
    - dmso: float. default=0
        Concentration of DMSO in %.

    Returns: pandas.DataFrame
    - temperature: T
    - helicity: %
    - derivative_neg: -dF/dT
    '''
    url_main = 'https://www.dna-utah.org/umelt/quartz/'
    with requests.Session() as s:
        r = s.get(url_main + 'um.php')
        token = re.search(r'&token=" \+ "(\S+)";', r.text).groups()[0]
        url = url_main + f'request.php?seq={seq}&density={resolution}&thermo=1&rs=0&cation={monovalent_cations}&mg={mg}&dmso={dmso}&token={token}'
        r2 = s.get(url)

    temperature = np.array(list(map(float, re.search('<temperature>([^<]+)</temperature>', r2.text).groups()[0].split(' '))))
    helicity = np.array(list(map(float, re.search('<helicity>([^<]+)</helicity>', r2.text).groups()[0].split(' '))))
    derivative = np.append((helicity[1:] - helicity[0:-1]) / (temperature[1:] - temperature[0:-1]), np.nan)
    return pd.DataFrame(dict(temperature=temperature, helicity=helicity, derivative_neg=-derivative))

# endregion --- Miscellaneous
