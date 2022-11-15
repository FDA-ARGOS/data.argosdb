#!/usr/bin/env python3
"""Assembly QC NCBI

This script will parse an XML for a single feature and return a list via terminal
or file, if an output is supplied.
"""
__version__ = "1.0.0"
__status__ = "PRODUCTION"

import csv
import os
import sys
import xml.etree.ElementTree as ET
from argparse import SUPPRESS, ArgumentParser
from ftplib import FTP

import requests
from Bio import Entrez


def usr_args():
    """User Arguments

    User supplied arguments from command line for function

    Returns
    -------
        ArgumentParser objects to be digested by subsequent functions.
    """

    parser = ArgumentParser(
        add_help=False,
        prog='assembly_qc_ncbi',
        description="Assembly QC NCBI. "
            "Used to gather metrics from NCBI assemblies.")

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input',
        required=True,
        help="Either 'api' for an api query or a file and path for parsing.")

    optional.add_argument('-o', '--output',
        help="Output file. If no output is provided, the resulting table will \
            be output to the terminal.")
    optional.add_argument('-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)
    optional.add_argument('-h', '--help',
        action='help',
        default=SUPPRESS,
        help='show this help message and exit')

    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def argosdb_api(assembly_stats, samples):
    """ARGOSDB API Query

    Will query the ARGOS DB API for data tables based on the BCO ID. 

    Parameters
    ----------
    assembly_stats: dict
        Final results dictionary. Mostly taken from the NCB genome stats reports.
    samples: dict
        Dictionary of values parsed from the genome XML file.

    """

    api_url = "https://beta-api.argosdb.org/records/search"
    results = []
    assemblies = []
    os.system('mkdir home/assembly')
    data = [{
        "bcoid": "ARGOS_000012",
        "offset": 1,
        "limit": 10000
    },
    {
        "bcoid": "ARGOS_000022",
        "offset": 1,
        "limit": 10000
    }]
    for item in data:
        response = requests.post(api_url, json=item)
        results.append(response.json())
        for record in response.json()['recordlist']:
            if record['genome_assembly_id'] not in assemblies:
                assemblies.append(record['genome_assembly_id'])
        ids = ','.join(assemblies)
        print(ids)
        os.system(f"efetch -db assembly -id {ids} -format docsum \
            > home/assembly/genome_assembly.xml")
    parse_xml('home/assembly/genome_assembly.xml', assembly_stats, samples)

def parse_xml(xml_file, assembly_stats, samples):
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    assembly_stats: dict
        A python dictionary of the assembly statitistics from NCBI with the
        genome_assembly_id as the index value.

    Returns
    -------
    assembly_stats: dict
        A python dictionary of the assembly statitistics from NCBI with the
        genome_assembly_id as the index value.
    """

    organism_name = genome_assembly_id \
        = lineage = taxonomy_id = bco_id = schema_version = analysis_platform \
        = analysis_platform_object_id = sra_run_id = ngs_read_file_source \
        = num_chromosomes = num_genes = assembly_type \
        = assembly_level = '-'

    analysis_platform = 'NCBI'
    bco_id = 'ARGOS_000038'
    schema_version = 'v1.0'
    root = ET.parse(xml_file).getroot()
    for item in root:
        if item.tag == "DbBuild":
            continue
        for summary in item.findall('./'):
            if summary.tag == 'FtpPath_Stats_rpt':
                print(summary.text)
                stats_report(summary.text, assembly_stats)

    for item in root:
        if item.tag == "DbBuild":
            continue
        for summary in item.findall('./'):
            if summary.tag == 'Synonym':
                for synonym in summary.findall('./'):
                    if synonym.tag == 'Genbank':
                        genome_assembly_id = synonym.text
            if summary.tag == 'Id':
                analysis_platform_object_id = summary.text
            if summary.tag == 'SpeciesName':
                organism_name = summary.text
            if summary.tag == 'Taxid':
                taxonomy_id = summary.text
                lineage = get_lineage(taxonomy_id)
            if summary.tag == 'AssemblyType':
                assembly_type = summary.text
            if summary.tag == 'FtpPath_GenBank':
                ngs_read_file_source = summary.text + '_genomic.gbff.gz'
            if summary.tag == 'Meta':
                for meta in summary.findall('./'):
                    if meta.tag == 'Stats':
                        for stat in meta.findall('./'):
                            if stat.attrib['category'] == 'chromosome_count':
                                num_chromosomes = stat.text
                            # if stat.attrib['category'] == 'total_length':
                            #     length = stat.text
                    if meta.tag == 'assembly-status':
                        assembly_level = meta.text

        samples[genome_assembly_id] = [organism_name, lineage, taxonomy_id,
            bco_id, schema_version, analysis_platform,
            analysis_platform_object_id, sra_run_id, ngs_read_file_source,
            num_chromosomes, num_genes, assembly_type,
            assembly_level]

    return assembly_stats, samples
def merge_results(assembly_stats, samples):
    """Merge Results

    Parameters
    ----------
    assembly_stats: dict
        Final results dictionary. Mostly taken from the NCB genome stats reports.
    samples: dict
        Dictionary of values parsed from the genome XML file.
    """

    for sequence in assembly_stats.items():
        assembly_stats[sequence[0]][0] = samples[assembly_stats[sequence[0]][2]][0]
        assembly_stats[sequence[0]][3] = samples[assembly_stats[sequence[0]][2]][1]
        assembly_stats[sequence[0]][4] = samples[assembly_stats[sequence[0]][2]][2]
        assembly_stats[sequence[0]][5] = samples[assembly_stats[sequence[0]][2]][3]
        assembly_stats[sequence[0]][6] = samples[assembly_stats[sequence[0]][2]][4]
        assembly_stats[sequence[0]][7] = samples[assembly_stats[sequence[0]][2]][5]
        assembly_stats[sequence[0]][8] = samples[assembly_stats[sequence[0]][2]][6]
        assembly_stats[sequence[0]][9] = samples[assembly_stats[sequence[0]][2]][7]
        assembly_stats[sequence[0]][10] = samples[assembly_stats[sequence[0]][2]][8]
        assembly_stats[sequence[0]][12] = samples[assembly_stats[sequence[0]][2]][9]
        assembly_stats[sequence[0]][13] = samples[assembly_stats[sequence[0]][2]][10]
        assembly_stats[sequence[0]][25] = samples[assembly_stats[sequence[0]][2]][11]
        assembly_stats[sequence[0]][26] = samples[assembly_stats[sequence[0]][2]][12]

    return assembly_stats, samples

def stats_report(stats_ftp, assembly_stats):
    """Assembly Stats Report

    Downloads and then parses the NCBI Assembly Stats Report

    Parameters
    ----------
    stats_ftp: str
        ftp url for NCBI assembly stats files
    assembly_stats: dict
        A dictionary for stroing results and creating final table
    """

    infraspecific_name = '-'
    refseq_assembly_id = '-'
    genbank_assembly_id = '-'
    count = 0
    # if os.path.exists('home/assembly/stats') is True:
    #     stats_file = 'GCF_000865725.1_ViralMultiSegProj15521_assembly_stats.txt'
    #     report_file = 'GCF_000865725.1_ViralMultiSegProj15521_assembly_report.txt'
    # else:
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    stats_ftp = stats_ftp.split('://')[1]
    stats_file = stats_ftp.split('/')[-1]
    report_file = stats_file.replace('assembly_stats','assembly_report')
    genome_dir = stats_ftp.split(stats_file)[0]
    genome_dir = genome_dir.replace('ftp.ncbi.nlm.nih.gov/', '')
    ftp.cwd(genome_dir)
    os.system('mkdir home/assembly/stats')

    with open('home/assembly/stats/' + stats_file, 'wb') as file:
        ftp.retrbinary(f'RETR {stats_file}', file.write )

    with open('home/assembly/stats/' + report_file, 'wb') as file:
        ftp.retrbinary(f'RETR {report_file}', file.write )

    with open('home/assembly/stats/' + report_file, encoding='utf-8') as  reportfile:
        report_tsv = csv.reader(reportfile, delimiter='\t')
        for row in report_tsv:
            if row[0].startswith('# GenBank assembly accession: '):
                genbank_assembly_id = row[0].replace('# GenBank assembly accession: ', '')

            if row[0].startswith('# RefSeq assembly accession:'):
                refseq_assembly_id= row[0].replace('# RefSeq assembly accession: ', '')

            if row[0].startswith('# Infraspecific name:'):
                infraspecific_name = row[0].replace('# Infraspecific name:  ', '')
            if row[0].startswith('#'):
                continue
            if row[6] == 'na':
                count += 1
                assembly_stats[row[0]] = ['-', infraspecific_name, genbank_assembly_id, '-', '-',
                    '-', '-', '-', '-', '-', '-', row[0], '-', '-', '-', '-', '-',
                    '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
            else:
                count += 1
                assembly_stats[row[6]] = ['-', infraspecific_name, genbank_assembly_id, '-', '-',
                    '-', '-', '-', '-', '-', '-', row[0], '-', '-', '-', '-', '-',
                    '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', refseq_assembly_id]

    with open('home/assembly/stats/' + stats_file, 'r', encoding='utf8') as statfile:
        genome_tsv = csv.reader(statfile, delimiter="\t")
        for row in genome_tsv:
            if row[0].startswith('# GenBank assembly accession: '):
                assembly_id = row[0].replace('# GenBank assembly accession: ', '')
            if row[0].startswith('#'):
                continue
            for genome in assembly_stats:
                segment = assembly_stats[genome][11]
                if assembly_id == assembly_stats[genome][2]:
                    assembly_stats[genome][14] = count
                    if row[0] == 'all' and row[3] == 'all':
                        if row[4] == 'scaffold-N75':
                            assembly_stats[genome][18] = row[5]
                        if row[4] == 'scaffold-N90':
                            assembly_stats[genome][19] = row[5]
                        if row[4] == 'scaffold-L50':
                            assembly_stats[genome][20] = row[5]
                    if row[1] == segment and row[3] == 'all':
                        if row[4] == 'gc-perc':
                            assembly_stats[genome][15] = row[5]
                        if row[4] == 'total-length':
                            assembly_stats[genome][16] = row[5]
                        if row[4] == 'scaffold-N50':
                            assembly_stats[genome][17] = row[5]
                        if row[4] == 'scaffold-N75':
                            assembly_stats[genome][18] = row[5]
                        if row[4] == 'scaffold-N90':
                            assembly_stats[genome][19] = row[5]
                        if row[4] == 'scaffold-L50':
                            assembly_stats[genome][20] = row[5]
                    if row[1] == 'na' and row[3] == 'all':
                        if row[4] == 'gc-perc':
                            assembly_stats[genome][15] = row[5]
                        if row[4] == 'total-length':
                            assembly_stats[genome][16] = row[5]
                        if row[4] == 'scaffold-N50':
                            assembly_stats[genome][17] = row[5]
                        if row[4] == 'scaffold-N75':
                            assembly_stats[genome][18] = row[5]
                        if row[4] == 'scaffold-N90':
                            assembly_stats[genome][19] = row[5]
                        if row[4] == 'scaffold-L50':
                            assembly_stats[genome][20] = row[5]
    print('count: ', count)
    return assembly_stats

def get_lineage(taxonomy_id):
    """Get lineage

    Retrieves the NCBI Taxonomic lineage given a tax id

    Parameters
    ----------
    taxonomy_id: str
        NCBI Taxonomy identifier
    """
    Entrez.email = "hadley_king@gwu.edu"
    Entrez.api_key = os.getenv('NCBI_API_KEY')
    handle = Entrez.efetch(db="taxonomy", id=taxonomy_id)
    record = Entrez.read(handle)
    return record[0]['Lineage']

def sample_output( samples, header, output):
    """Sample Output

    If an output file is supplied in the user arguments the samples dictionary
    will be output in TSV to the supplied file. If no output is provided, the
    samples dictionary will print to the terminal in a TSV format.

    Parameters
    ----------
    header: lst of str
        List of column headers for the output
    samples: dict
        dictionary of samples
    outputs: str, optional
        file path/name to output data to
    """

    if output:
        sample_file = os.path.abspath(output)
        with open(sample_file, 'w',  encoding='utf8') as file:
            writer = csv.writer(file, header, delimiter='\t')
            writer.writerow(header)
            for key in samples:
                row = [key]
                for item in range(len(samples[key])):
                    row.append(samples[key][item])
                writer.writerow(row)
        print('Samples written to ', sample_file)

    else:
        print('\t'.join(item for item in header))
        for run in samples:
            print(run, '\t', '\t'.join(str(item) for item in samples[run]))

def main():
    """Main"""

    assembly_stats = {}
    samples = {}
    header = ['ref_genome_acc', 'organism_name', 'infraspecific_name', \
        'genome_assembly_id', 'lineage', 'taxonomy_id', 'bco_id', \
        'schema_version', 'analysis_platform', 'analysis_platform_object_id', \
        'sra_run_id', 'ngs_read_file_source', 'genomic_section', \
        'num_chromosomes', 'num_genes', 'num_segments', 'assembly_gc_content', \
        'length', 'n50', 'n75', 'n90', 'l50', 'l75', \
        'query_coverage_against_reference', \
        'percent_identity_against_reference', 'percent_reads_unaligned', \
        'assembly_type', 'assembly_level', 'refseq_assembly_id']

    args = usr_args()
    if args.input == 'api':
        argosdb_api(assembly_stats, samples)
    else:
        parse_xml(args.input, assembly_stats, samples)
    merge_results(assembly_stats, samples)
    sample_output(assembly_stats, header, args.output)

if __name__ == "__main__":
    main()
