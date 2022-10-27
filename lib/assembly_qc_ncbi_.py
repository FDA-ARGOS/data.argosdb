#!/usr/bin/env python3
"""Assembly QC NCBI

This script will parse an XML for a single feature and return a list via terminal
or file, if an output is supplied.
"""
__version__ = "1.0.0"
__status__ = "PRODUCTION"

import os
import sys
import csv
from ftplib import FTP
import xml.etree.ElementTree as ET
from argparse import ArgumentParser, SUPPRESS
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

def argosdb_api(samples):
    """ARGOSDB API Query

    """

    api_url = "https://api.argosdb.org/records/search"
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
        os.system(f"efetch -db assembly -id {ids} -format docsum \
            > home/assembly/genome_assembly.xml")
    parse_xml('home/assembly/genome_assembly.xml', samples)

def parse_xml(xml_file, samples):
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    output: str, optional
        file path/name for data to be output to

    Returns
    -------
    samples: dict
        A python dictionary of the assembly statitistics from NCBI with the
        genome_assembly_id as the index value.
    """
    genome_stats = {}
    genome_assembly_id = lineage = taxonomy_id = bco_id = schema_version \
        = analysis_platform = analysis_platform_object_id = bioproject \
        = biosample = sra_run_id = ngs_read_file_name = ngs_read_file_source \
        = organism_name = ref_genome_acc = genomic_section \
        = num_chromosomes = num_genes = num_segments  = assembly_gc_content \
        = length = n50  = n75 = n90 = l50  = l75 = number_of_n \
        = percent_assembly_greater_20x = query_coverage_against_reference \
        = percent_identity_against_reference = reads_unaligned = id_method \
        = assembly_type = assembly_level = assembly_score = '-'

    analysis_platform = 'NCBI'
    bco_id = 'ARGOS_000038'
    schema_version = 'v0.9'
    genomic_section = 'all'

    root = ET.parse(xml_file).getroot()
    for item in root:
        if item.tag == "DbBuild":
            continue
        for summary in item.findall('./'):
            if summary.tag == 'AssemblyAccession':
                genome_assembly_id = summary.text
            if summary.tag == 'Id':
                analysis_platform_object_id = summary.text
            if summary.tag == 'SpeciesName':
                organism_name = summary.text
            if summary.tag == 'Taxid':
                taxonomy_id = summary.text
                lineage = get_lineage(taxonomy_id)
            if summary.tag == 'AssemblyType':
                assembly_type = summary.text
            if summary.tag == 'GB_BioProjects':
                if summary[0][0].tag == 'BioprojectAccn':
                    bioproject = summary[0][0].tag
            if summary.tag == 'RS_BioProjects':
                if summary[0][0].tag == 'BioprojectAccn':
                    bioproject = summary[0][0].text
            if summary.tag == 'BioSampleAccn':
                biosample = summary.text
            if summary.tag == 'FtpPath_GenBank':
                ngs_read_file_source = summary.text + '_genomic.gbff.gz'
                ngs_read_file_name = summary.text.split('/')[-1] \
                    + '_genomic.gbff.gz'
            if summary.tag == 'Synonym':
                for synonm in summary:
                    if synonm.tag == 'RefSeq':
                        ref_genome_acc = synonm.text
            if summary.tag == 'Meta':
                for meta in summary.findall('./'):
                    if meta.tag == 'Stats':
                        for stat in meta.findall('./'):
                            if stat.attrib['category'] == 'chromosome_count':
                                num_chromosomes = stat.text
                            if stat.attrib['category'] == 'total_length':
                                length = stat.text
                            if stat.attrib['category'] == 'contig_n50':
                                n50 = stat.text
                            if stat.attrib['category'] == 'contig_l50':
                                l50 = stat.text
                    if meta.tag == 'assembly-status':
                        assembly_level = meta.text
            if summary.tag == 'FtpPath_Stats_rpt':
                print(summary.text)
                genome_stats = stats_report(summary.text, genome_stats)
                num_segments = genome_stats['top-level-count']
                assembly_gc_content = genome_stats['gc-perc']
                n75 = genome_stats['scaffold-N75']
                n90 = genome_stats['scaffold-N90']

        samples[genome_assembly_id] = lineage , taxonomy_id , bco_id , schema_version \
            , analysis_platform , analysis_platform_object_id , bioproject \
            , biosample , sra_run_id , ngs_read_file_name , ngs_read_file_source \
            , organism_name , ref_genome_acc , genomic_section \
            , num_chromosomes , num_genes , num_segments  , assembly_gc_content \
            , length , n50  , n75 , n90 , l50  , l75 , number_of_n \
            , percent_assembly_greater_20x , query_coverage_against_reference \
            , percent_identity_against_reference , reads_unaligned , id_method \
            , assembly_type , assembly_level , assembly_score
    print(samples)
    return samples

def stats_report(stats_ftp, genome_stats):
    """Assembly Stats Report

    Downloads and then parses the NCBI Assembly Stats Report



    """

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
            if row[0].startswith('#'):
                continue
            genome_stats[row[6]] = [row[0], row[1], row[2], row[3], row[4],\
                row[5], row[7], row[8], row[9]]

    with open('home/assembly/stats/' + stats_file, 'r', encoding='utf8') as statfile:
        genome_tsv = csv.reader(statfile, delimiter="\t")
        for row in genome_tsv:
            if row[0].startswith('#'):
                continue
            if row[0] == 'all' and row[4] == 'gc-perc':
                genome_stats['gc-perc'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-count':
                genome_stats['top-level-count'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-N75':
                genome_stats['scaffold-N75'] = row[5]
            if row[0] == 'all' and row[4] == 'scaffold-N90':
                genome_stats['scaffold-N90'] = row[5]
    import pdb; pdb.set_trace()
    return genome_stats

def get_lineage(taxonomy_id):
    """Get lineage

    Retrieves the NCBI Taxonomic lineage given a tax id

    Parameters
    ----------
    taxonomy_id: str
        NCBI Taxonomy identifier
    """
    Entrez.email = "hadley_king@gwu.edu"
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

    samples = {}
    header = ['genome_assembly_id', 'lineage', 'taxonomy_id', 'bco_id', 'schema_version', \
        'analysis_platform', 'analysis_platform_object_id', 'bioproject', \
        'biosample',  'sra_run_id', 'ngs_read_file_name', 'ngs_read_file_source', \
        'organism_name', 'ref_genome_acc', 'genomic_section', \
        'num_chromosomes', 'num_genes', 'num_segments', 'assembly_gc_content', \
        'length', 'n50', 'n75', 'n90', 'l50', 'l75', 'number_of_n',  \
        'percent_assembly_greater_20x', 'query_coverage_against_reference', \
        'percent_identity_against_reference',  'reads_unaligned', 'id_method',\
        'assembly_type', 'assembly_level', 'assembly_score']

    args = usr_args()
    if args.input == 'api':
        argosdb_api(samples)
    else:
        parse_xml(args.input, samples)
    sample_output(samples, header, args.output)

if __name__ == "__main__":
    main()
