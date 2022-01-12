#!/usr/bin/env python3
""" SRA Run Analysis Parser

This script will parse the XML result file from the SRA Trace site(thehome of sra runs). It takes a
directory as an input and will iterate through all files, parsing each one to first a dictionary, and then
outputting to the terminal in TSV format, or to a file if an output file is provided. 

"""

import os
import sys
import csv
import requests
import argparse
import xml.etree.ElementTree as ET
__version__ = "0.1.0"
__status__ = "Beta"

samples = {}
header = ['Run (sra_run_id)', 'Experiment (sra_experiment_id)', 'Project (sra_project_id)', 'Biosample (sra_biosample_id)', \
    '# of bases (num_of_bases)', 'Size (file_size)', 'Published', 'Source (source)', 'Strategy (strategy)', 'Layout (layout)', \
    'Library Name', 'Library Selection', 'Instrument (instrument)', 'File type (file_type)', \
        'Unidentified reads (unidentified_reads)', 'Identified reads (identified_reads)', 'Organism of interest', \
            'Lineage for organism of interest', '% organism of interest', 'G/C Content']

def usr_args():
    """Program Arguments

    All arguments for process are defined here.
    """
    
        # initialize parser
    parser = argparse.ArgumentParser()

    # set usages options
    parser = argparse.ArgumentParser(
        prog='sra_trace_xml.py',
        usage='%(prog)s [options]')

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)

    parser.add_argument('-d', '--directory',
        required=True,
        help="Directory for list of input files. The input files should be SRR XML files.")
    parser.add_argument('-o', '--output',
        help="Output file. If no output is provided, the resulting table will be output to the terminal.")
    if len(sys.argv) <= 1:
        sys.argv.append('--help')
    
    return parser.parse_args()

def parseXML( xmlFile, samples):
    """Parse XML file
    
    Parameters
    ----------
    xmlFile: str
        file path/name to be parsed
    samples: dict
        dictionary of samples
    
    Returns
    -------
    samples: dict
        dictionary of samples
    
    """

    # print(file)
    base_total, base_a, base_c, base_g, base_t, base_n = '-', '-', '-', '-', '-', '-'
    sra_experiment_id = sra_run_id = sra_project_id = sra_biosample_id = num_of_bases = file_size \
        = published = source = strategy = layout = library_name = selection = instrument = file_type \
        = unidentified_reads = identified_reads = tax_id = lineage = percentIdentified = gc_content = '-'
    
    # create element tree object and get root element
    try:
        root = ET.parse(xmlFile).getroot()
        sra_run_id = root.attrib['run']
        for run in root.findall('./RUN'):
            num_of_bases = run.attrib['total_bases']
            file_size = run.attrib['size']
            published = run.attrib['published']
            file_type = run.findall('./SRAFiles/')[0].attrib['semantic_name']
            for run_tag in run.findall('./'):
                if run_tag.tag == 'EXPERIMENT_REF':
                    sra_experiment_id = run_tag.attrib['accession']
                if run_tag.tag == 'Bases':
                    base_total = int(run_tag.attrib['count'])
                    for base in run_tag:
                        if base.attrib['value'] == 'A':
                            base_a = int(base.attrib['count'])
                        if base.attrib['value'] == 'C':
                            base_c = int(base.attrib['count'])
                        if base.attrib['value'] == 'G':
                            base_g = int(base.attrib['count'])
                        if base.attrib['value'] == 'T':
                            base_t = int(base.attrib['count'])
                        if base.attrib['value'] == 'N':
                            base_n = int(base.attrib['count'])
                if run_tag.tag == 'tax_analysis':
                    total_spot_count = int(run_tag.attrib['total_spot_count'])
                    analyzed_spot_count = int(run_tag.attrib['analyzed_spot_count'])
                    try:
                        identified_reads = int(run_tag.attrib['identified_spot_count'])
                    except KeyError:
                        identified_reads = int(run_tag[0].attrib['total_count'])
                if run_tag.tag == 'Lineage':
                    for line in run_tag:
                        if lineage == '-':
                            lineage = str(line.text)
                        else:
                            lineage = lineage + '|'+ line.text
                        # lineage.append(line.text)
                        if line.attrib != {}: 
                            tax_id = line.attrib['taxid']
        for experiment in root.findall('./EXPERIMENT'):
            for platform in experiment.findall('./PLATFORM/'):
                for instriment in platform:
                    instrument = instriment.text
            for descriptor in experiment.findall('./DESIGN/LIBRARY_DESCRIPTOR'):
                for lib_name in descriptor.findall('./LIBRARY_NAME'): 
                    library_name = lib_name.text
                for lib_stratagy in descriptor.findall('./LIBRARY_STRATEGY'): 
                    strategy = lib_stratagy.text
                for lib_source in descriptor.findall('./LIBRARY_SOURCE'): 
                    source = lib_source.text
                for lib_selection in descriptor.findall('./LIBRARY_SELECTION'): 
                    selection = lib_selection.text
                for lib_layout in descriptor.findall('./LIBRARY_LAYOUT/'): 
                    layout = lib_layout.tag
        for biosample in root.findall('./SAMPLE/'):
            if biosample.tag == 'SAMPLE_NAME' and lineage == '-': 
                lineage = biosample[1].text
            if biosample.tag == 'IDENTIFIERS':
                sra_biosample_id = biosample[1].text
        for study in root.findall('./STUDY'):
            sra_project_id = study.attrib['alias']
        
        try:
            # print('Math check', base_total == int(num_of_bases))
            gc_content = float(base_g + base_c)/float(base_total)
            unidentified_reads = analyzed_spot_count - identified_reads
            percentIdentified = float(identified_reads)/float(analyzed_spot_count)
        except ValueError:
            pass
        # print(base_total, base_a, base_c, base_g, base_t, base_n, gc_content, lineage)
        samples[sra_run_id] = [sra_experiment_id, sra_project_id, sra_biosample_id, num_of_bases, \
            file_size, published, source, strategy, layout,library_name, selection, instrument, file_type, unidentified_reads, \
            identified_reads, tax_id, lineage, percentIdentified, gc_content]

    except ET.ParseError:
        print(xmlFile, 'not well-formed')
    
    return samples
    

def sampleOutput( samples, header, output):
    """Sample Output

    If an output file is supplied in the user arguments the samples dictionary will be output in TSV to the supplied file.
    If no output is provided, the samples dictionary will print to the terminal in a TSV format. 

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
        with open(sample_file, 'w') as file:
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
    args = usr_args()
    for item in os.listdir(args.directory):
        xmlFile = os.path.join(args.directory, item)
        parseXML(xmlFile, samples)
    sampleOutput(samples, header, args.output)
      
if __name__ == "__main__":
    main()