#!/usr/bin/env python3
import os, re, csv
import argparse

from Bio import Entrez
import xml.etree.ElementTree as ET
import pprint

pp = pprint.PrettyPrinter(indent=4)

def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem


parser = argparse.ArgumentParser(description='Extract BioSample metadata from NCBI Entrez.')

parser.add_argument('-o', '--output', metavar=("out"), required = True,
                    help='Output table presenting lookup results, if file already exists will update the file or append with new samples')
parser.add_argument('-i', '--input', metavar=("in"),required = False,
                    help='Input file of sample names')

parser.add_argument('-s', '--sample', metavar=('samples'), required = False, nargs='+',
                    help='Input file of sample names')

parser.add_argument('-e', '--email', required = True,
                    help='Input your email address for Entrez queries')
parser.add_argument('-u', '--update', required = False,
                    help='Input sample names as list')

parser.add_argument('--sra', required = False,
                    help='Sample names are SRR IDs not BioSamples')


args = parser.parse_args()

Entrez.email = args.email

biosamples = {}

separator = "\t"
query = set()
header_set = set()

if args.output != "-":
    if re.match(r'\.(csv|CSV)$',args.output):
        separator = ","
    if os.path.exists(args.output):
        with open(args.output,"rt") as fh:
            csvin = csv.reader(fh,delimiter=separator)
            header = next(csvin)
            for h in header:
                header_set.add(h)
            for row in csvin:
                sampleid = row[0]
                biosamples[sampleid] = {}
                # skip column 1 since it is the ID column
                for i in range(1,len(header)):
                    biosamples[sampleid][header[i]] = row[i]


if args.input:
    with open(args.input,"rt") as fh:
        for line in fh:
            line = line.strip()
            query.add(line)

if args.sample:
    for r in args.sample:
        query.add(r)

sampmatch = re.compile(r'SAMN(\d+)')
for sampid in query:
    qname = sampid
    m = sampmatch.match(sampid)
    if m:
        qname = m.group(1)
    handle = Entrez.efetch(db="biosample", id=qname)
    tree = ET.parse(handle)
    root = tree.getroot()

    for sample in root:
        BIOSAMPLE = sample.attrib['accession']
        if BIOSAMPLE not in biosamples:
            biosamples[BIOSAMPLE] = {}
        for attributes in root.iter('Attributes'):
            for metadata in attributes:
                keyname = metadata.attrib['attribute_name']
                if 'harmonized_name' in metadata.attrib:
                    keyname = metadata.attrib['harmonized_name']
                header_set.add(keyname)
                biosamples[BIOSAMPLE][keyname] = metadata.text
#                indent(metadata)
#                ET.dump(metadata)
#    records = Entrez.parse(handle)
#    for record in records:
#        print(record)

with open(args.output, 'wt') as outfh:
    outcsv = csv.writer(outfh, delimiter=separator)
    outheader = ['BioSample']
    outheader.extend(sorted(header_set))
    outcsv.writerow(outheader)
    for sample in sorted(biosamples):
        outrow = [sample]
        for h in header_set:
            if h in biosamples[sample]:
                outrow.append(biosamples[sample][h])
            else:
                outrow.append('')
        outcsv.writerow(outrow)
