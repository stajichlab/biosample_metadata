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

parser.add_argument('-u', '--update', required = False, action='store_true',
                    help='Input sample names as list')

parser.add_argument('--sra', required = False, action='store_true',
                    help='Sample names are SRR IDs not BioSamples')


parser.add_argument('--debug', required = False, action='store_true',
                    help='Debug this running')


args = parser.parse_args()

Entrez.email = args.email

biosamples = {}

separator = "\t"
query = set()
header_set = set()

if args.output != "-":
    if re.search(r'\.(csv|CSV)$',args.output):
        separator = ","
    if os.path.exists(args.output):
        with open(args.output,"rt") as fh:
            csvin = csv.reader(fh,delimiter=separator)
            header = next(csvin)
            for h in header[1:]:
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

if args.sra:
    newquery = set()
    for sraid in query:
        handle = Entrez.efetch(db="sra", id=sraid)
        tree = ET.parse(handle)
        root = tree.getroot()

        for sd in root.iter('SAMPLE'):
            if args.debug:
                indent(sd)
                ET.dump(sd)
            if "accession" in sd.attrib:
                if args.debug:
                    print("SRA BioSample Accession is {}".format(sd.attrib["accession"]))
                newquery.add(sd.attrib["accession"])

#            for identifier in sd:
#                if identifier
#                for extid in identifier:
#                    if (extid.tag == "EXTERNAL_ID" and
#                        "namespace" in extid.attrib and
#                        extid.attrib["namespace"].lower() == "biosample"):
#                        newquery.add(extid.text)
    query = newquery # update the query set with this

sampmatch = re.compile(r'SAM[A-Z](\d+)')
sampidquery = set()
for qname in query:
    if args.debug:
        print("query is {}".format(qname))
    handle = Entrez.esearch(db="biosample", term=qname)
    tree = ET.parse(handle)
    root = tree.getroot()
    if args.debug:
        indent(root)
        ET.dump(root)

    for idlist in root.iter("IdList"):
        for id in idlist:
            if args.debug:
                print("ID found {}",format(id))
                indent(id)
                ET.dump(id)
            if id.tag == "Id":
                sampidquery.add(id.text)

for sampid in sampidquery:
    handle = Entrez.efetch(db="biosample", id=sampid)
    if args.debug:
        print("sampid is {}".format(sampid))
    tree = ET.parse(handle)
    root = tree.getroot()
    for sample in root:
        if args.debug:
            indent(sample)
            ET.dump(sample)
        BIOSAMPLE = sample.attrib['accession']
        if BIOSAMPLE not in biosamples:
            biosamples[BIOSAMPLE] = {}
        for attributes in root.iter('Attributes'): # sample to root?
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
    sorted_header_set = sorted(header_set)
    outheader.extend(sorted_header_set)
    outcsv.writerow(outheader)
    for sample in sorted(biosamples):
        outrow = [sample]
        for h in sorted_header_set:
            if h in biosamples[sample]:
                outrow.append(biosamples[sample][h])
            else:
                outrow.append('')
        outcsv.writerow(outrow)
