#!/usr/bin/env python3
import os, re, args
import argparse

from Bio import Entrez

parser = argparse.ArgumentParser(description='Extract BioSample metadata from NCBI Entrez.')

parser.add_argument('-o', '--out', required = True,
                    help='Output table presenting lookup results, if file already exists will update the file or append with new samples')
parser.add_argument('-i', '--in', required = False,
                    help='Input file of sample names')
