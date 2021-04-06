import sys
import os
import argparse

import re
import csv

python_path = sys.argv[0]
app_dir = os.path.dirname(os.path.realpath(python_path))
sys.path.append(os.path.abspath(app_dir))

from Apriori import *

parser = argparse.ArgumentParser()

database_infile = None
output_dir = None

parser.add_argument('-i', action='store', dest='database_infile',
                    help='transaction database file as input. (i.e. filename.tsv)')
parser.add_argument('-o', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

database_infile = results.database_infile
output_dir = results.output_dir

if(database_infile == None):
    print('\n')
    print('error: please use the -i option to specify the transaction database file as input. (i.e. filename.tsv)')
    print('database_infile =' + ' ' + str(database_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the -o option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)


if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get the basename of the transaction database input file.
basename = os.path.basename(database_infile)

# Get the filename without the extension so we can make the association rule metrics output filename.
filename = os.path.splitext(basename)[0]

# The association rule metrics output file.
association_rule_metrics_outfile = os.path.join(output_dir, "_".join([filename, "association_rule_metrics.tsv"]))

# Minimum support count is 2.
min_support_count = 2

# Minimum confidence is 60 %.
min_confidence = 0.60

Apriori(database_infile, association_rule_metrics_outfile, min_support_count, min_confidence)
