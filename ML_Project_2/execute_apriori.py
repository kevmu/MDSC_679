'''

Name: Kevin Muirhead
UCID#: 00502756

execute_apriori.py - Executes the Apriori Algorithm, calculates association rule metrics
and prints the association rule metrics support, confidence, lift, leverage and conviction.

'''
import sys
import os
import argparse

import re
import csv

# Getting the path of the execute_apriori.py script.
python_path = sys.argv[0]
app_dir = os.path.dirname(os.path.realpath(python_path))
sys.path.append(os.path.abspath(app_dir))

# Import the Apriori Class.
from Apriori import *

parser = argparse.ArgumentParser()

database_infile = None
min_support_count = 2
min_confidence = 0.60
output_dir = None

# Usage example
# python $HOME/software/MDSC_679/ML_Project_2/execute_apriori.py --input_file $HOME/software/MDSC_679/ML_Project_2/test_datasets/test_dataset1/transaction_database1.tsv --min_support_count 2 --min_confidence 0.60 --output_dir  $HOME/software/MDSC_679/ML_Project_2/test_dataset1_output_dir

parser = argparse.ArgumentParser(description='Perform the AprioriTID algorithm using the Apriori.py class. Executes the Apriori Algorithm, calculates association rule metrics and prints the association rule metrics support, confidence, lift, leverage and conviction.' )

parser.add_argument('--input_file', action='store', dest='database_infile',
                    help='The transaction database file as input. (i.e. $HOME/filename.tsv)')
parser.add_argument('--min_support_count', action='store', dest='min_support_count',
                    help='The minimum support count as input. Can be any integer value greater than zero. min_support_count > 0. Default: 2')
parser.add_argument('--min_confidence', action='store', dest='min_confidence',
                    help='The minimum confidence as input. Can be in the range of 0.00-1.00. Default: 0.60')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='The output directory as input. (i.e. $HOME/output_dir)')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

database_infile = results.database_infile
min_support_count = results.min_support_count
min_confidence = results.min_confidence
output_dir = results.output_dir

if(database_infile == None):
    print('\n')
    print('error: please use the -i option to specify the transaction database file as input. (i.e. filename.tsv)')
    print('database_infile =' + ' ' + str(database_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(min_support_count == None):
    print('\n')
    print('error: please use the --min_support_count option to specify the minimum support count as input. Can be any integer value greater than zero. min_support_count > 0. Default: 2')
    print('database_infile =' + ' ' + str(min_support_count))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(min_confidence == None):
    print('\n')
    print('error: please use the --min_confidence option to specify the minimum confidence as input. Can be in the range of 0.00-1.00. Default: 0.60')
    print('database_infile =' + ' ' + str(min_confidence))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

# Create the output_dir directory if it does not already exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Get the basename of the transaction database input file.
basename = os.path.basename(database_infile)

# Get the filename without the extension so we can make the association rule metrics output filename.
filename = os.path.splitext(basename)[0]

# The association rule metrics output file.
association_rule_metrics_outfile = os.path.join(output_dir, "_".join([filename, "association_rule_metrics.tsv"]))

# Run the Apriori agorithm and print out the association rule metrics to a file.
Apriori(database_infile, association_rule_metrics_outfile, min_support_count, min_confidence)
