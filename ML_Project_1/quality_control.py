#!/usr/bin/env python

import csv
import sys
import os

# Load the Pandas libraries with alias 'pd'
#import pandas as pd
#from decimal import Decimal

python_path = sys.argv[0]
app_dir = os.path.dirname(os.path.realpath(python_path))
sys.path.append(os.path.abspath(app_dir))

emmax_kin = "/home/kevin.muirhead/emmax-beta-07Mar2010/emmax-kin"
emmax = "/home/kevin.muirhead/emmax-beta-07Mar2010/emmax"

#print(python_path)
#print(app_dir)

#sys.exit()

### TO DOs

# 1. Make parameter options parser for script.

# 2. Make function to gunzip the .gz files first before processing.

# 3. Split each file parser into functions using def. Make this file a quality control class for a main workflow file.
'''

return uncompressed_outfile - The uncompressed output file for parsing.
'''
def gunzip_input_files(compressed_infile):
    
    uncompressed_outfile = compressed_infile.replace("\.gz","")
    print(compressed_infile)
    print(uncompressed_outfile)
    sys.exit()
    
    os.system("gunzip -c {compressed_infile} > {uncompressed_outfile}".format(compressed_infile=compressed_infile, uncompressed_outfile=uncompressed_outfile))
    return(uncompressed_outfile)
    
'''

'''
def parse_phenotypes_file(phenotypes_infile):

    # Counter for header and entry.
    row_counter = 0

    # Phenotypes data structure.
    phenotypes_dict = {}

    # Number of phenotypes counter.
    num_phenotypes = 0
    
    # Parse the phenotype FT10.txt file.
    with open(phenotypes_infile, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        
        # Iterate over each row in the file.
        for row in csv_reader:
        
            # If the line is not the header.
            if(row_counter != 0):
                
                #print(row)

                # The genotype id of the phenotype.
                genotype_id = row[0]
                
                # The flowering time phenotype value.
                flowering_time = row[1]
                
                #print(genotype_id, flowering_time)
                
                # Filter out phenotypes that have missing data in the form of "NA" and count the number of phenotypes.
                if(flowering_time != "NA"):
                    phenotypes_dict[genotype_id] = flowering_time
                    num_phenotypes = num_phenotypes + 1
            row_counter = row_counter + 1
    csvfile.close()
    
    return(phenotypes_dict)


'''

'''
def parse_genotypes_file(genotypes_infile,phenotypes_infile,maf_threshold):

    # Get the phenotypes dictionary data structure so we can filter the genotypes by phenotype.
    phenotypes_dict = parse_phenotypes_file(phenotypes_infile)
    
#    # The output directory of the phenotypes file.
#    output_dir = os.path.dirname(phenotypes_infile)
    
    
#    # Get the output file name of the phenotypes file.
#    phenotypes_outfile = os.path.join(output_dir, "phenotypes.tsv")
#
#    # Writing the phenotypes to a TSV file for input into the RAINBOWR R program for SNP association tests.
#    tsvfile = open(phenotypes_outfile, 'w')
#    tsvwriter = csv.writer(tsvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
#
#    # Write the header for the phenotypes TSV file.
#    tsvwriter.writerow(['genotype_id', 'phenotype'])
#
    # Counter for header and data entry.
    row_counter = 0

    # Genotypes data structure.
    genotypes_dict = {}

    # Parse genotypes.csv file.
    with open(genotypes_infile, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        
        # Array for genotype ids in header line. So we can map SNPs to genotype ids.
        genotype_ids = []
        
        # Iterate over each row in the file.
        for row in csv_reader:
            if(row_counter == 0): # If row is header line.
            
                # Get the genotype_ids in index 2 (column 3) through end of row.
                genotype_ids = row[2:]
                #print(genotype_ids)
                #print(row_counter)
                #sys.exit()
            else: # Else row is data entry.
                
                # The chromosome id.
                chromosome_id = row[0]
                #print(chromosome_id)
                
                # The SNP position id.
                position_id = row[1]
                #print(position_id)
                
                # Get the alleles in index 2 (column 3) through end of row.
                alleles = row[2:]
                
                # Dictionary for each allele of the variant at a given genotype id.
                genotype_list = {}
                
                # Counter for the number of alleles.
                total_num_alleles = 0
                
                # Index for the genotype_ids and alleles arrays.
                index_counter = 0
                
                # The genotype_ids and the alleles data structures are not the same length.
                if(len(genotype_ids) != len(alleles)):
                    print("Length of genotype_ids and length of alleles are not the same.")
                    sys.exit()
                    
                # Iterate over the genotype_ids and alleles data structures while the index_counter is less than the length of the alleles data structure.
                while(index_counter < len(alleles)):
                    #print(index_counter)
                    genotype_id = genotype_ids[index_counter]
                    allele = alleles[index_counter]
                    
                    #print(genotype_id)
                    #print(allele)
                    
                    # If the genotype_id is in the phenotypes_dict then store the allele in the genotype_list.
                    # Filtering for genotype ids that have flowering time metadata.
                    if(genotype_id in phenotypes_dict):
                    
                        genotype_list[genotype_id] = allele
                        total_num_alleles = total_num_alleles + 1
                        #flowering_time = phenotypes_dict[genotype_id]
                        
                        #tsvwriter.writerow([genotype_id,flowering_time])

                        
                    index_counter = index_counter + 1
                #print(genotype_list)
            
                #sys.exit()
                #print(genotype_list)
                #print(len(genotype_list))
                #print(total_num_alleles)
                
                # Counter for counting allele bases.
                allele_base_counter = {}
                
                # Array list containing the allele bases.
                allele_base_list = []
                
                # Iterate over the genotype_list to count alleles and to check for biallelic SNPs.
                for genotype_id in genotype_list:
                    #print(genotype_id)
                    #print(genotype_list[genotype_id])
                    
                    # Initializing allele_base_counter with 0 for the first time.
                    if(not(genotype_list[genotype_id] in allele_base_counter)):
                        allele_base_counter[genotype_list[genotype_id]] = 0
                        
                    # After allele_base_counter is initialized start counting the bases.
                    if(genotype_list[genotype_id] in allele_base_counter):
                        allele_base_counter[genotype_list[genotype_id]] = allele_base_counter[genotype_list[genotype_id]] + 1
                    
                    allele_base_list.append(genotype_list[genotype_id])
                    #sys.exit()
                #print(allele_base_counter)
                #print(allele_base_list)
                
                # Obtain unique list of allele bases.
                unique_allele_set = set(allele_base_list)
                unique_allele_list = list(unique_allele_set)
                
                #print(unique_allele_list)
                
                # Filter out loci that do not contain biallellic snps.
                if (len(unique_allele_list) == 2):
                
                    # Calculate the minor and major allele frequencies and which is the minor and major allele.
                    minor_allele_count = min([allele_base_counter[unique_allele_list[0]], allele_base_counter[unique_allele_list[1]]])
                    major_allele_count = max([allele_base_counter[unique_allele_list[0]], allele_base_counter[unique_allele_list[1]]])

                    #print("Minor allele count:" + str(minor_allele_count))
                    #print("Major allele count:" + str(major_allele_count))

                    minor_allele_base = ""
                    major_allele_base = ""

                    # Figure out which base is the minor and major allele.
                    if(minor_allele_count == allele_base_counter[unique_allele_list[0]]):
                        minor_allele_base = unique_allele_list[0]
                        major_allele_base = unique_allele_list[1]

                    else:
                        minor_allele_base = unique_allele_list[1]
                        major_allele_base = unique_allele_list[0]

                    #print("Minor allele base:" + str(minor_allele_base))
                    #print("Major allele base:" + str(major_allele_base))

                    #print(total_num_alleles)
                    
                    if(total_num_alleles != (minor_allele_count + major_allele_count)):
                        print("total_num_alleles != (minor_allele_count + major_allele_count)")
                        sys.exit()
                    
                    # Calculate the minor allele frequency.
                    minor_allele_frequency = float(minor_allele_count) / float(total_num_alleles)

                    # Calculate the major allele frequency.
                    major_allele_frequency = float(major_allele_count) / float(total_num_alleles)
            
                    #print(genotype_list)
                    #sys.exit()
                    
                    # Filter minor allele frequency (MAF) >= maf_threshold.
                    if(minor_allele_frequency >= maf_threshold):
                    
                        # The genotype_metadata dictionary for storing minor allele count and base and major allele count and base.
                        genotype_metadata = {}
                        
                        # Store the minor allele count, minor allele base character, and the minor allele frequency in the genotype_metadata dictionary
                        genotype_metadata["minor_allele_count"] = minor_allele_count
                        genotype_metadata["minor_allele_base"] = minor_allele_base
                        genotype_metadata["minor_allele_frequency"] = minor_allele_frequency
                        
                        # Store the major allele count, major allele base character, and the major allele frequency in the geonotype_metadata dictionary
                        genotype_metadata["major_allele_count"] = major_allele_count
                        genotype_metadata["major_allele_base"] = major_allele_base
                        genotype_metadata["major_allele_frequency"] = major_allele_frequency
                        
                        # Store the total number of alleles in the genotype_metadata dictionary
                        genotype_metadata["total_num_alleles"] = total_num_alleles
                        
                        #print(genotype_metadata)

                        # Building the genotypes dictionary to have chromosome_id as the first key and position_id as
                        # the second key containing a dictionary of genotype ids and the corresponding SNP at that position
                        # for that genotype id.
                        if(not(chromosome_id in genotypes_dict)):
                            genotypes_dict[chromosome_id] = {}
                        if(not(position_id in genotypes_dict[chromosome_id])):
                            genotypes_dict[chromosome_id][position_id] = {}
                            genotypes_dict[chromosome_id][position_id]["genotype_list"] = genotype_list
                            genotypes_dict[chromosome_id][position_id]["genotype_metadata"] = genotype_metadata
                            
                        #print(genotypes_dict[chromosome_id][position_id]["genotype_list"])
                        #print(genotypes_dict[chromosome_id][position_id]["genotype_metadata"])
                        
                        #sys.exit()
                    
            row_counter = row_counter + 1


   
    #tsvfile.close()
    return(genotypes_dict,phenotypes_dict)
 
'''
 
'''
def generate_files_for_plink(phenotype_infile, genotypes_infile, maf_threshold, output_dir):

    # Get the genotypes dictionary data structure so we can iterate through the genotypes.
    (genotypes_dict,phenotypes_dict) = parse_genotypes_file(genotypes_infile, phenotypes_infile, maf_threshold)

    # The genotype id list data structure.
    genotype_list = []
    
    # Get the genotype ids list from the first genotypes_dict entry to sort contents of genotypes_dict at each chromosome_id and position_id pair.
    for chromosome_id in genotypes_dict:
        for position_id in genotypes_dict[chromosome_id]:
            genotype_list = genotypes_dict[chromosome_id][position_id]["genotype_list"]
            break
        break
    
    genotype_ids = sorted([int(genotype_id) for genotype_id in genotype_list.keys()])

    #print(genotype_ids)
    
    # The plink output directory.
    plink_output_dir = os.path.join(output_dir, "plink_format_files")

    # Create the plink output directory if it does not exist.
    if not os.path.exists(plink_output_dir):
        os.makedirs(plink_output_dir)
    
    genotype_ped_outfile = os.path.join(plink_output_dir, "genotypes.ped")

    #print(genotype_ped_outfile)

    # Writing the genotype scores to a TSV file for input into the rainbowr.r program for PLINK input files.
    tsvfile1 = open(genotype_ped_outfile, 'w')
    tsvwriter1 = csv.writer(tsvfile1, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    ## Write the header for the genotype counts TSV file.
    ##tsvwriter1.writerow(['marker_id'] + genotype_ids)

    genotype_map_outfile = os.path.join(plink_output_dir, "genotypes.map")
    print(genotype_map_outfile)
    # genotype map format for PLINK script
    
    # Writing the genotype map to a TSV file for input into the PLINK program for PLINK input files.
    tsvfile2 = open(genotype_map_outfile, 'w')
    tsvwriter2 = csv.writer(tsvfile2, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    ## Write the header for the genotype counts TSV file.
    ##tsvwriter2.writerow(['marker_id', 'chromosome_id', 'position_id'])


                
    # Individual genotypes dictionary for storing the genotypes of each individual for the PED formatted file for PLINK.
    individual_genotypes = {}
    
    # Count for marker ids.
    counter = 1
    for chromosome_id in genotypes_dict:
        for position_id in genotypes_dict[chromosome_id]:

            #print(chromosome_id)
            #print(position_id)
            genotype_list = genotypes_dict[chromosome_id][position_id]["genotype_list"]
            genotype_metadata = genotypes_dict[chromosome_id][position_id]["genotype_metadata"]

            #print(genotype_metadata["minor_allele_base"])
            #print(genotype_metadata["major_allele_base"])

            minor_allele_base = genotype_metadata["minor_allele_base"]
            major_allele_base = genotype_metadata["major_allele_base"]

            for genotype_id in genotype_ids:

                allele_base = genotype_list[str(genotype_id)]
                genotype = allele_base
                
                # For FarmCPU numeric format.
                # Initialize the genotype to a missing value by default.
                # Genotypes should be either 0, 1, or 2 for homozygote 1, heterozygote, and homozygote 2 genotypes.
#                encoded_genotype = "-9"
#                if(allele_base == major_allele_base):
#                    encoded_genotype = "0"
#                if(allele_base == minor_allele_base):
#                    encoded_genotype = "2"
                
                if(not(genotype_id in individual_genotypes)):
                    individual_genotypes[genotype_id] = []
                    individual_genotypes[genotype_id].append(genotype)
                    individual_genotypes[genotype_id].append(genotype)
                elif(genotype_id in individual_genotypes):
                    individual_genotypes[genotype_id].append(genotype)
                    individual_genotypes[genotype_id].append(genotype)

            marker_id = "_".join(["Chr" + chromosome_id, "Pos" + position_id])
            print(marker_id)
            
            # PLINK MAP File format.
            # chromosome (1-22, X, Y or 0 if unplaced)
            # rs# or snp identifier
            # Base-pair position (bp units)
            tsvwriter2.writerow([chromosome_id, marker_id, position_id])

    # PLINK PED File format.
    # Family ID
    # Individual ID
    # Paternal ID
    # Maternal ID
    # Sex (1=male; 2=female; other=unknown)
    # Phenotype
    for genotype_id in individual_genotypes:
        #print(genotype_id)
        
        if(str(genotype_id) in phenotypes_dict):
            phenotype = phenotypes_dict[str(genotype_id)]
            genotype_list = individual_genotypes[genotype_id]
            #if(float(phenotype) < 50):
                #tsvwriter1.writerow(["Family_ID", genotype_id, "Paternal_ID", "Maternal_ID", "other", 1] + genotype_list)
            #elif(float(phenotype) >= 50):
                            
                #tsvwriter1.writerow(["Family_ID", genotype_id, "Paternal_ID", "Maternal_ID", "other", 2] + genotype_list)
            tsvwriter1.writerow(["Family_ID", genotype_id, "Paternal_ID", "Maternal_ID", "other", phenotype] + genotype_list)
    #print(individual_genotypes)
    tsvfile1.close()
    tsvfile2.close()

    return (genotype_ped_outfile, genotype_map_outfile)



'''
Function run_SNPHWE_marker_tests(genotypes_counts_infile, output_dir)

Wrapper function for the snp_hwe_marker_tests.r script adapted from the http://csg.sph.umich.edu/abecasis/Exact/r_instruct.html website that  implements a fast exact Hardy-Weinberg Equilibrium test for SNPs as described in Wigginton, et al. (2005). The script makes use of the snp_hwe.r source file (http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r).

Input:

genotypes_counts_infile - The genotype counts TSV file obtained from the generate_genotypes_counts_file function.

snp_hwe_results_outfile - The SNPHWE marker tests output file.


'''
def run_emmax_association_tests(genotype_ped_file, genotype_map_file, emmax_output_dir):

    out_prefix = os.path.join(os.path.dirname(genotype_ped_file), "genotypes")
    #os.system("source ~/.bash_profile")
    #os.system("conda activate plink_env")
    # plink --ped ../genotypes.ped --map ../genotypes.map --recode12 --output-missing-genotype 0 --transpose --allow-no-sex --out ../genotypes
    os.system(("plink --ped {genotype_ped_file} --map {genotype_map_file} --recode12 --output-missing-genotype 0 --transpose --allow-no-sex --out {out_prefix}").format(genotype_ped_file=genotype_ped_file,genotype_map_file=genotype_map_file, out_prefix=out_prefix))
    
    genotypes_tfam_file = "".join([out_prefix, ".tfam"])
    emmax_phenotypes_infile = os.path.join(emmax_output_dir, "emmax_phenotypes.txt")

    emmax_outfile_prefix = os.path.join(emmax_output_dir, "emmax")

    os.system(("{emmax_kin} -v -d 10 {out_prefix}").format(emmax_kin=emmax_kin,out_prefix=out_prefix))
    
    os.system(("awk -F' ' '{{print $1,$2,$6}}' < {genotypes_tfam_file} > {emmax_phenotypes_infile}").format(genotypes_tfam_file=genotypes_tfam_file,emmax_phenotypes_infile=emmax_phenotypes_infile))
    
    bn_kinf_matrix = ".".join([out_prefix, "BN.kinf"])
    
    os.system(("{emmax} -v -d 10 -t {out_prefix} -p {emmax_phenotypes_infile} -k {bn_kinf_matrix} -o {emmax_outfile_prefix}").format(emmax=emmax,out_prefix=out_prefix,emmax_phenotypes_infile=emmax_phenotypes_infile,bn_kinf_matrix=bn_kinf_matrix,emmax_outfile_prefix=emmax_outfile_prefix))

def generate_genotypes_summary_file(genotypes_dict,genotypes_counts_outfile):

    # Writing the genotype counts to a TSV file for input into the snp-hwe.r program for snp exact tests.
    tsvfile = open(genotypes_counts_outfile, 'w')
    tsvwriter = csv.writer(tsvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    # Write the header for the genotype counts TSV file.
    tsvwriter.writerow(['MARKER_ID', 'p', 'q', 'HET', 'HOM1', 'HOM2'])
    #tsvwriter.writerow(['MARKER_ID', 'HET', 'HOM1', 'HOM2'])

    # Make HWE-SNP formatted file.
    for chromosome_id in genotypes_dict:
        for position_id in genotypes_dict[chromosome_id]:
            
            # Get the genotype_metadata dictionary.
            genotype_metadata = genotypes_dict[chromosome_id][position_id]["genotype_metadata"]
            
            # Get the minor allele count for this locus.
            minor_allele_count = genotype_metadata["minor_allele_count"]
            
            # Get minor allele base nucleotide character for this locus.
            minor_allele_base = genotype_metadata["minor_allele_base"]
                    
            # Get the major allele count for this locus.
            major_allele_count = genotype_metadata["major_allele_count"]
            
            # Get the major allele base nucleotide character for this locus.
            major_allele_base = genotype_metadata["major_allele_base"]
        
            # Get the total number of alleles for this locus.
            total_num_alleles = genotype_metadata["total_num_alleles"]
                    
            # Calculate the minor allele frequency.
            minor_allele_frequency = float(minor_allele_count) / float(total_num_alleles)

            # Calculate the major allele frequency.
            major_allele_frequency = float(major_allele_count) / float(total_num_alleles)

            #print(minor_allele_frequency)
            #print(major_allele_frequency)
            
            ## Calculate the Hardy Weinberg Equilibrium (HWE) equation.
            
            # The major allele frequency p.
            p = major_allele_frequency
            
            # The minor allele frequency q.
            q = minor_allele_frequency

            # Test to check that the major allele frequency p and minor allele frequency when added together equal 1. (p + q = 1)
            if((p + q) == 1):
                #print("pass")
                #print(p + q)
            
#                # The AA homozygote 1 genotype frequency (p ** 2).
#                homozygote_1_freq = (p ** 2)
#
#                # The aa homozygote 2 genotype frequency (q ** 2).
#                homozygote_2_freq = (q ** 2)
#
#                # The Aa heterozygote fequency (2pq).
#                heterozygote_freq = 2 * (p * q)

                #print(homozygote_1_freq + homozygote_2_freq + heterozygote_freq)
                
                #print(homozygote_1_freq, homozygote_2_freq, heterozygote_freq)
                
                #print(homozygote_1_freq * total_num_alleles, homozygote_2_freq * total_num_alleles, heterozygote_freq * total_num_alleles)
                
                # Generate a marker id for this locus.
                marker_id = '_'.join([':'.join(['chromosome',chromosome_id]), ':'.join(['position',position_id])])
                
#                # The AA homozygote 1 genotype count.
#                homozygote_1_count = homozygote_1_freq * total_num_alleles
#
#                # The aa homozygote 2 genotype count.
#                homozygote_2_count = homozygote_2_freq * total_num_alleles
#
#                # The Aa heterozygote count.
#                heterozygote_count = heterozygote_freq * total_num_alleles
                # The AA homozygote 1 genotype count.
                homozygote_1_count = minor_allele_frequency * total_num_alleles

                # The aa homozygote 2 genotype count.
                homozygote_2_count = major_allele_frequency * total_num_alleles

                # The Aa heterozygote count.
                heterozygote_count = 0 * total_num_alleles
                
                # Write the genotype counts entry.
                #'MARKER_ID', 'p','q','HET', 'HOM1', 'HOM2'
                tsvwriter.writerow([marker_id,p,q,heterozygote_count,homozygote_1_count,homozygote_2_count])
                #tsvwriter.writerow([marker_id,heterozygote_count,homozygote_1_count,homozygote_2_count])
                
            else:
                print("p + q ", str(p + q), " not equal to 1.")
                sys.exit()
    tsvfile.close()
    
'''
Function format_and_filter_SNPHWE_tests(snp_hwe_results_outfile):


Format and filter SNPHWE test output file for input into the generate_qqmat_plots function.

Plot the Quantile-Quantile (QQ) Plot for visualization of association tests and genotype quality (type I error rate) and Manhattan Plots to visualize significant SNP variants across the entire genome in the locations that were genotyped.

Input:

snp_hwe_results_infile -

maf_threshold - The Minor Allele Frequency threshold for filting SNP variants. MAF values greater than or equal to this value with be retained and MAF values less than this value will be filtered out of the dataset. The default for MAF is 0.05 or 5%. Default: 0.05

alpha_value - The alpha value threshold for filtering SNP variants. If the P-value for the SNPHWE tests are greater than or equal to the alpha value it is retained and if the P-value is less than the alpha value will be filtered out of the dataset. The default for alpha is 0.001 or 10^-3. Default: 0.001
.
Output:

Returns before filtering file and after filtering file for
'''


'''
Function generate_qqmat_plots():

Plot the Quantile-Quantile (QQ) Plot for visualization of association tests and genotype quality (type I error rate) and Manhattan Plots to visualize significant SNP variants across the entire genome in the locations that were genotyped.

Input:

 -

output_dir -

'''
#def generate_qqmat_plots():

'''
Function format_and_filter_genotypes(, )

Filters genotypes using the minor allele frequency (MAF) and HWE-SNP P-Value.

genotypes_counts_outfile -

genotypes_dict -
  
'''
#def format_and_filter_genotypes():


'''
Function encode_genotypes(genotypes_dict):

Depreciated
Encodes the major allele with 0 and the minor allele with 1.

Input:

genotypes_dict -

output_dir -

'''
#def encode_genotypes(genotypes_dict):





#gunzip -c gene_model.gff.gz > gene_model.gff
#gunzip -c genotype.csv.gz > genotype.csv

maf_threshold = 0.01

alpha_value = 0.001

#phenotypes_infile = "/Users/kevin.muirhead/Desktop/macbook_air/MDSC_679/ML_Project_1/FT10.txt"
phenotypes_infile = "/home/kevin.muirhead/FT10.txt"

#genotypes_infile = "/Users/kevin.muirhead/Desktop/macbook_air/MDSC_679/ML_Project_1/genotype.csv"
genotypes_infile = "/home/kevin.muirhead/genotype.csv"

#output_dir = "/Users/kevin.muirhead/Desktop/GWAS_output_dir"
output_dir = "/home/kevin.muirhead/GWAS_output_dir"

# Create the output directory if it does not exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
(genotype_ped_outfile, genotype_map_outfile) = generate_files_for_plink(phenotypes_infile, genotypes_infile, maf_threshold, output_dir)

# The emmax output directory.
emmax_output_dir = os.path.join(output_dir, "emmax_output_dir")

# Create the emmax output directory if it does not exist.
if not os.path.exists(emmax_output_dir):
    os.makedirs(emmax_output_dir)
    
run_emmax_association_tests(genotype_ped_outfile, genotype_map_outfile, emmax_output_dir)

##### NEED TO INCORPORATE THE ABOVE FUNCTIONS IN THIS FUNCTION THEN AFTER THAT ADD PHENOTYPES TO DICTIONARY.
#generate_genotypes_counts_file(genotypes_dict,genotypes_counts_outfile)

# Create the snp_hwe_results_output_dir if it does not exist.
#snp_hwe_results_output_dir = os.path.join(output_dir, "SNPHWE_results")
#if not os.path.exists(snp_hwe_results_output_dir):
#    os.makedirs(snp_hwe_results_output_dir)
#
# The path to the snp_hwe_results_outfile file.
#snp_hwe_results_outfile = os.path.join(snp_hwe_results_output_dir, "genotype_count_SNPHWE_tests.tsv")
#run_SNPHWE_marker_tests(genotypes_counts_outfile, snp_hwe_results_outfile)

#format_and_filter_SNPHWE_tests(snp_hwe_results_outfile, maf_threshold, alpha_value)
