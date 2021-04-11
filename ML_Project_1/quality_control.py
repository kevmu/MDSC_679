#!/usr/bin/env python

import csv
import sys
import os
import argparse

### Sample command for quality_filtering.py python script. Make sure that the ML_Project_1_env conda environment is activated.
# conda activate ML_Project_1_env
# python /Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/quality_control.py --phenotypes_infile /Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/INPUT_FILES/FT10.txt --genotypes_infile /Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/INPUT_FILES/genotype.csv.gz --gff_infile /Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/INPUT_FILES/gene_model.gff.gz --alpha_value 0.05 --maf_threshold 0.01 --output_dir /Users/kevin.muirhead/Desktop/GWAS_OUTPUT_DIR

# Get the path of the script so that we can use the R scripts in the directory.
python_path = sys.argv[0]
app_dir = os.path.dirname(os.path.realpath(python_path))
sys.path.append(os.path.abspath(app_dir))

parser = argparse.ArgumentParser()

phenotype_infile = None
genotypes_infile = None
gff_infile = None
maf_threshold = None
alpha_value = None
output_dir = None

parser = argparse.ArgumentParser(description='Perform quality filtering of the phenotypes infile and genotypes infile of missing values, biallecic SNPs, minor allele frequency (MAF), association mapping of filtered variants, visualization using QQ-Plots, Manhattan Plots, SNP density plots, and quality filtering of SNP variant using the alpha_value. Make sure that the ML_Project_1_env conda environment is activated. conda activate ML_Project_1_env')
parser.add_argument('--phenotypes_infile', action='store', dest='phenotypes_infile',
                    help='The phenotypes file as input. (i.e. FT10.txt)')
parser.add_argument('--genotypes_infile', action='store', dest='genotypes_infile',
                    help='The genotypes file as input. (i.e. genotype.csv.gz)')
parser.add_argument('--gff_infile', action='store', dest='gff_infile',
                    help='The gene model annotation GFF format file as input. (i.e. gene_model.gff.gz)')
parser.add_argument('--maf_threshold', action='store', dest='maf_threshold',
                    help='The minor allele frequency (MAF) threshold for filtering genotypes. (i.e. Default: 0.01)')
parser.add_argument('--alpha_value', action='store', dest='alpha_value',
                    help='The alpha value as input for filtering adjusted pvalues of the rMVP association tests by the False Disovery Rate (FDR) (i.e. Default: 0.05)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME/output_dir)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

# Parse parameter options results.
results = parser.parse_args()

# Parse the phenotypes input file.
phenotypes_infile = results.phenotypes_infile

# Parse the genotypes input file.
genotypes_infile = results.genotypes_infile

# Parse the gff annotation input file.
gff_infile = results.gff_infile

# Parse the minor allele frequency (MAF). Make sure that the maf threshold is a float.
maf_threshold = float(results.maf_threshold)

# Parse the alpha value. Make sure that the alpha value is a float.
alpha_value = float(results.alpha_value)

# Parse the output directory.
output_dir = results.output_dir

if(phenotypes_infile == None):
    print('\n')
    print('error: please use the --phenotypes_infile option to specify the phenotypes file. (i.e. FT10.txt)')
    print('phenotypes_infile =' + ' ' + str(phenotypes_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(genotypes_infile == None):
    print('\n')
    print('error: please use the --genotypes_infile option to specify the genotype file. (i.e. genotype.csv.gz)')
    print('genotypes_infile =' + ' ' + str(genotypes_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(gff_infile == None):
    print('\n')
    print('error: please use the --gff_infile option to specify the gene model annotation GFF format file as input. (i.e. gene_model.gff.gz)')
    print('gff_infile =' + ' ' + str(gff_infile))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(maf_threshold == None):
    print('\n')
    print('error: please use the --maf_threshold option to specify the minor allele frequency (MAF) threshold for filtering genotypes. (i.e. Default: 0.01)')
    print('maf_threshold =' + ' ' + str(maf_threshold))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(alpha_value == None):
    print('\n')
    print('error: please use the --alpha_value option to specify the alpha value as input for filtering adjusted pvalues of the rMVP association tests by the False Disovery Rate (FDR) ')
    print('alpha_value =' + ' ' + str(alpha_value))
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

# Create the output_dir directory if it does not exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



'''

Uncompress the compressed inout file and return the file path of the uncompressed file.

Input:

    compressed_infile - The compressed file for input.

    output_dir - The outpue directory for writting the uncompressed file.

Output:

    uncompressed_outfile - The uncompressed output file for parsing.
    
'''
def gunzip_input_files(compressed_infile, output_dir):
    
    #The basename of the input file.
    basename = os.path.basename(compressed_infile)
    
    # The filename.
    filename = os.path.splitext(basename)[0]
    
    #print(filename)
    uncompressed_outfile = os.path.join(output_dir, filename)
   
    # Run the gunzip -c command to uncompress the gzipped file.
    os.system("gunzip -c {compressed_infile} > {uncompressed_outfile}".format(compressed_infile=compressed_infile, uncompressed_outfile=uncompressed_outfile))

    return(uncompressed_outfile)
    
'''

Function parse_phenotypes_file(phenotypes_infile)

Parse the phenotypes file, filter out missing data and return the phenotype dictionary.

Input:

    phenotypes_infile - The phenotypes input file.

Output:

    phenotype_dict - The filtered phenotypes dictionary data structure where keys are genotype_ids and the value is the phenotype.

'''
def parse_phenotypes_file(phenotypes_infile):

    # If the phenotypes_infile is a compressed file then uncompress the file and return the file path of the uncompressed file.
    if(".gz" in phenotypes_infile):
        phenotypes_infile = gunzip_input_files(phenotypes_infile, output_dir)

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
            
    # Close the phenotype input file.
    csvfile.close()
    
    return(phenotypes_dict)


'''

Function parse_genotypes_file(genotypes_infile,phenotypes_infile,maf_threshold)

Parse the genotypes file using the phenotypes file, filter using the minor allele frequency (MAF) threshold and output the genotypes and phenotypes dictionary data structures.

Input:

    genotypes_infile - The genotypes input file.

    phenotypes_infile - The phenotypes input file.

    maf_threshold - The Minor Allele Frequency threshold for filtering SNP variants. MAF values greater than or equal to this value will be retained and MAF values less than this value will be filtered out of the dataset.

Output:

    phenotype_dict - The filtered phenotypes dictionary data structure where keys are genotype_ids and the value is the phenotype.

    genotypes_dict - The filtered genotypes dictionary data structure where keys are chromosome_id and position_id of each SNP and the value is the genotype_list and genotype_metadata dictionaries. The genotype list contains the genotype information of which genotype is at the position and the genotype metadata contains the minor/major allele counts, frequency, and the allele bases. As well as the total number of alleles.

'''
def parse_genotypes_file(genotypes_infile,phenotypes_infile,maf_threshold):

    # Get the phenotypes dictionary data structure so we can filter the genotypes by phenotype.
    phenotypes_dict = parse_phenotypes_file(phenotypes_infile)
    
    # If the genotypes_infile is a compressed file then uncompress the file and return the file path of the uncompressed file.
    if(".gz" in genotypes_infile):
        genotypes_infile = gunzip_input_files(genotypes_infile, output_dir)


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
   
    return(genotypes_dict,phenotypes_dict)
 
'''
Function generate_files_for_association_mapping(phenotype_dict, genotypes_dict, maf_threshold, output_dir)

Generates the plink Pedigree (PED) format and plink genotype map (MAP) format files and the phenotypes infile for the rMVP R script.

Input:

    phenotype_dict - The filtered phenotypes dictionary data structure where keys are genotype_ids and the value is the phenotype.

    genotypes_dict - The filtered genotypes dictionary data structure where keys are chromosome_id and position_id of each SNP and the value is the genotype_list and genotype_metadata dictionaries. The genotype list contains the genotype information of which genotype is at the position and the genotype metadata contains the minor/major allele counts, frequency, and the allele bases. As well as the total number of alleles.

     output_dir - The output directory to write the output files.
 
 Output:
 
     plink_genotype_ped_outfile - The Plink Pedigree (PED) format file.
     
     plink_genotype_map_outfile - The Plink Pedigree (PED) format file.
     
     mvp_phenotype_outfile - The rMVP phenotypes file.

'''
def generate_files_for_association_mapping(phenotype_dict, genotypes_dict, output_dir):

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
    
    # The association mapping output directory.
    association_mapping_output_dir = os.path.join(output_dir, "FILES_FOR_ASSOCIATION_MAPPING")

    # Create the plink output directory if it does not exist.
    if not os.path.exists(association_mapping_output_dir):
        os.makedirs(association_mapping_output_dir)
        
    # The plink genotypes ped output file.
    plink_genotype_ped_outfile = os.path.join(association_mapping_output_dir, "plink.genotype.ped.txt")

    # Writing the genotype scores to a TSV file for input into the plink program.
    tsvfile1 = open(plink_genotype_ped_outfile, 'w')
    tsvwriter1 = csv.writer(tsvfile1, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    # The plink genotypes map output file.
    plink_genotype_map_outfile = os.path.join(association_mapping_output_dir, "plink.genotype.map.txt")
    print(plink_genotype_map_outfile)
    
    # genotype map format for PLINK script
    # Writing the genotype map to a TSV file for input into the PLINK program for PLINK input files.
    tsvfile2 = open(plink_genotype_map_outfile, 'w')
    tsvwriter2 = csv.writer(tsvfile2, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    # The mvp phenotype output file.
    mvp_phenotype_outfile = os.path.join(association_mapping_output_dir, "mvp.phenotype.txt")
    print(mvp_phenotype_outfile)
    
    # Writing the mvp phenotype to a TSV file for input into the mvp program.
    tsvfile3 = open(mvp_phenotype_outfile, 'w')
    tsvwriter3 = csv.writer(tsvfile3, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    ## Write the header for the mvp phenotype TSV file.
    tsvwriter3.writerow(['genotype_id', 'phenotype'])

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
                
                    
                # Append the genotype twice to each genotype id. The plink format takes the genotype AA and delimits each allele of the genotype with a tab.
                if(not(genotype_id in individual_genotypes)):
                    individual_genotypes[genotype_id] = []
                    individual_genotypes[genotype_id].append(genotype)
                    individual_genotypes[genotype_id].append(genotype)
                elif(genotype_id in individual_genotypes):
                    individual_genotypes[genotype_id].append(genotype)
                    individual_genotypes[genotype_id].append(genotype)

            # Create the marker id for this SNP using the chromosome_id and position_id as unique identifiers.
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
    # Sex (1=male; 2=female; other=unknown)x
    # Phenotype
    for genotype_id in individual_genotypes:
        #print(genotype_id)
        
        if(str(genotype_id) in phenotypes_dict):
            phenotype = phenotypes_dict[str(genotype_id)]
            genotype_list = individual_genotypes[genotype_id]
            
            tsvwriter1.writerow(["Family_ID", genotype_id, "Paternal_ID", "Maternal_ID", "other", phenotype] + genotype_list)
    #print(individual_genotypes)
    
    # Iterate over the individual_genotypes list and output the phenotypes to the mvp phenotypes output file.
    for genotype_id in individual_genotypes:
        #print(genotype_id)
        tsvwriter3.writerow([str(genotype_id), phenotypes_dict[str(genotype_id)]])
    
    # Close the PLINK PED format file.
    tsvfile1.close()
    
    # Close the PLINK MAP format file.
    tsvfile2.close()
    
    # Close the phenotypes file.
    tsvfile3.close()
    
    return (plink_genotype_ped_outfile, plink_genotype_map_outfile, mvp_phenotype_outfile)

'''
Function run_mvp_association_tests(plink_genotype_ped_infile, plink_genotype_map_infile, mvp_phenotype_infile, association_mapping_output_dir)

Wrapper function for running plink to generate VCF files for the rMVP R script (rMVP_marker_tests.R) that performs a Generalized Linear Model (GLM), Mixed Linear Model (MLM) and the FarmCPU Model association tests on each SNP variant for each genotype. Plots the Quantile-Quantile (QQ) Plot for visualization of association tests and genotype quality (type I error rate) and Manhattan Plots to visualize significant SNP variants across the entire genome in the locations that were genotyped. Also generates SNP density plots and many other types of plots. Then runs the calculate_adjusted_pvalues.R script for calculating the Bonferroni correction and q-value (FDR) adjusted pvalues and returns the adjusted pvalues file as output.


Input:

    plink_genotype_ped_infile - The Pedigree format (PED) input file for plink.

    plink_genotype_map_infile - The genotype map format (MAP) input file for plink.

    mvp_phenotype_infile - The phenotype list input file for rMVP.

    association_mapping_output_dir - The association mapping output directory for plink and rMVP output files.


Output:

    adjusted_pvalues_outfile - The adjusted pvalues output file.

'''
def run_mvp_association_tests(plink_genotype_ped_infile, plink_genotype_map_infile, mvp_phenotype_infile, association_mapping_output_dir):

    # The output file prefix for plink command.
    out_prefix = os.path.join(os.path.dirname(plink_genotype_ped_infile), "mvp.genotype")
    
    # Execute the plink command and recode for the VCF file format and PCA output files.
    # plink --ped plink.genotype.ped.txt --map plink.genotype.map.txt --pca --recode vcf  --allow-no-sex --out mvp.genotype
    os.system(("plink --ped {plink_genotype_ped_infile} --map {plink_genotype_map_infile} --pca --recode vcf  --allow-no-sex --out {out_prefix}").format(plink_genotype_ped_infile=plink_genotype_ped_infile,plink_genotype_map_infile=plink_genotype_map_infile, out_prefix=out_prefix))
    
    # The original VCF file.
    mvp_genotypes_vcf_file = "".join([out_prefix, ".vcf"])
    
    # The fixed VCF file.
    mvp_fixed_genotypes_vcf_file = "".join([out_prefix, ".fixed.vcf"])
    
    # Have to remove the "Family_ID_" code from the vcf file generated from plink as it interferes with the rMVP script.
    os.system("sed 's/Family_ID_//g' < {mvp_genotypes_vcf_file} > {mvp_fixed_genotypes_vcf_file}".format(mvp_genotypes_vcf_file=mvp_genotypes_vcf_file, mvp_fixed_genotypes_vcf_file=mvp_fixed_genotypes_vcf_file))
    
    # Run the rMVP_marker_tests.R for performing the GLM, MLM, and FarmCPU association tests.
    # Rscript rMVP_marker_tests.R -i ~/Desktop/GWAS_output_dir2/FILES_FOR_ASSOCIATION_MAPPING/mvp.genotype.fixed.vcf -p ~/Desktop/GWAS_output_dir2/FILES_FOR_ASSOCIATION_MAPPING/mvp.phenotype.txt -o ~/Desktop/GWAS_output_dir2/ASSOCIATION_MAPPING_OUTPUT_DIR
    os.system(("Rscript rMVP_marker_tests.R -i {mvp_fixed_genotypes_vcf_file} -p {mvp_phenotype_infile} -o {association_mapping_output_dir}").format(mvp_fixed_genotypes_vcf_file=mvp_fixed_genotypes_vcf_file, mvp_phenotype_infile=mvp_phenotype_infile, association_mapping_output_dir=association_mapping_output_dir))
    
    # The mvp phenotype association test results Mixed Linear Model (MLM) file.
    mvp_phenotype_association_file = os.path.join(association_mapping_output_dir, "phenotype.MLM.csv")
    
    # The adjusted pvalues output file.
    adjusted_pvalues_outfile = os.path.join(association_mapping_output_dir, "phenotype.MLM.adjusted.pvalues.tsv")
        
    # Run the calculate_adjusted_pvalues.R script to calculate the boneferroni correction and qvalue (FDR) adjusted pvalues.
    os.system("Rscript calculate_adjusted_pvalues.R -i {mvp_phenotype_association_file} -o {adjusted_pvalues_outfile} ".format(mvp_phenotype_association_file=mvp_phenotype_association_file, adjusted_pvalues_outfile=adjusted_pvalues_outfile))
    
    return(adjusted_pvalues_outfile)

    
'''

Function parse_multiple_corrected_tests(adjusted_pvalues_infile,genotypes_dict,phenotypes_dict,alpha_value, parsed_genotypes_output_dir)

Parses the genotypes from the multiple corrected tests for adjusted pvalues using the alpha value as the cut-off. Retains qvalues <= 0.05 alpha_value for a 5% False discovery Rate (FDR). Outputs the encoded genotypes and the apriori transaction genotypes database files.

Input:

    adjusted_pvalues_infile - The adjusted pvalues input file.

    phenotype_dict - The filtered phenotypes dictionary data structure where keys are genotype_ids and the value is the phenotype.

    genotypes_dict - The filtered genotypes dictionary data structure where keys are chromosome_id and position_id of each SNP and the value is the genotype_list and genotype_metadata dictionaries. The genotype list contains the genotype information of which genotype is at the position and the genotype metadata contains the minor/major allele counts, frequency, and the allele bases. As well as the total number of alleles.

    alpha_value - The alpha value threshold for filtering SNP variants. If the P-value for the assocation test are greater than or equal to the alpha value it is retained and if the P-value is less than the alpha value will be filtered out of the dataset.

    parsed_genotypes_output_dir - The parsed genotypes output directory for writing output files.

'''
def parse_multiple_corrected_tests(adjusted_pvalues_infile,genotypes_dict,phenotypes_dict,alpha_value, parsed_genotypes_output_dir):

    # The genotype id list data structure.
    genotype_list = []
    
    # Get the genotype ids list from the first genotypes_dict entry to sort contents of genotypes_dict at each chromosome_id and position_id pair.
    for chromosome_id in genotypes_dict:
        for position_id in genotypes_dict[chromosome_id]:
            genotype_list = genotypes_dict[chromosome_id][position_id]["genotype_list"]
            break
        break
    
    genotype_ids = sorted([int(genotype_id) for genotype_id in genotype_list.keys()])
    
    # Filtered genotypes data structure.
    filtered_genotypes_dict = {}

    # Counter for header and row entries.
    row_counter = 0
    
    # Parse adjusted pvalues file.
    with open(adjusted_pvalues_infile, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        
        # Iterate over each row in the file.
        for row in csv_reader:
            print(row)
            # If the current row is not the header line.
            if(row_counter != 0):
                    
                #print(row)
                #print(row_counter)
                (marker_id,filtered_chromosome_id,filtered_position_id,ref_allele,alt_allele,effect,se,p_value,bonf_corr_pvalue,q_value) = row
                
                print(marker_id,filtered_chromosome_id,filtered_position_id,ref_allele,alt_allele,effect,se,p_value,bonf_corr_pvalue,q_value)
                

                # Split the marker_id into the chromosome_id and position_id components for accessing the genotypes_dict contents.
                # Chr5_Pos26963862
                chromosome_id = marker_id.split("_")[0].replace("Chr","")
                position_id = marker_id.split("_")[1].replace("Pos","")
                
                # If the p_value is less than or equal to the alpha_value.
                #print(float(p_value) <= float(alpha_value))
                if(float(p_value) <= float(alpha_value)):
                    filtered_genotypes_dict[marker_id] = genotypes_dict[chromosome_id][position_id]
                
                    #print(genotypes_dict[chromosome_id][position_id])
                
                    #sys.exit()
            row_counter = row_counter + 1
    
    
    #print(len(filtered_genotypes_dict))
    #sys.exit()
    
   
    # The encoded genotypes output file.
    encoded_genotypes_outfile = os.path.join(parsed_genotypes_output_dir, "encoded_genotypes.tsv")
        
    # Writing the encoded genotypes to a TSV file.
    tsvfile1 = open(encoded_genotypes_outfile, 'w')
    tsvwriter1 = csv.writer(tsvfile1, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    #filtered_genotypes_dict[marker_id]
    
    marker_ids = sorted([str(marker_id) for marker_id in filtered_genotypes_dict.keys()])
        
    # Write the header for the encoded genotypes TSV file.
    tsvwriter1.writerow(['genotype_id'] + marker_ids + ['phenotype'])

    # The apriori_genotype_pattern_output_dir output directory.
    apriori_genotype_pattern_output_dir = os.path.join(parsed_genotypes_output_dir, "apriori_genotype_pattern_files")

    # Create the apriori_genotype_pattern output directory if it does not exist.
    if not os.path.exists(apriori_genotype_pattern_output_dir):
        os.makedirs(apriori_genotype_pattern_output_dir)
    
    # Number of items per chromosome per genotype for Apriori algorithm.
    num_items = 24
    
    # The apriori genotypes file part number. There will be num_genotypes/num_items files. There will be at most num_items per file. Since patterns don't really matter the number of items per "transaction database" doesn't matter as we are trying to find the most frequent item sets.
    part_num = 1
    
    # The apriori input file data structure for organizing data in the proper transaction database format.
    apriori_genotype_files_dict = {}
    apriori_genotype_files_dict[str(part_num)] = {}
    
    # Row counter for filtered_genotypes_dict iteration.
    row_counter = 1
    
    # The encoded individual genotypes list. For writting genotypes to a file.
    encoded_individual_genotypes = {}
        
    # Iterate over the sorted marker_ids list so that we can make files for training, testing and evaluating our models and to test the Apriori algorithm on a genomic dataset.
    for marker_id in marker_ids:
        print(marker_id)
        genotype_list = filtered_genotypes_dict[marker_id]["genotype_list"]
        genotype_metadata = filtered_genotypes_dict[marker_id]["genotype_metadata"]
        
        print(genotype_list)
    
        # If the number of rows equals a multiple of the num_items. This is for organizing the apriori genotype transaction files. Have to partition files by number of items.
        if((row_counter % num_items) == 0):
            part_num = part_num + 1
            apriori_genotype_files_dict[str(part_num)] = {}

            #sys.exit()
        
        # Iterate over the sorted genotype ids list so that all genotypes for each row are the same and sorted.
        for genotype_id in genotype_ids:
        
            # Get the SNP variant of the current genotype id.
            allele_base = genotype_list[str(genotype_id)]
            print(genotype_metadata)
            #sys.exit()
            
            # The minor allele base for this marker_id.
            minor_allele_base = genotype_metadata["minor_allele_base"]
            
            # The major allele base for this marker_id.
            major_allele_base = genotype_metadata["major_allele_base"]
                        
            print(allele_base + " == " + major_allele_base)
            print(allele_base + " == " + minor_allele_base)
            
            # Compare the SNP variant of the current genotype id to the major and minor alleles. If SNP variant is the major allele then encode as 0 and if the SNP variant is the minor allele then encode as 1.
            if(allele_base == major_allele_base):
                encoded_genotype = 0
                if(not(str(genotype_id) in encoded_individual_genotypes)):
                    encoded_individual_genotypes[str(genotype_id)] = {}
                if(not(str(marker_id) in encoded_individual_genotypes[str(genotype_id)])):
                    encoded_individual_genotypes[str(genotype_id)][str(marker_id)] = encoded_genotype
                else:
                    encoded_individual_genotypes[str(genotype_id)][str(marker_id)] = encoded_genotype
                            
            else:
                encoded_genotype = 1
                if(not(str(genotype_id) in encoded_individual_genotypes)):
                    encoded_individual_genotypes[str(genotype_id)] = {}
                if(not(str(marker_id) in encoded_individual_genotypes[str(genotype_id)])):
                    encoded_individual_genotypes[str(genotype_id)][str(marker_id)] = encoded_genotype
                else:
                    encoded_individual_genotypes[str(genotype_id)][str(marker_id)] = encoded_genotype
                
            # Genotype item pattern for genotype_id from the marker_id and allele base for generating apriori item patterns much like AprioriGWAS (2014) paper.
            # Chr5_Pos26963862_A
            apriori_genotype_item_pattern = "_".join([marker_id, allele_base])
            
            # If the genotype_id is not present in the apriori_genotype_files_dict data structure then initialize an empty array and append the first item.
            if(not(str(genotype_id) in apriori_genotype_files_dict[str(part_num)])):
                apriori_genotype_files_dict[str(part_num)][str(genotype_id)] = []
                apriori_genotype_files_dict[str(part_num)][str(genotype_id)].append(apriori_genotype_item_pattern)
            # Else if the genotype_id is present in the apriori_genotype_files_dict data structure than append the item.
            elif(str(genotype_id) in apriori_genotype_files_dict[str(part_num)]):
                apriori_genotype_files_dict[str(part_num)][str(genotype_id)].append(apriori_genotype_item_pattern)
            
        # Row counter for figuring out what is header and what is a row entry.
        row_counter = row_counter + 1
    
    # Iterate over the genotype_ids list so that we can add the phenotype data to the encoded_individual_genotypes data structure.
    for genotype_id in genotype_ids:
        phenotype = phenotypes_dict[str(genotype_id)]
        if(str(genotype_id) in encoded_individual_genotypes):
            encoded_individual_genotypes[str(genotype_id)]['phenotype'] = phenotype
            
    # Iterate over the genotype_ids list and marker_ids list so that we can write the contents to a file.
    for genotype_id in sorted(genotype_ids):
        genotype_data_list = []
        for entry_key in sorted(encoded_individual_genotypes[str(genotype_id)]):
            genotype_data_list.append(encoded_individual_genotypes[str(genotype_id)][str(entry_key)])
        print(genotype_data_list)
        tsvwriter1.writerow([str(genotype_id)] + genotype_data_list)
       
    # Close the encoded genotype file
    tsvfile1.close()
        
    print(apriori_genotype_files_dict)
    for part_num in apriori_genotype_files_dict:
        print(str(part_num))
        
        apriori_genotypes_outfile = os.path.join(apriori_genotype_pattern_output_dir, "apriori_genotypes_part" + str(part_num) + ".tsv")

        # Writing the transaction genotypes database to a TSV file.
        tsvfile2 = open(apriori_genotypes_outfile, 'w')
        tsvwriter2 = csv.writer(tsvfile2, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

        # Write the header for the transaction genotypes database TSV file.
        tsvwriter2.writerow(["transaction_ids", "transaction_values"])
        for genotype_id in apriori_genotype_files_dict[str(part_num)]:
            transaction_values = ",".join(apriori_genotype_files_dict[str(part_num)][str(genotype_id)])
            tsvwriter2.writerow([str(genotype_id), transaction_values])
        #sys.exit()
        
        # Close the Apriori genotype transaction database.
        tsvfile2.close()


'''
annotations.gff.gz
'''
def parse_annotation_gff_file(gff_infile):

    # If the gff_infile is a compressed file then uncompress the file and return the file path of the uncompressed file.
    if(".gz" in gff_infile):
        gff_infile = gunzip_input_files(gff_infile, output_dir)

    # Counter for header and entry.
    row_counter = 0

    # Annotations data structure.
    annotations_dict = {}

    # Number of annoations counter.
    num_annotations = 0
    
    # Parse the annotation gff file.
    with open(gff_infile, 'r') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
        
        # Iterate over each row in the file.
        for row in csv_reader:
            if(not("#" in row)):
                print(row)
                sys.exit()
#           # If the line is not the header.
#            if(row_counter != 0):
#
#                #print(row)
#
#                # The genotype id of the phenotype.
#                genotype_id = row[0]
#
#                # The flowering time phenotype value.
#                flowering_time = row[1]
#
#                #print(genotype_id, flowering_time)
#
#                # Filter out phenotypes that have missing data in the form of "NA" and count the number of phenotypes.
#                if(flowering_time != "NA"):
#                    phenotypes_dict[genotype_id] = flowering_time
#                    num_annotations = num_annotations + 1
            row_counter = row_counter + 1
            
    # Close the annotations gff input file.
    csvfile.close()
    
    return(annotations_dict)


### MAIN FUNCTION ###

## Parse the genotypes input file using the phenotypes file and MAF threshold for filtering genotypes.
## Get the phenotypes and genotypes dictionary data structures for quick look up by genotype id.
#(genotypes_dict,phenotypes_dict) = parse_genotypes_file(genotypes_infile, phenotypes_infile, maf_threshold)
#
#
#(plink_genotype_ped_outfile, plink_genotype_map_outfile, mvp_phenotype_outfile) = generate_files_for_association_mapping(phenotypes_dict, genotypes_dict, output_dir)
#
## The association mapping output directory.
#association_mapping_output_dir = os.path.join(output_dir, "ASSOCIATION_MAPPING_OUTPUT_DIR")
#
## Create the emmax output directory if it does not exist.
#if not os.path.exists(association_mapping_output_dir):
#    os.makedirs(association_mapping_output_dir)
#
## Run the run_mvp_association_tests and obtain the adjusted_pvalues_infile.
#adjusted_pvalues_infile = run_mvp_association_tests(plink_genotype_ped_outfile, plink_genotype_map_outfile, mvp_phenotype_outfile, association_mapping_output_dir)
#
## The parsed_genotypes_output_dir output directory.
#parsed_genotypes_output_dir = os.path.join(output_dir, "PARSED_GENOTYPES_OUTPUT_DIR")
#
## Create the parsed genotypes output directory if it does not exist.
#if not os.path.exists(parsed_genotypes_output_dir):
#    os.makedirs(parsed_genotypes_output_dir)
#
## Parse the multiple corrected assocation test adjusted pvalues file to obtain the encoded genotypes file and and the apriori transaction genotypes database file.
#parse_multiple_corrected_tests(adjusted_pvalues_infile, genotypes_dict, phenotypes_dict, alpha_value, parsed_genotypes_output_dir)

parse_annotation_gff_file(gff_infile)
