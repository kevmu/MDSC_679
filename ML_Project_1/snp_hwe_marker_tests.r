# The library list of R package dependences.
library('getopt');

# Adapted from "R version of SNP-HWE" (http://csg.sph.umich.edu/abecasis/Exact/r_instruct.html).
source("snp_hwe.r")

# The snp_hwe_marker_tests.r program usage example.
# rscript snp_hwe_marker_tests.r -i ~/Desktop/macbook_air/MDSC_679/ML_Project_1/genotypes_counts.tsv -o ~/Desktop/macbook_air/MDSC_679/ML_Project_1/SNPHWE_Results/genotype_count_SNP-HWE_tests.tsv

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'infile','i', 1, "character",
    'outfile','o', 1, "character",
    'help','h', 0, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if(!is.null(opt$help)){
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

# if no parameter was given for any of these print
# a friendly message and exit with a non-zero error
# code
if(is.null(opt$infile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}
if(is.null(opt$outfile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

## Initialize directory and file name variables.

# The input file parameter variable infile.
infile = opt$infile;

# The output file parameter variable snp_hwe_results_outfile.
snp_hwe_results_outfile = opt$outfile;

# Read in the genotype counts file as a dataframe containing rows with a marker_id, heterozygote count, homozygote 1 count and homozygote 2 count.
dataframe = read.table(infile, header=TRUE)

# Get the number of markers in the genotype counts dataframe.
n_markers = nrow(dataframe)

# Write the header of the snp_hwe_results_outfile file.
cat(paste('MARKER_ID', 'p', 'q', 'HET', 'HOM_1', 'HOM_2', 'p_value', sep='\t'), file=snp_hwe_results_outfile, sep='\n')

# Iterate over the genotype counts dataframe.
for (i in 1:n_markers) {
   
    # Get the marker id at column 1.
    marker_id  = dataframe[i, 1]
    
    # Get the major allele frequency (p) at column 2.
    p = dataframe[i, 2]
    
    # Get the minor allele frequency (q) at column 3.
    q = dataframe[i, 3]
    
    # Get the heterozygote count at column 4.
    hets = dataframe[i, 4]
    
    # Get the homozygote 1 count at column 5.
    hom_1 = dataframe[i, 5]
    
    # Get the homozygote 2 count at column 6.
    hom_2 = dataframe[i, 6]

    # Perform the exact SNP test of Hardy-Weinberg Equilibrium (HWE) to obtain a P-value for the marker.
    p_value = SNPHWE(hets, hom_1, hom_2)
    
    # Write the marker entry and results of the SNPHWE function to the snp_hwe_results_outfile file.
    cat(paste(marker_id,hets,hom_1,hom_2,p_value,sep='\t'),file=snp_hwe_results_outfile,sep='\n',append=TRUE)

}




