# The library list of R package dependences.
library('getopt');
#library('qqman')
library('qvalue')
# The qqman_plot.r program usage example.

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'infile','i', 1, "character",
    'snp_hwe_results_outfile','o', 1, "character",
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
if(is.null(opt$snp_hwe_results_outfile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

## Initialize directory and file name variables.

# The input file parameter variable infile.
infile = opt$infile;

# The output directory parameter variable outfile_dir.
snp_hwe_results_outfile = opt$snp_hwe_results_outfile;

# generate base output directory if it does not exist.
#dir.create(output_dir, showWarnings = FALSE);


dataframe = read.table(infile, header=TRUE)

#head(dataframe)

#print(dataframe[,3])
qvalues = qvalue(p = dataframe[,3])

print(qvalues)

q_values = qvalues$qvalues

#print(qvalues$lfdr)

#summary(qvalues)
#print(qvalues$qvalues)
#head(qvalues)

# Get the number of markers in the genotype counts dataframe.
n_markers = nrow(dataframe)

# Write the header of the snp_hwe_results_outfile file.
cat(paste('marker_id', 'beta', 'p_value', 'q_value', sep='\t'), file=snp_hwe_results_outfile, sep='\n')
#cat(paste('MARKER_ID', 'HET', 'HOM_1', 'HOM_2', 'p_value', sep='\t'), file=snp_hwe_results_outfile, sep='\n')

# Iterate over the genotype counts dataframe.
for (i in 1:n_markers) {

    # Get the marker id at column 1.
    marker_id  = dataframe[i, 1]

    # Get the major allele frequency (p) at column 2.
    beta = dataframe[i, 2]

    # Get the minor allele frequency (q) at column 3.
    p_value = dataframe[i, 3]

    # Perform the exact SNP test of Hardy-Weinberg Equilibrium (HWE) to obtain a P-value for the marker.
    q_value = q_values[i]

    # Write the marker entry and results of the SNPHWE function to the snp_hwe_results_outfile file.
    cat(paste(marker_id,beta,p_value,q_value,sep='\t'),file=snp_hwe_results_outfile,sep='\n',append=TRUE)
    #cat(paste(marker_id,hets,hom_1,hom_2,p_value,sep='\t'),file=snp_hwe_results_outfile,sep='\n',append=TRUE)
}





