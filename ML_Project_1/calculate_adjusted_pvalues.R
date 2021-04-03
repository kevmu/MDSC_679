# The library list of R package dependences.
library('getopt');
library('qvalue')

# The R program usage example.
# Rscript adjust_pvalues.R -i emmax.ps -o emmax_adjusted_pvalues.txt

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'infile','i', 1, "character",
    'adjusted_pvalue_results_outfile','o', 1, "character",
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
if(is.null(opt$adjusted_pvalue_results_outfile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

## Initialize directory and file name variables.

# The input file parameter variable infile.
infile = opt$infile;

# The output directory parameter variable outfile_dir.
adjusted_pvalue_results_outfile = opt$adjusted_pvalue_results_outfile;

dataframe = read.table(infile, header=TRUE)

#head(dataframe)

#print(dataframe[,3])

# Get the pvalues for calculating adjusted pvalues
pvalues = dataframe[,3]

# Calculate the Bonferroni correction to adjust pvalues.
bonf_corr_pvalues = p.adjust(p = pvalues, method = "bonferroni")

# Calculate the FDR qvalue adjusted pvalues.
qvalues = qvalue(p = pvalues)

#print(qvalues)

q_values = qvalues$qvalues

#print(qvalues$lfdr)

#summary(qvalues)
#print(qvalues$qvalues)
#head(qvalues)

# Get the number of markers in the genotype counts dataframe.
n_markers = nrow(dataframe)

# Write the header of the adjusted_pvalue_results_outfile file.
cat(paste('marker_id', 'beta', 'p_value', 'bonf_corr_pvalue', 'q_value', sep='\t'), file=adjusted_pvalue_results_outfile, sep='\n')

# Iterate over the genotype counts dataframe.
for (i in 1:n_markers) {

    # Get the marker id at column 1.
    marker_id  = dataframe[i, 1]

    # Get the emmax algorithm beta value for the marker at column 2.
    beta = dataframe[i, 2]

    # Get the emmax algorithm pvalue for the marker at column 3.
    p_value = dataframe[i, 3]

    # Get the Bonferroni corrected pvalue for the marker.
    bonf_corr_pvalue = bonf_corr_pvalues[i]

    # Get the FDR adjusted pvalue (qvalue) for the marker.
    q_value = q_values[i]

    # Write the marker entry and results of the adjusted pvalues to the adjusted_pvalue_results_outfile file.
    cat(paste(marker_id,beta,p_value,bonf_corr_pvalue,q_value,sep='\t'),file=adjusted_pvalue_results_outfile,sep='\n',append=TRUE)
}





