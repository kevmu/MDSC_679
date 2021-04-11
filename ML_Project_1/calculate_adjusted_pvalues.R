# The library list of R package dependences.
library('getopt');
#library('qvalue')

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

# Read the input file into a dataframe.
dataframe = read.csv(infile, header=TRUE)

# Get the pvalues for calculating adjusted pvalues from column 8.
pvalues = dataframe[, 8]

# Calculate the Bonferroni correction to adjust pvalues.
bonf_corr_pvalues = p.adjust(p = pvalues, method = "bonferroni")

# Calculate the FDR to adjust pvalues using the BY method. Benjamini & Yekutieli (2001) ("BY").
qvalues = p.adjust(p = pvalues, method = "BY")

#print(qvalues)

# Get the number of markers in the genotype counts dataframe.
num_markers = nrow(dataframe)

# Write the header of the adjusted_pvalue_results_outfile file.
cat(paste('SNP', 'CHROM', 'POS', 'REF', 'ALT', 'Effect', 'SE', 'p_value', 'bonf_corr_pvalue', 'q_value', sep='\t'), file=adjusted_pvalue_results_outfile, sep='\n')

# Iterate over the genotype counts dataframe.
for (i in 1:num_markers) {

    # SNP    CHROM    POS    REF    ALT    Effect    SE    phenotype.MLM
    # Get the SNP id at column 1.
    snp_id  = dataframe[i, 1]
    
    # Get the chromosome id at column 2.
    chromosome_id  = dataframe[i, 2]

    # Get the position id at column 3.
    position_id = dataframe[i, 3]

    # Get the reference allele at column 4.
    ref_allele = dataframe[i, 4]

    # Get the alternative allele at column 5.
    alt_allele = dataframe[i, 5]
    
    # Get the mvp MLM effect at column 6.
    effect = dataframe[i, 6]
    
    # Get the mvp MLM se at column 7.
    se = dataframe[i, 7]
    
    # Get the mvp MLM se at column 8.
    p_value = dataframe[i, 8]
    
    # Get the Bonferroni corrected pvalue for the marker.
    bonf_corr_pvalue = bonf_corr_pvalues[i]

    # Get the FDR adjusted pvalue (qvalue) for the marker.
    q_value = qvalues[i]

    # Write the marker entry and results of the adjusted pvalues to the adjusted_pvalue_results_outfile file.
    cat(paste(snp_id,chromosome_id,position_id,ref_allele,alt_allele,effect,se,p_value,bonf_corr_pvalue,q_value,sep='\t'),file=adjusted_pvalue_results_outfile,sep='\n',append=TRUE)
}





