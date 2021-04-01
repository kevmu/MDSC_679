# The library list of R package dependences.
library('getopt');
#library('qqman')
library('gaston')
# The qqman_plot.r program usage example.

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'infile','i', 1, "character",
    'output_dir','o', 1, "character",
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
if(is.null(opt$output_dir)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

## Initialize directory and file name variables.

# The input file parameter variable infile.
infile = opt$infile;

# The output directory parameter variable outfile_dir.
output_dir = opt$output_dir;

# generate base output directory if it does not exist.
dir.create(output_dir, showWarnings = FALSE);


dataframe = read.table(infile, header=TRUE)

summary(dataframe)
png("question3_plot.png")
print(sort(as.vector(dataframe$P)))
#qq(sort(as.vector(dataframe$P)))
#manhattan(dataframe)

qqplot.pvalues(dataframe$P, col.abline = "red", CB = TRUE, col.CB = "gray80", CB.level = 0.95)
dev.off()

