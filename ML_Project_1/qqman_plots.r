# The library list of R package dependences.
library('getopt');
library(qqman)

# The qqman_plot.r program usage example.

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

# The output file parameter variable outfile.
outfile = opt$outfile;

dataframe = read.table(infile, header=TRUE)

