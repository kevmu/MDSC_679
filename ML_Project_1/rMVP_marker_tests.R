# Name: Kevin Muirhead
# UCID#: 00502756

# rMVP_marker_tests.R - Performs the rMVP association tests for SNP variant markers.


# The library list of R package dependences.
library('getopt');
library('rMVP')

# The R program usage example.
# Rscript rMVP_marker_tests.R -i genotypes.vcf -p phenotype_infile -o association_mapping_output_dir

# Get options, using the spec as defined by the enclosed list.
# We read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'vcf_infile','i', 1, "character",
    'phenotype_infile','p', 1, "character",
    'association_mapping_output_dir','o', 1, "character",
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
if(is.null(opt$vcf_infile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}
if(is.null(opt$phenotype_infile)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}
if(is.null(opt$association_mapping_output_dir)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

## Initialize directory and file name variables.

# The VCF input file parameter variable vcf_infile.
vcf_infile = opt$vcf_infile;

# The phenotype input file parameter variable phenotype_infile.
phenotype_infile = opt$phenotype_infile;

# The output directory parameter variable association_mapping_output_dir.
association_mapping_output_dir = opt$association_mapping_output_dir;

out_prefix=paste(association_mapping_output_dir,"mvp.genotype.vcf", sep="/")


# Generate the mvp data files for input into the MVP function.
MVP.Data(

    fileVCF=vcf_infile,
    filePhe=phenotype_infile,

    sep.map="\t",
    sep.phe="\t",
    fileKin=FALSE,
    filePC=FALSE,
    out=out_prefix
)

#genotype <- attach.big.matrix("/Users/kevin.muirhead/Desktop/GWAS_output_dir1/FILES_FOR_ASSOCIATION_MAPPING/mvp.genotype.vcf.geno.desc")
#phenotype <- read.table("/Users/kevin.muirhead/Desktop/GWAS_output_dir1/FILES_FOR_ASSOCIATION_MAPPING/mvp.genotype.vcf.phe",head=TRUE)
#map <- read.table("/Users/kevin.muirhead/Desktop/GWAS_output_dir1/FILES_FOR_ASSOCIATION_MAPPING/mvp.genotype.vcf.geno.map" , head = TRUE)

genotype <- attach.big.matrix(paste(association_mapping_output_dir,"mvp.genotype.vcf.geno.desc", sep="/"))
phenotype <- read.table(paste(association_mapping_output_dir,"mvp.genotype.vcf.phe", sep="/"),head=TRUE)
map <- read.table(paste(association_mapping_output_dir,"mvp.genotype.vcf.geno.map", sep="/"), head = TRUE)

imMVP <- MVP(
    phe=phenotype,
    geno=genotype,
    map=map,
#    nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
#    nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
#    nPC.FarmCPU=3,
#    priority="speed",   ##for Kinship construction
    ncpus=2,
#    vc.method="BRENT",  ##only works for MLM
    maxLoop=100,
    method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
    #permutation.threshold=TRUE,
    #permutation.rep=100,
#    p.threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU"),
#    method="FarmCPU",

#    method="MLM",

    file.output=TRUE,

#    outpath="/Users/kevin.muirhead/Desktop/GWAS_output_dir1/ASSOCIATION_MAPPING_OUTPUT_DIR/"
    outpath=association_mapping_output_dir

)


