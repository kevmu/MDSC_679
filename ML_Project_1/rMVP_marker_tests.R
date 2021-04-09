library(rMVP)


MVP.Data(
    #fileNum="/Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/demo_data/numeric.txt",

    #filePhe="/Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/demo_data/mdp_traits_validation.txt",
    #fileMap="/Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/demo_data/mdp_SNP_information.txt",
fileVCF"/Users/kevin.muirhead/Desktop/GWAS_output_dir1/plink_format_files/test.vcf",
fileMap="/Users/kevin.muirhead/Desktop/GWAS_output_dir1/plink_format_files/genotypes.map",
filePhe="/Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/demo_data/mdp_traits_validation.txt",

             sep.map="\t",
             sep.phe="\t",
             fileKin=FALSE,
             filePC=FALSE,
            auto_transpose=TRUE,
             #priority="memory",
             #maxLine=10000,
             out="mvp.num"
)

genotype <- attach.big.matrix("mvp.num.geno.desc")
phenotype <- read.table("mvp.num.phe",head=TRUE)
map <- read.table("mvp.num.geno.map" , head = TRUE)

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
#    method=c("GLM", "MLM", "FarmCPU"),
    method="FarmCPU",

#    method="MLM",

    file.output=TRUE,

    outpath="/Users/kevin.muirhead/Desktop/MDSC_679/ML_Project_1/demo_data/"

)


