# MDSC_679: Machine Learning Projects

Name: Kevin Muirhead
<br>
UCID #: 00502756
<br>
<br>
Contains the source code, documentation and Microsoft Word Documents for the MDSC 679 Course ML_Project_1 and ML_Project_2.

## Installation

To use this pipeline, clone this repository into your project directory using the following command:

```
project_dir="$HOME/software"
mkdir -p $project_dir
cd $project_dir
git clone https://github.com/kevmu/MDSC_679.git
```
Install the miniconda3 conda environment.

```
install_dir=$HOME
cd $install_dir

wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
sh $miniconda_install_file -b -p "$install_dir/Miniconda3"
$install_dir/Miniconda3/bin/conda init
export PATH="${install_dir}/Miniconda3/bin:$PATH"

```

## MDSC_679: ML_Project_1

Install the ML_Project_1 conda environment.

```
# Create the ML_Project_1_env conda environment
conda create --name ML_Project_1_env

# Activate the conda environment.
conda activate ML_Project_1_env

# Install the plink program conda package.
conda install -c bioconda plink

# Install the r-base package for R.
conda install -c conda-forge r-base

# Open the R console terminal.
R

## In R terminal

# Install rMVP.
install.packages("rMVP", dependencies=TRUE)

# Test library is installed.
library("rMVP")

# Exit the R console terminal.
q()

# Make the ML_Model_env conda environment
conda create --name ML_Model_env

# Activate the ML_Model_env conda environment
conda activate ML_Model_env

# Install the numpy package for python. Manipulate numpy arrays.
conda install -c conda-forge numpy

# Install the scikit-learn module for python. For machine learning libraries.
conda install -c anaconda scikit-learn

# Install the pandas module for python. Uploading table files and Manipulating dataframes.
conda install -c anaconda pandas

# Install the matplotlib module for python. Generate publication quality plots.
conda install -c conda-forge matplotlib

```
Make sure that you are in the ML_Project_1 directory to use the rMVP_marker_tests.R and calculate_adjusted_pvalues.R rscript files.

Program Usage quality_control.py:

```
# Activate the ML_Project_1_env conda environment.
conda activate ML_Project_1_env

# Execute the quality_control.py python script for filtering genotypes, perform association mapping using the rMVP R package script, quality control and formatting files for input into the machine learning models; the polgenetic linear model, LASSO L1 Regression, and RIDGE L2 Regression.
python $HOME/software/MDSC_679/ML_Project_1/quality_control.py --phenotypes_infile $HOME/software/MDSC_679/ML_Project_1/INPUT_FILES/FT10.txt --genotypes_infile $HOME/software/MDSC_679/ML_Project_1/INPUT_FILES/genotype.csv.gz --gff_infile $HOME/software/MDSC_679/ML_Project_1/INPUT_FILES/gene_model.gff.gz --alpha_value 0.05 --maf_threshold 0.01 --output_dir $HOME/GWAS_OUTPUT_DIR

```

Program Input Parameters:

| Parameters | Description |
| ---------- | ----------- |
| --phenotypes_infile | The file path containing list of sample names or ids (string) |
| --genotypes_infile | The directory path for input fastq files (string) |
| --gff_infile | The gene model annotation GFF format file (string) |
| --maf_threshold | The minor allele frequency (MAF) threshold for filtering genotypes. (i.e. Default: 0.01) (float) |
| --alpha_value | The alpha value as input for filtering adjusted pvalues of the rMVP association tests by using adjusted pvalues using Bonferroni correction (i.e. Default: 0.05) |
| --output_dir | The output directory to write the output directories and files. (string) |

```

Program Usage polygenetic_linear_regression.py:

```
# Activate the ML_Model_env conda environment.
conda activate ML_Model_env

# Execute the polygenetic_linear_regression.py python script for performing Polygenetic Linear Regression.
python $HOME/software/MDSC_679/ML_Project_1/polygenetic_linear_regression.py

```

Program Input Parameters:

| Parameters | Description |
| ---------- | ----------- |
| --phenotypes_infile | The file path containing list of sample names or ids (string) |
| --output_dir | The output directory to write the output directories and files. (string) |

```

Program Usage lasso_regression.py:

```
# Activate the ML_Model_env conda environment.
conda activate ML_Model_env

# Execute the lasso_regression.py python script for performing LASSO Regression (L1 Regularization).
python $HOME/software/MDSC_679/ML_Project_1/lasso_regression.py


```

Program Input Parameters:

| Parameters | Description |
| ---------- | ----------- |
| --phenotypes_infile | The file path containing list of sample names or ids (string) |
| --output_dir | The output directory to write the output directories and files. (string) |


```

Program Usage ridge_regression.py:

```
# Activate the ML_Model_env conda environment.
conda activate ML_Model_env

# Execute the ridge_regression.py python script for performing RIDGE Regression (L2 Regularization).
python $HOME/software/MDSC_679/ML_Project_1/ridge_regression.py

```

Program Input Parameters:

| Parameters | Description |
| ---------- | ----------- |
| --phenotypes_infile | The file path containing list of sample names or ids (string) |
| --output_dir | The output directory to write the output directories and files. (string) |


## MDSC_679: ML_Project_2

The Implementation of the AprioriTID algorithm from the following paper;
<br>
<br>
“Fast Algorithms for Mining Association Rules”, Agrawal, R., Ramakrishnan, S. (1994), Proc. 20th int. conf. very large data bases, VLDB 1215(pp. 487-499). doi: 10.1.1.40.7506
<br>
<br>
Contains the Design and Implementation Report Microsoft Word Document, source code and source code documentation.

Program Usage execute_apriori.py:

```

# Execute the execute_apriori.py python script for performing the Apriori Algorithm, calculates association rule metrics and prints the association rule metrics support, confidence, lift, leverage and conviction.

python $HOME/software/MDSC_679/ML_Project_2/execute_apriori.py --input_file $HOME/software/MDSC_679/ML_Project_2/test_datasets/test_dataset1/transaction_database1.tsv --min_support_count 2 --min_confidence 0.60 --output_dir  $HOME/software/MDSC_679/ML_Project_2/test_dataset1_output_dir

Make sure that you are in the ML_Project_2 directory to use the Apriori.py python class file.

Program Input Parameters:

| Parameters | Description |
| ---------- | ----------- |
| --database_infile | The transaction database file as input. (i.e. $HOME/filename.tsv) (string) |
| --min_support_count | The minimum support count as input. Can be any integer value greater than zero. min_support_count > 0. Default: 2 (integer) |
| --min_confidence | The minimum confidence as input. Can be in the range of 0.00-1.00. Default: 0.60 (float) |
| --output_dir | The output directory as input. (i.e. $HOME/output_dir) (string) |



