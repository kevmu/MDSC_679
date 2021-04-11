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

## MDSC_679: ML_Project_2

The Implementation of the AprioriTID algorithm from the following paper;
<br>
<br>
“An Improved Apriori Algorithm For Association Rules.”, Al-Maolegi, Mohammed & Arkok, Bassam. (2014). International Journal on Natural Language Computing. 3. 10.5121/ijnlc.2014.3103.
<br>
<br>
Contains the Design and Implementation Report Microsoft Word Document, source code and source code documentation.
