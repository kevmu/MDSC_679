

## Install EMMAX

# Download the linux emmax binary. 
wget http://csg.sph.umich.edu//kang/emmax/download/emmax-beta-07Mar2010.tar.gz

# Untar the emmax tarball.
tar xvzf emmax-beta-07Mar2010.tar.gz 

# Change directory to start testing to see if the command works on your system before proceeding.
cd emmax-beta-07Mar2010/


# Install the R conda environment for Linux.
conda env create --file r_env.yaml

## Install R library packages in R.

# Execute the R console command. 
R

# The install packages list of commands to install R package dependencies.
# Use to install R get options parameters "getopt" library.
install.packages("getopt", dependencies=TRUE)

# Use to install R QQ-plot and manhattan plot "qqman" library.
install.packages("qqman", dependencies=TRUE)

# Use to install R jdstorey/qvalue "qvalue" library.
install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")

