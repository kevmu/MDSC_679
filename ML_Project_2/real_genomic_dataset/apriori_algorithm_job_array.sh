#!/bin/bash
#SBATCH --partition=cpu2019
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem=38G
#SBATCH --array=1-4%4
#SBATCH --output=apriori_algorithm_job_array.%A_%a.out
#SBATCH --error=apriori_algorithm_job_array.%A_%a.err

# The list of apriori_genotype_pattern_files 
list="${HOME}/MDSC_679/ML_Project_2/real_genomic_dataset/apriori_genotype_pattern_file_list.txt"

output_dir="${HOME}/MDSC_679/ML_Project_2/real_genomic_dataset"
mkdir -p $output_dir


IFS=$'\n' array=($(<$list))

# The apriori_genotype_pattern_file.
apriori_genotype_pattern_file=${array[$SLURM_ARRAY_TASK_ID-1]}

python3 "${HOME}/MDSC_679/ML_Project_2/execute_apriori.py" -i $apriori_genotype_pattern_file -o $output_dir


