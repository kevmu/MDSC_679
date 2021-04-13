#!/bin/bash
#SBATCH --partition=cpu2019
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --mem=38G
#SBATCH --array=1-10%10
#SBATCH --output=apriori_test_dataset2_job_array.%A_%a.out
#SBATCH --error=apriori_test_dataset2_job_array.%A_%a.err

# The start time to print.
start_time=$(date)
echo "started at: ${start_time}"

# Make a logs directory for the logs files if it does not already exist.
#logs_dir="${HOME}/MDSC_679/ML_Project_2/test_datasets/test_dataset2/logs"
#mkdir -p $logs_dir

# The list of apriori_genotype_pattern_files 
list="${HOME}/MDSC_679/ML_Project_2/test_datasets/test_dataset2/test_dataset2_file_list.txt"

# The output directory.
output_dir="${HOME}/MDSC_679/ML_Project_2/test_datasets/test_dataset2/test_dataset2_output_dir"
mkdir -p $output_dir

IFS=$'\n' array=($(<$list))

# The apriori_genotype_pattern_file.
apriori_genotype_pattern_file=${array[$SLURM_ARRAY_TASK_ID-1]}

# Run the execute_apriori.py command using time to capture run time usage.
time python3 "${HOME}/MDSC_679/ML_Project_2/execute_apriori.py" --input_file $apriori_genotype_pattern_file --min_support_count 2 --min_confidence 0.60 --output_dir  $output_dir

# The end time to print.
end_time=$(date)
echo "finished with exit code $? at: ${end_time}"

