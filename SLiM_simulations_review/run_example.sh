#!/bin/sh

#SBATCH --ntasks=1 
#SBATCH --mem=32000
#SBATCH --cpus-per-task=1
#SBATCH -p hns,normal
#SBATCH --time=05:00:00
#SBATCH --job-name=slim

module load gsl
module load python/3.6.1
module load R

#usage:
#perl introgression_simulations_v7.pl config_file_name out_folder start_sim_num numreps admixture_prop numDMIs s-DMI s-ancestral-geno number-load s-load Ne

#example - load
perl introgression_simulations_v7.pl introgression_DMI_load_config_v4.slim slim_out 1 5 0.25 0 0.0 0.0 150 -0.001 10000

#example - DMI
perl introgression_simulations_v7.pl introgression_DMI_load_config_v4.slim slim_out_DMI 1 5 0.25 3 -0.3 0.0 0 0.0 10000
