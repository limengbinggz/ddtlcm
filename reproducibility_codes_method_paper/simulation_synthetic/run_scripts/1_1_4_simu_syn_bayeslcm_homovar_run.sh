#!/bin/bash
cd simulation_synthetic/run_scripts
seed_thetas=20 # number of class profiles to generate

write_slurm() {
    echo "#!/bin/bash
#SBATCH --job-name=indhomo_tree$1_N$2_seedY$3
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=3g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-${seed_thetas}
#SBATCH -o run_scripts/reports/%x_%A_%a.out # output records to /reports

cd simulation_synthetic
module load R
Rscript --verbose 1_1_4_simu_syn_bayeslcm_homovar_run.R $1 $2 \$SLURM_ARRAY_TASK_ID $3" > simu_1_tree$1_N$2_seedY$3.slurm
if $4
then
    sbatch simu_1_tree$1_N$2_seedY$3.slurm
fi
}

run=true #false if only write scripts without submitting jobs
tree_idxs=(1 2 3 4)
Ns=(100 200 400)
seed_Ys=5 # number of multivariate binary datasets to generate for each theta

for tree_idx in "${tree_idxs[@]}"; do
    for N in "${Ns[@]}"; do
        for (( seed_Y = 1; seed_Y <= seed_Ys; seed_Y++ )); do
            write_slurm ${tree_idx} ${N} ${seed_Y} ${run}
        done
    done
done