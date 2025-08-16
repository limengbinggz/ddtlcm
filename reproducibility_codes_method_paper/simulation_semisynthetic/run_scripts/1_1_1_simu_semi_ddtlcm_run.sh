#!/bin/bash
cd simulation_semisynthetic/run_scripts
seed_thetas=20 # number of class profiles to generate

write_slurm() {
    echo "#!/bin/bash
#SBATCH --job-name=dl_N$1_seedY$2_fix$3
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=3g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-${seed_thetas}
#SBATCH -o run_scripts/reports/%x_%A_%a.out # output records to /reports

cd simulation_semisynthetic/
module load R
Rscript --verbose 1_1_1_simu_semi_ddtlcm_run.R $1 \$SLURM_ARRAY_TASK_ID $2 $3" > simu_1_N$1_seedY$2_fix$3.slurm
if $4
then
    sbatch simu_1_N$1_seedY$2_fix$3.slurm
fi
}

run=true #false if only write scripts without submitting jobs
seed_Ys=5 # number of multivariate binary datasets to generate for each theta
fix_tree_ats=("FALSE" "truth" "misspecified")
Ns=(400 800)

for N in "${Ns[@]}"; do
    for (( seed_Y = 1; seed_Y <= seed_Ys; seed_Y++ )); do
        for fix_tree_at in "${fix_tree_ats[@]}"; do
            write_slurm ${N} ${seed_Y} ${fix_tree_at} ${run}
        done
    done
done