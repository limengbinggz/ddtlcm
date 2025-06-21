#!/bin/bash
cd simulation_semisynthetic/run_scripts
seed_thetas=20 # number of class profiles to generate

write_slurm() {
    echo "#!/bin/bash
#SBATCH --job-name=selectK_N$1_seedY$2_K$3_fold$4
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=3g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-${seed_thetas}
#SBATCH -o run_scripts/reports/%x_%A_%a.out # output records to /reports

cd simulation_semisynthetic/
module load R
Rscript --verbose 2_1_simu_semi_selectK_cv_run.R  $1 \$SLURM_ARRAY_TASK_ID $2 $3 $4" > sim_N$1_seedY$2_K$3_fold$4.slurm
if $5
then
    sbatch sim_N$1_seedY$2_K$3_fold$4.slurm
fi
}

# #true
run=true #false if only write scripts without submitting jobs
Ns=(400 800)
seed_Ys=5 # number of multivariate binary datasets to generate for each theta
K_candidates=(4 5 6 7 8)
fold_nums=5

for N in "${Ns[@]}"; do
    # for (( seed_theta = 1; seed_theta <= seed_thetas; seed_theta++ )); do
    for (( seed_Y = 1; seed_Y <= seed_Ys; seed_Y++ )); do
        for K_candidate in "${K_candidates[@]}"; do
            for (( fold_num = 1; fold_num <= fold_nums; fold_num++ )); do
                write_slurm ${N} ${seed_Y} ${K_candidate} ${fold_num} ${run}
            done
        done
    done
done