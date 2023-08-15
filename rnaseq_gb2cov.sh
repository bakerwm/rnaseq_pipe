

################################################################################
# GLOBAL VARIABLES
SLURM_PARTITION="CPU1"
N_CPU=1  # --ntasks
MEM=4000 # --mem MB
################################################################################

## make slurm script
function make_slurm() {
    local wk_dir=$1
    local bam=$2
    [[ $# -lt 2 ]] && echo "not enough args..." && return 1
    wk_dir=$(realpath -s ${wk_dir})
    bam=$(realpath -s ${bam})
    ## check run_pipe.sh
    gb2cov="${HOME}/biosoft/run_rnaseq/genebody_cov.sh"
    [[ ! -f ${gb2cov} ]] && echo "genebody_cov.sh not exists" && return 1

    cat << EOF
#!/bin/bash
#SBATCH --job-name=gb2cov        # create a short name for your job
#SBATCH --partition=${SLURM_PARTITION}         # node name
#SBATCH --ntasks=${N_CPU}               # total cores, 
#SBATCH --mem=${MEM}              # 30GB for STAR
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=5-01:00:00        # total run time limit (D-HH:MM:SS)
#SBATCH --output=log.${SLURM_PARTITION}.%j.out

module purge
# source /home/wangm/miniconda3/bin/activate
source "/home/wangm/miniconda3/etc/profile.d/conda.sh"
conda activate hiseq

cd ${wk_dir}
# bash ${gb2cov} ${wk_dir}/results/genebody_cov hg38 1 ${bam}
bash ${gb2cov} ${wk_dir}/results/06.genebody_cov hg38 1 ${bam}

EOF
}
export -f make_slurm 


function make_slurm_sh() {
    local wk_dir=$1
    local bam=$2
    [[ $# -lt 2 ]] && echo "missing args" && return 1
    local fname=$(basename ${bam%.bam})
    local job_dir=${wk_dir}/jobs
    [[ ! -d ${job_dir} ]] && mkdir -p ${job_dir}
    local job_file="${job_dir}/submit_${fname}_slurm.sh"
    # sbatch run_slurm.sh ${fq1}
    make_slurm ${wk_dir} ${bam} > ${job_file}
    echo ${job_file}    
}
export -f make_slurm_sh


function run_gb2cov() {
    [[ $# -lt 3 ]] && echo "Usage: run_gb2cov.sh <wk_dir> <bam_dir> <0|1>" && return 1
    local wk_dir=$1
    local bam_dir=$2
    local run=$3 # 0|1
    ##
    for bam in ${bam_dir}/*.bam
    do
        [[ ! -f ${bam} ]] && continue
        ## create job ##
        job_file=$(make_slurm_sh ${wk_dir} ${bam})
        ## status ##
        [[ $(grep -ci SBATCH ${job_file}) -gt 4 ]] && flag="ok" || flag="failed"
        [[ ${run} == 1 ]] && ss="yes" || ss="no"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${flag} "submit" ${ss} $(basename ${job_file})
        ## submit ##
        [[ ${run} == 1 && ${flag} == "ok" ]] && sbatch ${job_file}
    done
}
export -f run_gb2cov

# [[ $# -lt 3 ]] && echo "Usage: run_gb2cov.sh <wk_dir> <bam_dir> <0|1>" && exit 1
run_gb2cov $@
