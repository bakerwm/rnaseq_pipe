
# run rnaseq pipeline
# Usage
# rnaseq.sh config.txt


SRC_DIR=$(dirname $(realpath -s $0)) # script dir
# echo "SRC_DIR=${SRC_DIR}"

# load config
function load_config() {
    conf=$1
    if [[ -f ${conf} ]] 
    then 
        while read -r line
        do
            eval "$line"
        done < "$conf"
    else
        echo "config file not exists: ${conf}"
    fi
    # required global variables:
    [[ -z ${GENOME} ]] && GENOME="hg38"
    [[ -z ${RAW_DIR} ]] && RAW_DIR="./"
    [[ -z ${WK_DIR} ]] && WK_DIR="./"
    [[ -z ${N_CPU} ]] && N_CPU=8
    [[ -z ${N_JOB} ]] && N_JOB=1
}
export -f load_config


## make slurm script
function prep_slurm() {
    local wk_dir=$1
    local fq1=$2
    local fq2=$3
    [[ $# -lt 3 ]] && echo "not enough args..." && return 1
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    wk_dir=$(realpath -s ${wk_dir})
    fq1=$(realpath -s ${fq1})
    fq2=$(realpath -s ${fq2})
    ## global variables
    SRC_DIR=$(dirname $(realpath -s $0)) # script dir
    [[ -z ${SLURM_PARTITION} ]] && SLURM_PARTITION="CPU2"
    [[ -z ${N_CPU} ]] && N_CPU=8
    # ## prepare job files    
    # local fname=$(basename ${fq1%_1.fq.gz})
    # local job_dir=${wk_dir}/jobs
    # [[ ! -d ${job_dir} ]] && mkdir -p ${job_dir}
    # local job_file="${job_dir}/submit_${fname}_slurm.sh"
    ## check run_pipe.sh
    run_rnaseq="${SRC_DIR}/run_rnaseq.sh"
    [[ ! -f ${run_rnaseq} ]] && echo "run_rnaseq.sh not exists" && return 1

    cat << EOF
#!/bin/bash
#SBATCH --job-name=rna           # create a short name for your job
#SBATCH --partition=${SLURM_PARTITION}         # node name
#SBATCH --ntasks=${N_CPU}               # total cores, 
#SBATCH --mem=30000              # 30GB for STAR
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=5-01:00:00        # total run time limit (D-HH:MM:SS)
#SBATCH --output=log.${SLURM_PARTITION}.%j.out

module purge
# source ${HOME}/miniconda3/bin/activate
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate hiseq

cd ${wk_dir}
bash ${run_rnaseq} ${wk_dir} ${fq1} ${fq2}
EOF
}
export -f prep_slurm 


# script_dir
function make_slurm() {
    [[ $# -lt 3 ]] && echo "missing args" && return 1
    local wk_dir=$1
    local fq1=$2
    local fq2=$3
    [[ ! ${fq1} = *1.fq.gz ]] && echo "not fq1: ${fq1}" && return 1
    [[ ! ${fq2} = *2.fq.gz ]] && echo "not fq2: ${fq1}" && return 1
    # local fq2=${fq1/1.fq/2.fq}
    local fname=$(basename ${fq1%_1.fq.gz})
    local job_dir=${wk_dir}/jobs
    [[ ! -d ${job_dir} ]] && mkdir -p ${job_dir}
    local job_file="${job_dir}/submit_${fname}_slurm.sh"
    # sbatch run_slurm.sh ${fq1}
    prep_slurm ${wk_dir} ${fq1} ${fq2} > ${job_file}
    echo ${job_file}    
}
export -f make_slurm


## run RNAseq directly, on server
function rnaseq_direct() {
    [[ $# -lt 2 ]] && echo "Usage: rnaseq.sh <config.txt> [1|0]" && exit 1
    local config=$1
    local run=$2 # 0=not, 1=run
    load_config $1 # init arguments, 
    ## required: GENOME, RAW_DIR, WK_DIR, N_CPU
    ## show config
    echo "--------------------------------------------------------------------------------"
    echo "GENOME:          ${GENOME}"
    echo "N_CPU:           ${N_CPU}"
    echo "N_JOB:           ${N_JOB}"
    echo "slurm_partition: ${SLURM_PARTITION}"
    echo "RAW_DIR:         ${RAW_DIR}"
    echo "WK_DIR:          ${WK_DIR}"
    echo "job_dir:         ${WK_DIR}/jobs"
    echo "--------------------------------------------------------------------------------"
    echo ""

    # source ${HOME}/miniconda3/bin/activate
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate hiseq

    # prepare slurm jobs
    # replace string
    [[ ${run} == 1 ]] && cmd="bash" || cmd="echo"
    parallel -j ${N_JOB} --rpl '{/(.+?)/(.*?)} s/$$1/$$2/;' ${cmd} "${SRC_DIR}/run_rnaseq.sh" ${WK_DIR} {} {/_1.fq/_2.fq} ::: ${RAW_DIR}/*1.fq.gz
}
export -f rnaseq_direct


# check if partition exists
function is_valid_slurm_partition() {
    pt=$(sinfo -h -p $1)
    [[ -z ${pt} ]] && echo 1 || echo 0
}
export -f is_valid_slurm_partition


## run RNAseq on SLURM
function rnaseq_slurm() {
    [[ $# -lt 2 ]] && echo "Usage: rnaseq.sh <config.txt> [1|0]" && exit 1
    local config=$1
    local run=$2 # 0=not, 1=run
    load_config $1 # init arguments, 
    ## required: GENOME, RAW_DIR, WK_DIR, N_CPU
    ## show config
    echo "--------------------------------------------------------------------------------"
    echo "GENOME:          ${GENOME}"
    echo "N_CPU:           ${N_CPU}"
    echo "slurm_partition: ${SLURM_PARTITION}"
    echo "RAW_DIR:         ${RAW_DIR}"
    echo "WK_DIR:          ${WK_DIR}"
    echo "job_dir:         ${WK_DIR}/jobs"
    echo "--------------------------------------------------------------------------------"
    echo ""

    # prepare slurm jobs
    for fq1 in ${RAW_DIR}/*1.fq.gz
    do
        local fq2=${fq1/1.fq/2.fq}
        [[ -f ${fq1} && -f ${fq2} ]] || continue
        ## create job ##
        job_file=$(make_slurm ${WK_DIR} ${fq1} ${fq2})
        ## status ##
        [[ $(grep -ci SBATCH ${job_file}) -gt 4 ]] && job_flag="ok" || job_flag="failed"
        [[ ${run} == 1 ]] && run_flag="yes" || run_flag="no"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${job_flag} "submit" ${run_flag} $(basename ${job_file})
        ## submit ##
        [[ ${run} == 1 && ${job_flag} == "ok" ]] && sbatch ${job_file}
    done
}
export -f rnaseq_slurm


function main() {
    [[ $# -lt 2 ]] && echo "Usage: rnaseq.sh <config.txt> [1|0]" && exit 1
    local config=$1
    # local run=$2 # 0=not, 1=run
    load_config $1 # init arguments

    ## check if SLURM or not
    if [[ -z $(sinfo -h -p ${SLURM_PARTITION}) ]] 
    then 
        rnaseq_direct $@
    else 
        rnaseq_slurm $@
    fi
}
export -f main

main $@

