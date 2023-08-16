#!/usr/bin/bash 

# Script: rnaseq_gb2cov.sh
# Author: Ming Wang
# Version: 1.0
# Email: wangm08@hotmail.com
# Date: 2023-08-14
#
# Usage: rnaseq_gb2cov.sh <config.txt>

# Why this script?
# geneBody_coverage.py command from RSeQC takes much too long time. 2-7 housrs for 2-10 million reads
#
# How to:
# 1. run rnaseq pipeline, except genebody_cov, using multiple threads for each file
# 2. run this script alone, with single thread for each file


################################################################################
SRC_DIR=$(dirname $0) # script dir
GB_COV="${SRC_DIR}/genebody_cov.sh" # to-do
################################################################################

# load all global variables
# para: <config>
function load_config() {
    [[ -f $1 ]] || return 1
    source $1 # config
    # update resources
    WK_DIR=$(realpath -s ${WK_DIR})
    DATA_DIR=$(realpath -s ${DATA_DIR})
    [[ ${N_CPU} =~ ^[0-9]+$ ]] || N_CPU=8
    [[ ${N_JOB} =~ ^[0-9]+$ ]] || N_JOB=1
    [[ ${N_MEM} =~ ^[0-9]+$ ]] || N_MEM=30000 # MB, 30GB for Salmon,STAR alignment of hg38
    ## force
    N_CPU=1 #
    N_MEM=4000 #
}
export -f load_config


# check bam files in rnaseq project
# para: <bam_dir>
function check_bam_dir() {
    local bam_dir=$1
    local n_bam=$(ls ${bam_dir}/ | grep -c ".bam$") #
    if [[ ${n_bam} -ne 0 ]] 
    then
        echo "[${n_bam}] bam files found: ${bam_dir}" && return 0 # pass
    else
        echo "[0] no bam files found: ${bam_dir}" && return 1
    fi
}
export -f check_bam_dir


# run RNAseq_gb2cov to SLURM (HPC)
# para: <bam_list>
function make_slurm() {
    [[ $# -lt 1 ]] && echo "Usage: make_slurm <bam1> [bam2, ... bamN]" && return 1
    local bam_list="$@"
    local slurm_file="${WK_DIR}/run_rnaseq_gb2cov_slurm.sh"
    local out_dir="${WK_DIR}/results/06.genebody_cov"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # determine the total CPU/mem values
    local sum_CPU=$(echo ${N_CPU} ${N_JOB} | awk '{print $1*$2}')
    local sum_MEM=$(echo ${N_MEM} ${N_JOB} | awk '{print $1*$2}')

    cat << EOF > ${slurm_file}
#!/bin/bash
#SBATCH --job-name=gb2cov        # create a short name for your job
#SBATCH --partition=${SLURM_PARTITION}         # node name
#SBATCH --ntasks=${sum_CPU}               # total cores, 
#SBATCH --mem=${N_MEM}              # 30GB for STAR
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=5-01:00:00        # total run time limit (D-HH:MM:SS)
#SBATCH --output=log.${SLURM_PARTITION}.%j.out

module purge
# source /home/wangm/miniconda3/bin/activate
source "/home/wangm/miniconda3/etc/profile.d/conda.sh"
conda activate hiseq

cd ${WK_DIR}
bash ${GB_COV} ${out_dir} ${GENOME} ${N_JOB} ${bam_list}
EOF
    # output
    echo ${slurm_file}
}
export -f make_slurm 


# run RNAseq_gb2cov from terminal local, directly
# para: <bam_list>
function make_local() {
    [[ $# -lt 1 ]] && echo "Usage: make_local <bam1> [bam2, ... bamN]" && return 1
    local bam_list="$@"
    local cmd_file="${WK_DIR}/run_rnaseq_gb2cov_local.sh"
    local out_dir="${WK_DIR}/results/06.genebody_cov"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    echo "source $HOME/miniconda3/etc/profile.d/conda.sh" > ${cmd_file}
    echo "conda activate hiseq" >> ${cmd_file}
    echo "cd ${WK_DIR}" >> ${cmd_file}
    echo bash ${GB_COV} ${out_dir} ${GENOME} ${N_JOB} "${bam_list}" >> ${cmd_file}
    echo ${cmd_file}
}
export -f make_local


function rnaseq_gb2cov() {
    [[ $# -lt 2 ]] && echo "Usage: run_gb2cov <config> <run:0|1>" && return 1
    local config=$1
    local run=$2 # 0=not, 1=run
    load_config ${config} # global variables
    local bam_dir="${WK_DIR}/results/02.bam_files"
    check_bam_dir ${bam_dir} # global variables
    [[ $? -ne 0 ]] && echo "[error] - bam files not vaild: ${WK_DIR}/02.bam_files/" && return 1 # previous command failed
    local bam_list=$(ls ${bam_dir}/*bam | xargs)
    # prepare and run jobs
    [[ ${run} == 1 ]] && run_flag="yes" || run_flag="no" #
    if [[ -z $(sinfo -h -p ${SLURM_PARTITION}) ]] # global variables
    then
        # run on local server
        local cmd_file=$(make_local "${bam_list}")
        [[ $(grep -ci run_rnaseq.sh ${cmd_file}) -gt 0 ]] && job_flag="ok" || job_flag="failed"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${job_flag} "run" ${run_flag} $(basename ${cmd_file})
        [[ ${run} == 1 && ${job_flag} == "ok" ]] && bash ${cmd_file}
    else
        # run on HPC SLURM
        local slurm_file=$(make_slurm "${bam_list}")
        [[ $(grep -ci SBATCH ${slurm_file}) -gt 4 ]] && job_flag="ok" || job_flag="failed"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${job_flag} "run" ${run_flag} $(basename ${slurm_file})
        [[ ${run} == 1 && ${job_flag} == "ok" ]] && sbatch ${slurm_file}
    fi
}
export -f rnaseq_gb2cov


rnaseq_gb2cov $@
