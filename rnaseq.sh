#!/usr/bin/bash 

# Script: rnaseq.sh
# Author: Ming Wang
# Version: 1.0
# Email: wangm08@hotmail.com
# Date: 2023-08-14
#
# Usage: rnaseq.sh <config.txt>

################################################################################
# GLOBAL VARIABLES
SRC_DIR=$(dirname $(realpath -s $0)) # script dir
RNAseq_SH="${SRC_DIR}/run_rnaseq.sh" # <out_dir> <data_dir> [n_job]
RUN_SALMON="${SRC_DIR}/run_salmon.sh"
GB_COV="${SRC_DIR}/genebody_cov.sh" # to-do
################################################################################


# check fastq files in data_dir
# para: <data_dir>
function check_data_dir() {
    local data_dir=$1
    # check fq files
    local n_fq1=$(ls ${data_dir}/ | grep -c _1.fq.gz) #
    local n_fq2=$(ls ${data_dir}/ | grep -c _2.fq.gz) #
    if [[ ${n_fq1} -eq ${n_fq2} && ${n_fq1} -gt 0 ]] 
    then
        echo "[${n_fq1}] pair(s) of read1/read2 files found: ${data_dir}" && return 0 # pass
    else
        [[ ${n_fq1} -eq 0 ]] && echo "[0] no fastq files found: ${data_dir}" && return 1
        # check unpaired: fq1, fq2
        local tag=0 #
        for fq1 in ${data_dir}/*_1.fq.gz
        do
            fq2="${fq1%_1.fq.gz}_2.fq.gz"
            [[ ! -f ${fq2} ]] && echo "${fq1} -"
        done
        ##
        for fq2 in ${data_dir}/*_2.fq.gz
        do
            fq1="${fq1%_2.fq.gz}_1.fq.gz"
            [[ ! -f ${fq1} ]] && echo "- ${fq2}"
        done
        [[ ${tag} -ne 0 ]] && echo "[error] - above files not matched" && return 1
    fi
}
export -f check_data_dir


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
}
export -f load_config


# list config, global variables
# para: <>
function list_config() {
    # number of fastq files
    N_FQ=$(ls ${DATA_DIR}/ | grep -c _1.fq) #
    # show global variables:
    echo "--------------------------------------------------------------------------------"
    echo "WK_DIR:          ${WK_DIR}"
    echo "DATA_DIR:        ${DATA_DIR}"
    echo "GENOME:          ${GENOME}"
    echo "GENOME_DIR:      ${GENOME_DIR}"
    echo "BIN_SIZE:        ${BIN_SIZE}"
    echo "NORM_BY:         ${NORM_BY}"
    echo "N_CPU:           ${N_CPU}"
    echo "N_JOB:           ${N_JOB}"
    echo "N_MEM (MB):      ${N_MEM}"
    echo "slurm_partition: ${SLURM_PARTITION}"
    echo "n fastq files:   ${N_FQ}"
    echo ">> scripts <<"
    echo "run_rnaseq.sh    : ${RNAseq_SH}"
    echo "run_salmon.sh    : ${RUN_SALMON}"
    echo "genebody_cov.sh  : ${GB_COV}"
    echo "--------------------------------------------------------------------------------"
    echo ""
}
export -f list_config


# get prefix of fastq file
# para: <file>
function fx_prefix() {
    local fq=$1
    basename ${fq} | sed -Ee 's/(_[0-9]+)?.f(ast)?[aq](.gz)?$//i'
}
export -f fx_prefix


# info of SLURM
# para: <>
function slurm_info() {
    local part=$1 # partition
    # partition list
    local pt_list=$(sinfo -o "%P" -h | xargs | sed 's/*//')
    # node list
    local node_list=$(sinfo -o "%n" -h | xargs | sed 's/*//')
    # avail cpu (partition)
    local pt_cpus=$(${pt_list} | xargs -n 1 -I{} sinfo -o "%N %c" -p {} -h)
    # avail cpu (node)
    local node_cpus=$(${pt_list} | xargs -n 1 -I{} sinfo -o "%N %c" -n {} -h)
}
export -f slurm_info


# make slurm script for rnaseq analysis
# GLOBAL VARIABLES: SLURM_PARTITION, N_CPU, N_MEM, RNAseq_SH2
# para: <wk_dir> <data_dir>
function make_slurm() {
    [[ $# -lt 1 ]] && echo "Usage: make_slurm <config>" && return 1
    local config=$(realpath -s $1)
    load_config ${config} # global variables
    local slurm_file="${WK_DIR}/run_rnaseq_slurm.sh"    
    [[ ! -d ${WK_DIR} ]] && mkdir -p ${WK_DIR}
    # determine the total CPU/mem values
    local sum_CPU=$(echo ${N_CPU} ${N_JOB} | awk '{print $1*$2}')
    local sum_MEM=$(echo ${N_MEM} ${N_JOB} | awk '{print $1*$2}')

    cat << EOF > ${slurm_file}
#!/bin/bash
#SBATCH --job-name=rnaX          # create a short name for your job
#SBATCH --partition=${SLURM_PARTITION}         # node name
#SBATCH --ntasks=${sum_CPU}               # total cores, 
#SBATCH --mem=${sum_MEM}              # 30GB for STAR
#SBATCH --cpus-per-task=1         # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=5-01:00:00         # total run time limit (D-HH:MM:SS)
#SBATCH --output=log.${SLURM_PARTITION}.%j.out

module purge
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate hiseq

cd ${WK_DIR}
bash ${RNAseq_SH} ${config}
EOF
    # output
    echo ${slurm_file}
}
export -f make_slurm


# run RNAseq from terminal local, directly
# GLOBAL VARIABLES: WK_DIR, DATA_DIR, RNAseq_SH1, RNAseq_SH
# para: <config>
function make_local() {
    [[ $# -lt 1 ]] && echo "Usage: make_local <config>" && return 1
    local config=$(realpath -s $1)
    load_config ${config} # global variables
    local cmd_file="${WK_DIR}/run_rnaseq_local.sh"
    [[ ! -d ${WK_DIR} ]] && mkdir -p ${WK_DIR}
    echo "source $HOME/miniconda3/etc/profile.d/conda.sh" > ${cmd_file}
    echo "conda activate hiseq" >> ${cmd_file}
    echo "bash ${RNAseq_SH} ${config}" >> ${cmd_file}
    echo ${cmd_file}
}
export -f make_local


# RNAseq analysis
# para: <config> <run:0|1>
function rnaseq() {
    [[ $# -lt 2 ]] && echo "Usage: rnaseq <config.txt> [run:1|0]" && return 1
    local config=$1
    local run=$2 # 0=not, 1=run
    load_config ${config} # global variables
    list_config ${config} # init arguments, global variables
    check_data_dir "${DATA_DIR}" # global variables
    [[ $? -ne 0 ]] && echo "[error] - fq files vaild" && return 1 # previous command failed
    # prepare and run jobs
    [[ ${run} == 1 ]] && run_flag="yes" || run_flag="no" #
    if [[ -z $(sinfo -h -p ${SLURM_PARTITION}) ]] # global variables
    then
        # run on local server
        local cmd_file=$(make_local ${config})
        [[ $(grep -ci run_rnaseq.sh ${cmd_file}) -gt 0 ]] && job_flag="ok" || job_flag="failed"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${job_flag} "run" ${run_flag} $(basename ${cmd_file})
        [[ ${run} == 1 && ${job_flag} == "ok" ]] && bash ${cmd_file}
    else
        # run on HPC SLURM
        local slurm_file=$(make_slurm ${config})
        [[ $(grep -ci SBATCH ${slurm_file}) -gt 4 ]] && job_flag="ok" || job_flag="failed"
        printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${job_flag} "run" ${run_flag} $(basename ${slurm_file})
        [[ ${run} == 1 && ${job_flag} == "ok" ]] && sbatch ${slurm_file}
    fi
}
export -f rnaseq


rnaseq $@
