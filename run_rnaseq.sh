#!/usr/bin/bash

# Script: run_rnaseq.sh
# Author: Ming Wang
# Version: 1.0
# Email: wangm08@hotmail.com
# Date: 2023-08-14
#
# Usage: run_rnaseq.sh <out_dir> <data_dir> <n_job>

# Description
#
# Analysis for RNAseq samples
#
# Requirements
# 1. paird-end files, named by ".fq.gz"
# 2. `parallel` >20230522
#
# Getting Started
# 1. trimming.   Remove adapters, low quality bases at 3' end, polyA,T
# 2. mapping.    mapping reads to reference genome
# 3. bigWig.     Generate bigwig for each sample
# 4. annotation. Using Picard for annotation
# 6. genebody.   Using RseQC for genebody coverage analysis
# 7. strandness. Using featureCounts to count reads on sens and anti strand
# 5. TPM count.  Using Salmon to quantify gene-level TPM (Salmon >=1.10.0)
# 8. genecount.  Using Salmon, TPM to count genes
# 9. subsample.  Subsample fastq file by 0.1 to 1.0 M reads, for gene count
#
# Directory structure
# results/
# ├── 00.clean_data
# ├── 01.align
# ├── 02.bam_files
# ├── 03.bw_files
# ├── 04.anno
# ├── 05.quant
# ├── 06.genebody_cov
# ├── 07.RNAseq_salmon
# └── 08.sub_data

################################################################################
## Global Variables, parsing from file: ./rnaseq.sh and config.txt
SRC_DIR=$(dirname $0) # script dir
RNAseq_SH="${SRC_DIR}/run_rnaseq.sh" # <out_dir> <data_dir> [n_job]
RUN_SALMON="${SRC_DIR}/run_salmon.sh"
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
    REF_FLAT="${GENOME_DIR}/${GENOME}/annotation_and_repeats/${GENOME}.refFlat.txt" # hg38.refFlat.txt
}
export -f load_config


# list config, global variables
# para: <>
function list_config() {
    # number of fastq files
    N_FQ=$(ls ${DATA_DIR}/ | grep -c _1.fq) #
    # show global variables:
    echo "--------------------------------------------------------------------------------"
    echo "WK_DIR           : ${WK_DIR}"
    echo "DATA_DIR         : ${DATA_DIR}"
    echo "GENOME           : ${GENOME}"
    echo "GENOME_DIR       : ${GENOME_DIR}"
    echo "BIN_SIZE         : ${BIN_SIZE}"
    echo "NORM_BY          : ${NORM_BY}"
    echo "N_CPU            : ${N_CPU}"
    echo "N_JOB            : ${N_JOB}"
    echo "N_MEM (MB)       : ${N_MEM}"
    echo "slurm_partition  : ${SLURM_PARTITION}"
    echo "SCRIPT_DIR       : ${SCRIPT_DIR}"
    echo "run_rnaseq.sh    : ${RNAseq_SH}"
    echo "run_salmon.sh    : ${RUN_SALMON}"
    echo "genebody_cov.sh  : ${GB_COV}"
    echo "N_fastq_files    : ${N_FQ}"
    echo "--------------------------------------------------------------------------------"
    echo ""
}
export -f list_config


# current time
# para: <>
function get_curr_time() {
    # formatted_time=$(date +"%Y-%m-%d %H:%M:%S")
    date +"%Y-%m-%d %H:%M:%S"
}
export -f get_curr_time


## Genome info for featureCounts
function get_gtf() {
    local group=$1 # exon,intron,3u,5u
    # [[ -z ${group} ]] && group="exon"
    case ${group} in
        exon|intron|3u|5u)
            group=${group} ;;
        *)
            group="exon" ;;
    esac
    echo ${GENOME_DIR}/${GENOME}/annotation_and_repeats/refseq_gtf/${GENOME}.refseq.exon.gtf
}
export -f get_gtf


# prefix
# para: <file>
function fx_prefix() {
    local fq=$1
    # local prefix="$(basename ${fq/.gz})"
    # prefix=$(echo ${prefix} | sed -Ee 's/(_\d+)?.f(ast)?q(.gz)?//')
    basename ${fq%.gz} | sed -Ee 's/(_[0-9]+)?.f(ast)?[aq](.gz)?//i'
}
export -f fx_prefix


# check fastq files
# para: <data_dir>
function check_data_dir() {
    local data_dir=$1
    local n_fq1=$(ls ${data_dir}/ | grep -c _1.fq.gz) #
    local n_fq2=$(ls ${data_dir}/ | grep -c _2.fq.gz) #
    if [[ ${n_fq1} -eq ${n_fq2} && ${n_fq1} -gt 0 ]]
    then
        echo "[${n_fq1}] pair(s) of read1/read2 files found: ${data_dir}" && return 0 # pass
    else
        [[ ${n_fq1} -eq 0 ]] && echo "[0] no fastq files found: ${data_dir}" && return 1
        # check unpaired: fq1, fq2
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
        echo "[error] - above files not matched" && return 1
    fi
}
export -f check_data_dir


# trimming reads
# para: <out_dir> <fq1> <fq2>
function trim_ad() {
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    local fq_name=$(fx_prefix ${fq1})
    local cmd="${out_dir}/${fq_name}.cmd.sh"
    local trim_out="${out_dir}/${fq_name}.trim.stdout"
    local trim_err="${out_dir}/${fq_name}.trim.stderr"
    # auto detect adapters, by hiseq (TruSeq, Nextera, smallRNA)
    echo "hiseq trim --rm-polyN A --times 4 -p ${N_CPU} -j 1 \
        -m 20 -o ${out_dir} -1 ${fq1} -2 ${fq2} \
        1>${trim_out} 2>${trim_err}" > ${cmd}
    # run command
    bash ${cmd}
    ## specify adapters, for XP seq, NSR+Tn5
    ## hiseq trim -n 4 -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTG -p ${N_CPU} -j 1 --cut-after-trim 7,-7 -m 20 -o ${clean_dir} -1 ${fq1} -2 ${fq2}
    ## hiseq qc -i data/clean_data/${fq_name}/*gz -o ${clean_dir}/qc -p ${N_CPU} -j 1 -c fastqc
    ## hiseq qc -i ${fq1} ${fq2} -o ${raw_dir}/qc -p 4 -j 1 -c fastqc
    # local clean_fq1="${clean_dir}/${fq_name}/${fq_name}_1.fq.gz"
    # local clean_fq2="${clean_dir}/${fq_name}/${fq_name}_2.fq.gz"
}
export -f trim_ad


# run STAR for single fastq file
# para: <out_dir> <genome> <fq1> <fq2>
function align_star() {
    local out_dir=$1
    local genome=$2
    local fq1=$3
    local fq2=$4
    local fname=$(fx_prefix ${fq1})
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    local cmd="${out_dir}/${fname}.cmd.sh"
    local align_out="${out_dir}/${fname}.star.stdout"
    local align_err="${out_dir}/${fname}.star.stderr"
    local bam_star="${out_dir}/${fname}/${fname}.bam"
    [[ -f ${bam_star} ]] && tag_bam="# " || tag_bam=""
    echo "${tag_bam}hiseq align -a STAR -o ${align_dir} --to-rRNA \
        -p ${N_CPU} -j 1 -g ${GENOME} \
        -1 ${clean_fq1} -2 ${clean_fq2} \
        1>${align_out} 2>${align_err}" > ${cmd}
    # run command
    bash ${cmd}
    [[ ! -f ${bam_star} ]] && echo "    [error] - bam not file found"
}
export -f align_star


# bam to bigWig
# para: <out_dir> <genome> <bam>
function bam2bw() {
    local out_dir=$1
    local genome=$2
    local bam=$3
    # local scale=$4
    local bname=$(basename ${bam%.bam}) # trim from right
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local log_out="${out_dir}/${bname}.deeptools.stdout"
    local log_err="${out_dir}/${bname}.deeptools.stderr"
    ## case
    case ${genome} in
        hg38|human|hs)
            gsize="2913022398" ;;
        mm10|mouse|ms)
            gsize="2652783500" ;;
        dm6|fly|dm)
            gsize="142573017" ;;
        *)
            gsize=1
            echo "Unknown genome, expect [dm6,mm10,hg38]," ;;
    esac
    # local gsize="2652783500" # mm10
    # local gsize="142573017"  # dm6
    # local gsize="2913022398" # hg38
    local cmd="${out_dir}/${bname}.cmd.sh"
    ## fwd
    local bw_fwd="${out_dir}/${bname}_fwd.bigWig"
    [[ -f ${bw_fwd} ]] && tag_fwd="# " || tag_fwd=""
    echo "${tag_fwd}bamCoverage -p ${N_CPU} -b ${bam} -o ${bw_fwd} \
        --filterRNAstrand forward --binSize ${BIN_SIZE} \
        --effectiveGenomeSize ${gsize} --scaleFactor 1.0 \
        --normalizeUsing ${NORM_BY} --skipNAs \
        --centerReads --smoothLength 150 \
        1>${log_out} 2>${log_err}" > ${cmd}
    # [[ -f ${bw_fwd} ]] && echo "# file exists: ${bw_fwd}" > ${cmd} || echo ${cmd_fwd} > ${cmd}
    ## rev
    local bw_rev="${out_dir}/${bname}_rev.bigWig"
    [[ -f ${bw_rev} ]] && tag_rev="# " || tag_rev=""
    echo "${tag_rev}bamCoverage -p ${N_CPU} -b ${bam} -o ${bw_rev} \
        --filterRNAstrand reverse --binSize ${BIN_SIZE} \
        --effectiveGenomeSize ${gsize} --scaleFactor 1.0 \
        --normalizeUsing ${NORM_BY} --skipNAs \
        --centerReads --smoothLength 150 \
        1>>${log_out} 2>>${log_err}" >> ${cmd}
    # [[ -f ${bw_rev} ]] && echo "# file exists: ${bw_rev}" >> ${cmd} || echo ${cmd_rev} >> ${cmd}
    ## run command
    bash ${cmd}
}
export -f bam2bw


# see: http://broadinstitute.github.io/picard/index.html#GettingHelp
# para: <out_dir> <bam>
function anno() {
    local out_dir=$1
    local bam=$2
    local bname=$(basename ${bam%.bam})
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local anno_stat=${out_dir}/${bname}.anno.stat
    local log_out="${out_dir}/${bname}.picard.stdout"
    local log_err="${out_dir}/${bname}.picard.stderr"
    # REF_FLAT, see global variable
    local cmd="${out_dir}/${bname}.cmd.sh"
    # rRNA="${genome}.rRNA.interval_list"
    [[ ! -f ${REF_FLAT} ]] && echo "file not exists: ${REF_FLAT}" && return 1
    [[ -f ${anno_stat} ]] && tag="# " || tag=""
    echo "${tag}picard CollectRnaSeqMetrics -I ${bam} -O ${anno_stat} \
        -REF_FLAT ${REF_FLAT} -STRAND FIRST_READ_TRANSCRIPTION_STRAND \
        1>${log_out} 2>${log_err}" > ${cmd}
    # run command
    bash ${cmd}
}
export -f anno


# count reads using featureCounts
# para: <out_dir> <genome> <bam>
function quant() {
    local out_dir=$1
    local genome=$2
    local bam=$3
    local bname=$(basename ${bam%.bam})
    local group="exon"
    local gtf=$(get_gtf ${genome} ${group})
    ## parameters: paired-end
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local para="-p -T ${N_CPU} -M -O --fraction -F GTF -t gene -g gene_id"
    local cmd="${out_dir}/${bname}.cmd.sh"
    local prefix="${out_dir}/${bname}.${group}"
    local log_out="${out_dir}/${bname}.featureCounts.stdout"
    local log_err="${out_dir}/${bname}.featureCounts.stderr"
    ## sense
    local quant_sens="${prefix}.sens.txt"
    # local quant_log_sens="${prefix}.sens.featureCounts.stderr"
    [[ -f ${quant_sens} ]] && tag_sens="# " || tag_sens=""
    echo "${tag_sens}featureCounts -s 1 ${para} -a ${gtf} \
        -o ${quant_sens} ${bam} 1>${log_out} 2>${log_err}" > ${cmd}
    ## anti
    local quant_anti="${prefix}.anti.txt"
    local quant_log_anti="${prefix}.anti.featureCounts.stderr"
    [[ -f ${quant_anti} ]] && tag_anti="# " || tag_anti=""
    echo "${tag_anti}featureCounts -s 2 ${para} -a ${gtf} \
        -o ${quant_anti} ${bam} 1>>${log_out} 2>>${log_err}" >> ${cmd}
    # run command
    bash ${cmd}
}
export -f quant


# run salmon for single fastq file: fq1
# para: <run_salmon> <out_dir> <genome> <fq1> <fq2>
function align_salmon() {
    [[ $# -lt 6 ]] && echo "Usage: align_salmon <run_salmon> <out_dir> <genome> <fq1> <fq2> <n_cpu>" && return 1
    local run_salmon=$1 # script: run_salmon.sh
    local out_dir=$2
    local genome=$3
    local fq1=$4
    local fq2=$5
    local n_cpu=$6
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    local fname=$(fx_prefix ${fq1})
    local cmd="${out_dir}/${fname}.cmd.sh"
    local log_out="${out_dir}/${fname}.slamon.stdout"
    local log_err="${out_dir}/${fname}.slamon.stderr"
    local quant_sf="${out_dir}/${fname}/quant.sf"
    [[ -f ${quant_sf} ]] && tag_sf="# " || tag_sf=""
    echo "${tag_sf}bash ${run_salmon} ${out_dir} ${genome} ${fq1} ${fq2} ${n_cpu} \
        1>${log_out} 2>${log_err}" > ${cmd}
    # run command
    bash ${cmd}
}
export -f align_salmon


# run pipeline for single pair of read1/2
# para: <config> <fq1> <fq2>
function run_rnaseq_fq() {
    [[ $# -lt 3 ]] && echo "run_rnaseq_fq <config> <fq1> <fq2>" && return 1
    local config=$1  # for global variables
    local fq1=$2
    local fq2=$3
    local fq_name=$(fx_prefix ${fq1})
    [[ ! -f ${fq1} ]] && echo "File not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "File not exists: ${fq2}" && return 1
    #####################################################################
    # Error-1. could not detect global variables using parallel.        #
    # in case, to use parallel for running multiple samples             #
    # functions within paralle could not detect global variables.       #
    #                                                                   #
    # Solution. pass the config to this function, within parallel       #
    #####################################################################
    load_config ${config} # global variables !!!
    # global variables could not be reached by functions within parallel
    local run_salmon="${SCRIPT_DIR}/run_salmon.sh" # local var within this process
    local out_dir=${WK_DIR}/results

    ## 0. init env
    if [[ ! ${CONDA_DEFAULT_ENV} == "hiseq" ]]
    then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        conda activate hiseq
    fi

    echo "[$(get_curr_time) 1/10] - trimming adapters [${fq_name}]"
    local clean_dir="${out_dir}/00.clean_data"
    trim_ad ${clean_dir} ${fq1} ${fq2}
    local clean_fq1="${clean_dir}/${fq_name}/${fq_name}_1.fq.gz"
    local clean_fq2="${clean_dir}/${fq_name}/${fq_name}_2.fq.gz"

    echo "[$(get_curr_time) 2/10] - mapping (STAR) [${fq_name}]"
    local align_dir="${out_dir}/01.align"
    align_star ${align_dir} ${GENOME} ${clean_fq1} ${clean_fq2}

    echo "[$(get_curr_time) 3/10] - preparing bam file link [${fq_name}]"
    local bam_dir="${out_dir}/02.bam_files"
    [[ ! -d ${bam_dir} ]] && mkdir -p ${bam_dir}
    local bam="${bam_dir}/${fq_name}.bam"
    local bam_star="${align_dir}/${fq_name}/${fq_name}.bam"
    local bam_star_rel=$(realpath --relative-to=${bam_dir} ${align_dir}/${fq_name}/)
    ln -fs ${bam_star_rel}/${fq_name}.bam ${bam}
    [[ ! -f ${bam}.bai ]] && samtools index -@ {N_CPU} ${bam}

    echo "[$(get_curr_time) 4/10] - generating bigWig files [${fq_name}]"
    local bw_dir="${out_dir}/03.bw_files"
    bam2bw ${bw_dir} ${GENOME} ${bam}

    echo "[$(get_curr_time) 5/10] - running annotation (Picard) [${fq_name}]"
    local anno_dir="${out_dir}/04.anno/picard"
    anno ${anno_dir} ${bam}

    echo "[$(get_curr_time) 6/10] - quantifying sens,anti reads [${fq_name}]"
    local quant_dir="${out_dir}/05.quant"
    quant ${quant_dir} ${GENOME} ${bam}

    echo "[$(get_curr_time) 7/10] - generating genebody coverage (RSeQC), skipped ..."
    local gb_dir="${out_dir}/06.genebody_cov"
    ## bash ${GB_COV} ${gb_dir} ${GENOME} 1 ${bam} # skipped, run this step alone

    echo "[$(get_curr_time) 8/10] - mapping (Salmon) [${fq_name}]"
    local salmon_dir="${out_dir}/07.RNAseq_salmon"
    align_salmon ${run_salmon} ${salmon_dir} ${GENOME} ${clean_fq1} ${clean_fq2} ${N_CPU}

    echo "[$(get_curr_time) 9/10] - mapping subset reads (Salmon), for 0.1 to 1.0 M [${fq_name}]"
    local sub_dir="${out_dir}/08.sub_data"
    for i in 100000 200000 400000 600000 800000 1000000
    do
        i2=$(awk -v n=${i} 'BEGIN{printf "%02d",n/100000}') # fix name
        local sub_data_dir="${sub_dir}/sub_${i2}/clean_data" # suffix
        echo ${i2}
        hiseq sample -j 1 -n ${i} -o ${sub_data_dir} -i ${clean_fq1} ${clean_fq2}
        # run salmon
        local sub_align_dir="${sub_dir}/sub_${i2}/"
        local sub_fq1=${sub_data_dir}/${fq_name}_1.fq.gz
        local sub_fq2=${sub_data_dir}/${fq_name}_2.fq.gz
        align_salmon ${run_salmon} ${sub_align_dir} ${GENOME} ${sub_fq1} ${sub_fq2} ${N_CPU}
    done

    echo "[$(get_curr_time) 10/10] - finished! [${fq_name}]"
}
export -f run_rnaseq_fq


# run pipeline for multiple files
# para: <out_dir> <data_dir>
function run_rnaseq() {
    [[ $# -lt 1 ]] && echo "run_rnaseq <config>" && return 1
    echo "working on $$"
    local config=$1
    load_config ${config} # global variables
    list_config ${config} # init arguments, global variables
    tmp=$(check_data_dir "${DATA_DIR}") # global variables
    [[ $? -ne 0 ]] && echo "[error] - fq files vaild" && return 1 # previous command failed
    # prepare subdirectories
    # in-case, parallel create directories in the same time
    local out_dir="${WK_DIR}/results"
    mkdir -p ${out_dir}/00.clean_data  ${out_dir}/01.align  ${out_dir}/02.bam_files \
        ${out_dir}/03.bw_files  ${out_dir}/04.anno/picard  ${out_dir}/05.quant  \
        ${out_dir}/06.genebody_cov  ${out_dir}/07.RNAseq_salmon  \
        ${out_dir}/08.sub_data/sub_01
    # list_config ${config} # init arguments, global variables
    parallel -j ${N_JOB} --rpl '{/(.+?)/(.*?)} s/$$1/$$2/;' run_rnaseq_fq ${config} {} {/_1.fq/_2.fq} ::: ${DATA_DIR}/*1.fq.gz
    # for fq1 in ${DATA_DIR}/*1.fq.gz
    # do
    #     fq2="${fq1%_1.fq.gz}_2.fq.gz"
    #     run_rnaseq_fq ${config} ${fq1} ${fq2}
    #     break
    # done
}
export -f run_rnaseq

# main
run_rnaseq $@
