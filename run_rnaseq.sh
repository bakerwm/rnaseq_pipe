
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
#



## Global Variables
[[ -z ${N_CPU} ]] && N_CPU=8
[[ -z ${GENOME} ]] && GENOME=hg38 # mm10, hg38, dm6
[[ -z ${BIN_SIZE} ]] && BIN_SIZE=50 # for bigWig files
[[ -z ${NORM_BY} ]] && NORM_BY="CPM" # for bigWig
[[ -z ${GENOME_DIR} ]] && GENOME_DIR="${HOME}/data/genome"
SRC_DIR=$(dirname $0) # script dir
RUN_SALMON="${SRC_DIR}/run_salmon.sh"
GB_COV="${SRC_DIR}/genebody_cov.sh"


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


## prefix
function fx_prefix() {
    local fq=$1    
    # local prefix="$(basename ${fq/.gz})"
    # prefix=$(echo ${prefix} | sed -Ee 's/(_\d+)?.f(ast)?q(.gz)?//')
    basename ${fq%.gz} | sed -Ee 's/(_[0-9]+)?.f(ast)?[aq](.gz)?//i'
}
export -f fx_prefix


## trimming reads 
function trim_ad() {
    local fq1=$1
    local fq2=$2
    local wk_dir=$3 # working directory
    local clean_dir=${wk_dir}/data/clean_data
    # local raw_dir=$(dirname ${fq1})
    # local clean_dir="${wk_dir}/data/clean_data"
    # hiseq trim -n 4 -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTG -p ${N_CPU} -j 1 --cut-after-trim 7,-7 -m 20 -o ${clean_dir} -1 ${fq1} -2 ${fq2}
    hiseq trim --rm-polyN A --times 4 -p ${N_CPU} -p ${N_CPU} -m 20 -o ${clean_dir} -1 ${fq1} -2 ${fq2}
    # hiseq qc -i data/clean_data/${fq_name}/*gz -o ${clean_dir}/qc -p ${N_CPU} -j 1 -c fastqc
    # hiseq qc -i ${fq1} ${fq2} -o ${raw_dir}/qc -p 4 -j 1 -c fastqc
}
export -f trim_ad


## bam to bigWig
function bam2bw() {
    local out_dir=$1
    local bam=$2
    local genome=$3
    # local scale=$4
    local bname=$(basename ${bam%.bam}) # trim from right
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
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
    cmd="${out_dir}/${bname}.cmd.sh"
    ## fwd
    local bw_fwd="${out_dir}/${bname}_fwd.bigWig"
    [[ -f ${bw_fwd} ]] && tag_fwd="# " || tag_fwd=""
    cmd_fwd="${tag_fwd}bamCoverage -p ${N_CPU} -b ${bam} -o ${bw_fwd} --filterRNAstrand forward --binSize ${BIN_SIZE} --effectiveGenomeSize ${gsize} --scaleFactor 1.0 --normalizeUsing ${NORM_BY} --skipNAs --centerReads --smoothLength 150"
    # [[ -f ${bw_fwd} ]] && echo "# file exists: ${bw_fwd}" > ${cmd} || echo ${cmd_fwd} > ${cmd}    
    ## rev
    local bw_rev="${out_dir}/${bname}_rev.bigWig"
    [[ -f ${bw_rev} ]] && tag_rev="# " || tag_rev=""
    cmd_rev="${tag_rev}bamCoverage -p ${N_CPU} -b ${bam} -o ${bw_rev} --filterRNAstrand reverse --binSize ${BIN_SIZE} --effectiveGenomeSize ${gsize} --scaleFactor 1.0 --normalizeUsing ${NORM_BY} --skipNAs --centerReads --smoothLength 150"
    # [[ -f ${bw_rev} ]] && echo "# file exists: ${bw_rev}" >> ${cmd} || echo ${cmd_rev} >> ${cmd}
    ## run command 
    echo ${cmd_fwd} > ${cmd}
    echo ${cmd_rev} >> ${cmd}
    bash ${cmd}
}
export -f bam2bw


# http://broadinstitute.github.io/picard/index.html#GettingHelp
function anno() {
    local out_dir=$1
    local bam=$2
    local bname=$(basename ${bam%.bam})
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local anno_stat=${out_dir}/${bname}.anno.stat
    local anno_log=${out_dir}/${bname}.anno.log
    ## REF_FLAT, see global variables
    local REF_FLAT="${GENOME_DIR}/${GENOME}/annotation_and_repeats/${GENOME}.refFlat.txt" # hg38.refFlat.txt
    cmd="${out_dir}/${bname}.cmd.sh"
    # rRNA="${genome}.rRNA.interval_list"
    [[ ! -f ${REF_FLAT} ]] && echo "file not exists: ${REF_FLAT}" && exit 1
    [[ -f ${anno_stat} ]] && tag="# " || tag=""
    cmd_shell="${tag}picard CollectRnaSeqMetrics -I ${bam} -O ${anno_stat} -REF_FLAT ${REF_FLAT} -STRAND FIRST_READ_TRANSCRIPTION_STRAND >${anno_log}"
    echo ${cmd_shell} > ${cmd}
    bash ${cmd}
}
export -f anno


## count reads
function quant() {
    local out_dir=$1
    local genome=$2
    local bam=$3
    local bname=$(basename ${bam%.bam})
    # gtf=$1
    local group="exon"
    local gtf=$(get_gtf ${genome} ${group})
    ## parameters: paired-end
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local para="-p -T ${N_CPU} -M -O --fraction -F GTF -t gene -g gene_id"
    local cmd="${out_dir}/${bname}.cmd.sh"
    local prefix="${out_dir}/${bname}.${group}"
    ## sense
    local quant_sens="${prefix}.sens.txt"
    local quant_log_sens="${prefix}.sens.featureCounts.stderr"
    [[ -f ${quant_sens} ]] && tag_sens="# " || tag_sens=""
    cmd_sens="${tag_sens}featureCounts -s 1 ${para} -a ${gtf} -o ${quant_sens} ${bam} 2> ${quant_log_sens}"
    ## anti
    local quant_anti="${prefix}.anti.txt"
    local quant_log_anti="${prefix}.anti.featureCounts.stderr"
    [[ -f ${quant_anti} ]] && tag_anti="# " || tag_anti=""
    cmd_anti="${tag_anti}featureCounts -s 2 ${para} -a ${gtf} -o ${quant_anti} ${bam} 2> ${quant_log_anti}"
    # save to cmd
    echo ${cmd_sens} > ${cmd}
    echo ${cmd_anti} >> ${cmd}
    bash ${cmd}
}
export -f quant


## run for single fastq file: fq1
function run_salmon() {
    local out_dir=$1
    local genome=$2
    local fq1=$3
    local fq2=$4
    bash ${RUN_SALMON} ${out_dir} ${genome} ${fq1} ${fq2}
}
export -f run_salmon


# para: <wk_dir> <fq1> <fq2>
function run_rnaseq() {
    local wk_dir=$1 # working directory
    local fq1=$2
    local fq2=$3
    [[ $# -lt 3 ]] && echo "not enough args..." && return 1
    [[ ! -f ${fq1} ]] && echo "File not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "File not exists: ${fq2}" && return 1
    wk_dir=$(realpath -s ${wk_dir}) # absolute path
    local fq_name=$(basename ${fq1%_1.fq.gz})

    ## 0. check env
    if [[ ! ${CONDA_DEFAULT_ENV} == "hiseq" ]] 
    then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        conda activate hiseq
    fi
    
    ## 1. trim reads
    # read1: cut 7  (before trim)
    # read2: cut -7 (after trim)
    echo "[1/10] - trimming adapters"
    local raw_dir=$(dirname ${fq1})
    local clean_dir="${wk_dir}/results/00.clean_data"
    ## hiseq trim -n 4 -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTG -p ${N_CPU} -j 1 --cut-after-trim 7,-7 -m 20 -o ${clean_dir} -1 ${fq1} -2 ${fq2}
    hiseq trim --rm-polyN A --times 4 -p ${N_CPU} -j 1 -m 15 -o ${clean_dir} -1 ${fq1} -2 ${fq2}
    ## hiseq qc -i data/clean_data/${fq_name}/*gz -o ${clean_dir}/qc -p ${N_CPU} -j 1 -c fastqc
    ## hiseq qc -i ${fq1} ${fq2} -o ${raw_dir}/qc -p 4 -j 1 -c fastqc
    local clean_fq1="${clean_dir}/${fq_name}/${fq_name}_1.fq.gz"
    local clean_fq2="${clean_dir}/${fq_name}/${fq_name}_2.fq.gz"

    ## 2. align to genome, using STAR
    echo "[2/10] - mapping (STAR)"
    local align_dir="${wk_dir}/results/01.align"
    local bam_star="${align_dir}/${fq_name}/${fq_name}.bam"
    [[ ! -f ${bam_star} ]] && \
        hiseq align -a STAR -o ${align_dir} --to-rRNA -p ${N_CPU} -j 1 -g ${GENOME} -1 ${clean_fq1} -2 ${clean_fq2}
    [[ ! -f ${bam_star} ]] && echo "    [error] - bam not file found"

    ## 3. prepare bam files
    echo "[3/10] - copy bam files"
    local bam_dir="${wk_dir}/results/02.bam_files"
    local bam="${bam_dir}/${fq_name}.bam"
    local bam_star_rel=$(realpath --relative-to=${bam_dir} $(dirname ${bam_star}))
    [[ ! -d ${bam_dir} ]] && mkdir -p ${bam_dir}
    ln -fs ${bam_star_rel}/${fq_name}.bam ${bam}
    [[ ! -f ${bam}.bai ]] && samtools index -@ {N_CPU} ${bam}

    ## 4. prepare bw files
    echo "[4/10] - generating bigWig files"
    local bw_dir="${wk_dir}/results/03.bw_files"
    bam2bw ${bw_dir} ${bam} ${GENOME}

    ## 5. annotation bam files
    echo "[5/10] - annotation (Picard)"
    local anno_dir="${wk_dir}/results/04.anno/picard"
    anno ${anno_dir} ${bam}

    ## 6. quant
    echo "[6/10] - quantifying sens,anti reads"
    local quant_dir="${wk_dir}/results/05.quant"
    quant ${quant_dir} ${GENOME} ${bam}

    ## 7. genebody coverage
    echo "[7/10] - generating genebody coverage (RSeQC), skipped ..."
    local gb_dir="${wk_dir}/results/06.genebody_cov"
    # bash ${GB_COV} ${gb_dir} ${GENOME} 1 ${bam}

    ## 8. quant, TPM, gene_count, salmon
    echo "[8/10] - mapping (Salmon)"
    local salmon_dir="${wk_dir}/results/07.RNAseq_salmon"
    run_salmon ${salmon_dir} ${GENOME} ${clean_fq1} ${clean_fq2}

    ## 9. subsample: 0.1M to 1M reads by 0.2M step
    echo "[9/10] - subset reads by 0.1 to 1.0 M"
    local sub_dir="${wk_dir}/results/08.sub_data"
    for i in 100000 200000 400000 600000 800000 1000000
    do
        i2=$(awk -v n=${i} 'BEGIN{printf "%02d",n/100000}') # fix name
        local sub_data_dir="${sub_dir}/sub_${i2}/clean_data" # suffix
        echo ${i2}
        hiseq sample -j 1 -n ${i} -o ${sub_data_dir} -i ${clean_fq1} ${clean_fq2}
        # run salmon 
        local sub_fq1=${sub_data_dir}/${fq_name}_1.fq.gz
        local sub_fq2=${sub_data_dir}/${fq_name}_2.fq.gz
        local sub_align_dir="${sub_dir}/sub_${i2}/"
        run_salmon ${sub_align_dir} ${GENOME} ${sub_fq1} ${sub_fq2}
    done
    
    echo "[10/10] - finished!"
}
export -f run_rnaseq


[[ $# -lt 3 ]] && echo "bash run_rnaseq.sh <wk_dir> <fq1_1.fq.gz> <fq2_2.fq.gz>" && exit 1
echo "working on $$"
run_rnaseq $1 $2 $3

