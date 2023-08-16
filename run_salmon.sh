# Run Salmon for gene-level quantification
# return
# 1. quant.sf (TPM)
# 2. quant.gene.tsv (TPM)
#
# Date: 2023-03-18

################################################################################
# GLOBAL VARIABLES
# salmon="/home/wangm/biosoft/salmon/current/bin/salmon"
salmon="/data1/yuyang/wangm/biosoft/salmon/current/bin/salmon"
################################################################################


# get the salmon index
# para: <genome>
function get_index() {
    # gdir="/data/yulab/wangming/data/genome/"
    gdir="${HOME}/data/genome"
    case $1 in 
        dm6|mm10|hg38)
            echo ${gdir}/${1}/salmon_index/${1}
            ;;
        *)
            echo 1
            ;;
    esac
}
export -f get_index


# prepare tx2gene Rscript
# para: <>
function tx2gene_R() {
    # convert transcript level to gene level
    # using scripts
    cat << EOF
## quant gene-level TPM
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1){
  print("Usage: Rscript run_salmon.R <quant.sf>")
  print("")
  print("Option:")
  print("  quant.sf     The output of salmon")
  stop("arguments failed")
}
sf <- args[1]

# locate the info.json
info <- file.path(dirname(sf), "cmd_info.json")
da   <- jsonlite::read_json(info)
idx  <- da\$index
t2g_f <- file.path(idx, "tx2gene.csv")

if(! file.exists(t2g_f)) {
  warning(glue::glue("t2g file not exists: {t2g_f}"))
} else {
  suppressPackageStartupMessages(library(tximport))
  tx2gene <- readr::read_csv(t2g_f, show_col_types = FALSE)
  txi <- tximport::tximport(sf, type = "salmon", tx2gene = tx2gene,
                            abundanceCol = "TPM")
  # gene, length, count, tpm
  df1 <- cbind(as.data.frame(txi\$length),
               as.data.frame(txi\$abundance),
               as.data.frame(txi\$counts))
  df1 <- round(df1, 4)
  colnames(df1) <- c("length", "TPM", "count")
  df1 <- tibble::rownames_to_column(df1, "id")
  # output
  gene_tpm <- file.path(dirname(sf), "quant.gene.tsv")
  message(glue::glue("save gene counts to file: {gene_tpm}"))
  readr::write_tsv(df1, gene_tpm, col_names = TRUE)
}
EOF
}
export -f tx2gene_R


# get fastq name
# para: <file>
function fx_prefix() {
    local fq=$1    
    # local prefix="$(basename ${fq/.gz})"
    # prefix=$(echo ${prefix} | sed -Ee 's/(_\d+)?.f(ast)?q(.gz)?//')
    basename ${fq/.gz} | sed -Ee 's/(_[0-9]+)?.f(ast)?[aq](.gz)?//i'
}
export -f fx_prefix


# run salmon alignment
# para: <out_dir> <genome> <fq1> <fq2> [n_cpu:4]
function run_salmon() {
    [[ $# -lt 4 ]] && echo "Usage: run_salmon.sh <out_dir> <genome> <fq1> <fq2> [n_cpu:4]" && return 1
    local out_dir=$1
    local genome=$2
    local fq1=$3
    local fq2=$4
    local n_cpu=$5
    [[ -z $5 ]] && n_cpu=4 # default value
    [[ ${n_cpu} =~ ^[0-9]+$ ]] || n_cpu=4 # default value
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    [[ ${fq1} == ${fq2} ]] && fq2="" # for single-end
    local prefix=$(fx_prefix ${fq1})
    local smp_dir="${out_dir}/${prefix}"
    [[ ! -d ${smp_dir} ]] && mkdir -p ${smp_dir}
    smp_dir=$(realpath -s ${smp_dir})
    # check index
    local idx=$(get_index ${genome})
    [[ ! -d ${idx} ]] && echo "index not exists: ${idx}" && return 1
    # check exists
    sf="${smp_dir}/quant.sf"
    log="${smp_dir}/log.stdout"
    # 1. run salmon
    if [[ -f ${sf} ]] ; then
        echo "file exists: ${sf}"
    else
        if [[ -f ${fq2} ]] ; then
            # Paired-end mode
            ${salmon} quant --gcBias --validateMappings -l A -p ${n_cpu} -i ${idx} -o ${smp_dir} -1 ${fq1} -2 ${fq2} 2> ${log}
        else
            # Single-end mode
            ${salmon} quant --gcBias --validateMappings -l A -p ${n_cpu} -i ${idx} -o ${smp_dir} -r ${fq1} 2> ${log}
        fi
    fi
    # 2. run tx2gene 
    run_t2g="${smp_dir}/run_tx2gene.sh"
    t2g_R="${smp_dir}/tx2gene.R"
    gg="${smp_dir}/quant.gene.tsv"
    if [[ -f ${gg} ]] ; then
        echo "file exists: ${gg}"
    else
        tx2gene_R > ${t2g_R}
        # echo "Rscript $(basename ${t2g_R}) $(basename ${sf})" > ${run_t2g}
        echo "Rscript $(basename ${t2g_R}) $(basename ${sf})" > ${run_t2g}
        pre_dir=$(pwd)
        cd ${smp_dir}
        bash ${run_t2g}
        cd ${pre_dir}
    fi
}
export -f run_salmon


run_salmon $@
