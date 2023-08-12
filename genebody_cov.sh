
# Genebody coverage for BAM/bigWig files 
#
# Using subcommands from RSeQC:
# - BAM: geneBody_coverage.py 
# - BigWIG: geneBody_coverage2.py
#

## geneBody_coverage.py
## Usage: 
# Options:
#   --version             show program's version number and exit
#   -h, --help            show this help message and exit
#   -i INPUT_FILES, --input=INPUT_FILES
#                         Input file(s) in BAM format. "-i" takes these input:
#                         1) a single BAM file. 2) "," separated BAM files. 3)
#                         directory containing one or more bam files. 4) plain
#                         text file containing the path of one or more bam file
#                         (Each row is a BAM file path). All BAM files should be
#                         sorted and indexed using samtools.
#   -r REF_GENE_MODEL, --refgene=REF_GENE_MODEL
#                         Reference gene model in bed format. [required]
#   -l MIN_MRNA_LENGTH, --minimum_length=MIN_MRNA_LENGTH
#                         Minimum mRNA length (bp). mRNA smaller than
#                         "min_mRNA_length" will be skipped. default=100
#   -f OUTPUT_FORMAT, --format=OUTPUT_FORMAT
#                         Output file format, 'pdf', 'png' or 'jpeg'.
#                         default=pdf
#   -o OUTPUT_PREFIX, --out-prefix=OUTPUT_PREFIX
#                         Prefix of output files(s). [required]


function prep_bed() {
    cat << EOF
Pre-build gene models could bd downloaded from: https://sourceforge.net/projects/rseqc/files/BED

Prepare bed file for RSeQC, geneBody_coverage.py command, with kent utils: https://hgdownload.soe.ucsc.edu/admin/exe/

# 1. Download GTF files from ENSEMBL (eg: release-103)
$ wget https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz 

# 2. Convert GTF to BED
$ gtfToGenePred <(gunzip -c Homo_sapiens.GRCh38.103.gtf.gz) Homo_sapiens.GRCh38.103.genePred
$ genePredToBed Homo_sapiens.GRCh38.103.genePred Homo_sapiens.GRCh38.103.genePred.bed
## if I want to convert the chr name to UCSC format
$ chromToUcsc --get hg38 # download hg38.chromAlias.tsv
$ chromToUcsc -a hg38.chromAlias.tsv -i Homo_sapiens.GRCh38.103.genePred.bed -o Homo_sapiens.GRCh38.103.genePred.ucsc.bed

# 3. Check the bed files
$ head -n 2 Homo_sapiens.GRCh38.103.genePred.bed
1       11868   14409   ENST00000456328 0       +       14409   14409   0       3       359,109,1189,   0,744,1352,
1       12009   13670   ENST00000450305 0       +       13670   13670   0       6       48,49,85,78,154,218,    0,169,603,965,1211,1443,

$ head -n 2 Homo_sapiens.GRCh38.103.genePred.ucsc.bed
chr1    11868   14409   ENST00000456328 0       +       14409   14409   0       3       359,109,1189,   0,744,1352,
chr1    12009   13670   ENST00000450305 0       +       13670   13670   0       6       48,49,85,78,154,218,    0,169,603,965,1211,1443,

$ wc -l Homo_sapiens.GRCh38.103.genePred.ucsc.bed
234393 Homo_sapiens.GRCh38.103.genePred.ucsc.bed
EOF
}
export -f prep_bed 


function get_bed() {
    local genome=$1
    local bed=${HOME}/data/genome/${genome}/annotation_and_repeats/gene_bed/${genome}.gene.bed
    # [[ ! -f {bed} ]] && prep_bed && return 1
    echo ${bed}
}
export -f get_bed


function setup_env() {
    # create env for rseqc
    cat << EOF
# Prepare conda env for RSeQC
$ conda create -n rseqc -c bioconda rseqc
$ conda activate rseqc
$ geneBody_coverage.py -h
EOF
}
export -f setup_env


function check_env() {
    ## for conda env
    # Get the current Conda environment name
    conda_env_name="${CONDA_DEFAULT_ENV}"
    # Check if a Conda environment is activated
    if [[ ${conda_env_name} == "rseqc" ]] 
    then
        echo "Current Conda environment: ${conda_env_name}"
    else
        # source ${HOME}/miniconda3/bin/activate
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        conda activate rseqc
        # geneBody_coverage.py -h
    fi
    # commands:
    if [[ -z $(command -v geneBody_coverage.py) || -z $(command -v geneBody_coverage2.py) ]] 
    then 
        echo "RSeQC commands not found, follow the instructions to install RSeQC:"
        setup_env # usage
        return 1
    fi
}
export -f check_env

# for single bam
# para: <out_dir> <genome> <bam>
function gbcov() {
    [[ $# -lt 3 ]] && echo "Usage: gbcov <out_dir> <genome> <bam>" && return 1
    # check_env
    [[ ! -z $(check_env) ]] && return 1 #
    check_env # switch to rseqc
    # genebody coverage
    local out_dir=$1
    local genome=$2
    local bam=$3 # could be bam or bigWig
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    local bed=$(get_bed ${genome})
    [[ ! -f ${bed} ]] && echo "bed file not found: ${bed}" && return 1
    ## check if bam or bigWig
    if [[ ${bam} = *.bam ]] 
    then 
        local bname=$(basename ${bam/.bam})
        local gb_cmd="geneBody_coverage.py"
    elif [[ ${bam} = *.bigWig || ${bam}  = *.bw ]] 
    then 
        local bname=$(basename ${bam/.bigWig})
        banme=${bname/.bw}
        local gb_cmd="geneBody_coverage.py"
    else
        echo "unknown bam file: ${bam}"
        return 1
    fi
    out_dir=$(realpath -s ${out_dir})
    bam=$(realpath -s ${bam})
    ## prep command
    local prefix="${out_dir}/${bname}"
    local out_txt="${prefix}.geneBodyCoverage.txt"
    local cmd="${prefix}.cmd.sh"
    local log="${prefix}.rseqc.log"
    if [[ -f ${out_txt} ]] 
    then 
        echo "file exists: ${out_txt}"
    else
        echo "${gb_cmd} -i ${bam} -r ${bed} -o ${prefix} > ${log}"  > ${cmd}
        bash ${cmd}
        # echo "123" >> ${cmd}
    fi
}
export -f gbcov


# for multiple bam files in parallel
# para: <out_dir> <genome> <n_jobs> <bam1> ... <bamn>
function gbcovx() {
    [[ $# -lt 4 ]] && echo "Usage: gbcovx <out_dir> <genome> <n_jobs> <bam1> ... <bamN>" && return 1
    local out_dir=$1
    local genome=$2
    local n_jobs=$3 # number of jobs in parallel
    local args=($@)
    local bam=${args[@]:3}
    # run in parallel
    echo ${bam[@]} | xargs -n 1 | parallel -j ${n_jobs} gbcov ${out_dir} ${genome} {}
}
export -f gbcovx

gbcovx $@
# echo 'genebody_coverage.py --help'

