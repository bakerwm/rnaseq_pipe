
## global variables
## in-case, define the variables by config.txt file
SRC_DIR=$(dirname $(realpath -s $0)) # script dir
[[ -z ${PARTITION_NAME} ]] && PARTITION_NAME="CPU2"
[[ -z ${N_CPU} ]] && N_CPU=8


## make slurm script
function make_slurm() {
    local wk_dir=$1
    local fq1=$2
    local fq2=$3
    [[ $# -lt 3 ]] && echo "not enough args..." && return 1
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    wk_dir=$(realpath -s ${wk_dir})
    fq1=$(realpath -s ${fq1})
    fq2=$(realpath -s ${fq2})
    ## check run_pipe.sh
    run_rnaseq="${SRC_DIR}/run_rnaseq.sh"
    [[ ! -f ${run_rnaseq} ]] && echo "run_rnaseq.sh not exists" && return 1

    cat << EOF
#!/bin/bash
#SBATCH --job-name=rna           # create a short name for your job
#SBATCH --partition=${PARTITION_NAME}         # node name
#SBATCH --ntasks=${N_CPU}               # total cores, 
#SBATCH --mem=30000              # 30GB for STAR
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=5-01:00:00        # total run time limit (D-HH:MM:SS)
#SBATCH --output=log.${PARTITION_NAME}.%j.out

module purge
# source ${HOME}/miniconda3/bin/activate
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate hiseq

cd ${wk_dir}
bash ${run_rnaseq} ${wk_dir} ${fq1} ${fq2}
EOF
}
export -f make_slurm 


[[ $# -lt 3 ]] && echo "Usage: make_slurm.sh <wk_dir> <fq1> <fq2>" && exit 1
make_slurm $1 $2 $3
