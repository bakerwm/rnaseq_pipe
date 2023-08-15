
## submit
[[ $# -lt 2 ]] &&  echo "Usage: submit_jobs.sh <wk_dir> [1|0] # default [0]" && exit 1
wk_dir=$1
run=$2

## submit large jobs to slurm
# wk_dir="/data1/yuyang/wangm/work/2023/rnaseq_xp/20230630_B6xB6_xp"
data_dir=${wk_dir}/data/raw_data
cd ${wk_dir}

function make_slurm() {
    local fq1=$1
    local wk_dir=$2
    [[ $# -lt 2 ]] && echo "missing args" && return 1
    [[ ! ${fq1} = *1.fq.gz ]] && echo "not fq1: ${fq1}" && return 1
    local fq2=${fq1/1.fq/2.fq}
    local fname=$(basename ${fq1/_1.fq.gz})
    local job_dir=${wk_dir}/jobs
    [[ ! -d ${job_dir} ]] && mkdir -p ${job_dir}
    local job_file="${job_dir}/submit_${fname}_slurm.sh"
    # sbatch run_slurm.sh ${fq1}
    bash make_slurm.sh ${wk_dir} ${fq1} ${fq2} > ${job_file}
    echo ${job_file}    
}


for fq1 in ${data_dir}/*1.fq.gz
do
    ## create job ##
    job_file=$(make_slurm ${fq1} ${wk_dir})
    ## status ##
    [[ $(grep -ci SBATCH ${job_file}) -gt 4 ]] && flag="ok" || flag="failed"
    [[ ${run} == 1 ]] && ss="yes" || ss="no"
    printf "%8s: %-4s %8s: %-4s : %s\n" "file" ${flag} "submit" ${ss} $(basename ${job_file})
    ## submit ##
    [[ ${run} == 1 && ${flag} == "ok" ]] && sbatch ${job_file}
done

