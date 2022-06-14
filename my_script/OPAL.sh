#!/bin/bash
#PBS -N kranken test
#PBS -l nodes=node7:ppn=48
#PBS -l walltime=9000:00:00
##PBS -V
set +e
# 加载环境
source /data/home/wlzhang/Work_Environment.sh
python=/data/home/wlzhang/anaconda3/envs/bio/bin/python
# var
file_dir=/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/normal/data
amber_format=/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/normal/script/AMBER/amber_format.py

function run_once() {
    file=$1
    example_id=$(basename $file .fastq)
    echo "AMBER格式"
    $python $amber_format $file 0.0.1 $example_id
    # rm $statistics_name
}
cd ${file_dir}
for file in $(ls $file_dir/*.fastq); do
    run_once $file
done

echo "合并数据"
cd ${file_dir}
# rm weight.profiling.OPAL
# rm origin.profiling.OPAL
# rm binning.OPAL
# cat *.weight.profiling.OPAL >weight.profiling.OPAL
# cat *.origin.profiling.OPAL >origin.profiling.OPAL
# cat *.binning.OPAL >binning.OPAL
# opal=/data/home/wlzhang/classfication/software/statistics/OPAL/opal.py
python $opal -g binning.OPAL origin.profiling.OPAL weight.profiling.OPAL GCN.origin.OPAL GCN.weight.OPAL -l "Origin,Weight,GCN,Weight+GCN" -o output_dir/
