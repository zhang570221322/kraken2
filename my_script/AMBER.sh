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
rm weight.output.AMBER
rm origin.output.AMBER
rm binning.AMBER
cat *.weight.output.AMBER >weight.output.AMBER
cat *.origin.output.AMBER >origin.output.AMBER
cat *.binning.AMBER >binning.AMBER
amber=/data/home/wlzhang/classfication/software/statistics/AMBER/amber.py
taxonomy_db=$kraken2_db_cus/taxonomy
python $amber -g $file_dir/binning.AMBER -l "Origin,Weight " $file_dir/origin.output.AMBER $file_dir/weight.output.AMBER --ncbi_nodes_file $taxonomy_db/nodes.dmp --ncbi_names_file $taxonomy_db/names.dmp --ncbi_merged_file $taxonomy_db/merged.dmp -o output_dir/
