#!/bin/bash
#PBS -N kranken test
#PBS -l nodes=node7:ppn=48
#PBS -l walltime=9000:00:00
##PBS -V
set +e
 
opal=/data/home/wlzhang/classfication/software/statistics/OPAL/opal.py
cd /data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/normal/script/OPAL/OPAL_DATA
python $opal  --plot_abundances  --metrics_plot_re  c,p,l,w,f,t  -g  binning.OPAL  origin.profiling.OPAL cemtrofuge.OPAL  k2mem.OPAL mash.OPAL sourmash.OPAL motu.OPAL  GCN.weight.OPAL   -l  Kraken2,Centrifuge,k2Mem,Mash,SourMash,mOTU,Weight+GCN  -o output_dir/

