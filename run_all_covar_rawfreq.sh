#!/bin/bash

# WARNING: DO NOT RUN THIS BENCHMARK AND THE RAWFREQ BENCHMARK CONCURRENTLY.
# Can't parallelise script (currently) because of dependence on symlink to relevant arch. Boo
# Also TODO: get the python program to read in the inputs for each test protein on the fly (for a given architecture); that would avoid the time-consuming step of Theano checking whether it should recompile functions.

# UPDATE 21/01/2018 to include min and max sequence separation
# UPDATE 22/01/2018 additional param governs whether to process covar or rawfreq data

if [ -z "$3" ]; then
    echo "Usage: $0 'covar'|'rawfreq' min_seq_separation max_seq_separation" >&2
    exit 1
fi

ip_type=$1
topomin=$2
topomax=$3

if [ $ip_type != 'covar' ] && [ $ip_type != 'rawfreq' ]; then
    echo "argument 1 must be one of 'covar' or 'rawfreq'." >&2
    exit 1
fi

run_psicov150_bench=1

root_dir=${PWD}

# you may need to edit this path:
#psicov_150_dir=${PWD}/psicov150
psicov_150_dir=${HOME}/projects/DeepCov/psicov150

tgt_list=${psicov_150_dir}/target.lst.sorted

# directory with different CNN architecture definitions and params
if [ $ip_type == 'covar' ]; then
    arch_dir=${root_dir}/covar-model-all
    arch_subdir_prefix='covar'
else
    arch_dir=${root_dir}/pairfreq-model-all
    arch_subdir_prefix='pairfreq'
fi

bench_script=${root_dir}/run_psicov150_bench_singleRF.sh

# filenames; these are hardcoded in the NN predictor python script and are constant between covar and rawfreq models
nndef='nndef.py'
params='FINAL_fullmap_metapsicov_model.npz'

master_out_file_mean=${root_dir}/all_windowsize_results_MEAN_${ip_type}_min${topomin}_max${topomax}.txt
#master_out_file_median=${root_dir}/all_windowsize_results_MEDIAN_${ip_type}_min${topomin}_max${topomax}.txt

printf "window\tL\tL/2\tL/5\tL/10\tL=100\n" > $master_out_file_mean
#printf "window\tL\tL/2\tL/5\tL/10\tL=100\n" > $master_out_file_median

# rf = receptive field (window size)
for rf in $(cat ${arch_dir}/field_sizes.txt)
do
    echo "Processing RF size $rf ..."

    curr_arch_dir=${arch_dir}/${arch_subdir_prefix}${rf}

    out_dir=${curr_arch_dir}/bench_results_aggregated_min${topomin}_max${topomax}
    mkdir -p $out_dir
    
    if [ $run_psicov150_bench != 0 ]; then

	# create symlink to nndef.py and .npz file for current arch, removing any file in cwd with the same names
	ln -sfv ${curr_arch_dir}/$nndef $nndef
	ln -sfv ${curr_arch_dir}/$params $params
	
	echo "Running benchmark..."
	
	$bench_script $ip_type $topomin $topomax $psicov_150_dir >/dev/null 2>${out_dir}/bench1.err
    fi

    echo "Running postprocess scripts..."

    ${root_dir}/extract_results_from_benchmark.sh ${curr_arch_dir}/results_min${topomin}_max${topomax} $out_dir $topomin $topomax ${tgt_list} > /dev/null 2> ${out_dir}/bench2.err

    cd $out_dir

    printf "$rf\t" >> $master_out_file_mean
    # printf "$rf\t" >> $master_out_file_median

    Rscript ${root_dir}/bench_plots_and_stats_covar_pairfreq.R $ip_type $topomin $topomax ${root_dir}/src $root_dir 2>${out_dir}/bench3.err

    cd - > /dev/null
    echo "Done."

    # cleanup symlinks and one .pyc file
    rm -f $nndef $params ${nndef}c
done

echo "All done."
