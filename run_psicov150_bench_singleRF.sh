#!/bin/bash
# EDIT 21/01/2018 to include min and max sequence separation values

if [ -z "$4" ]; then 
    echo "Usage: $0 <'covar' or 'rawfreq'> <min seq separation> <max seq separation> <path to psicov150>" >&2
    echo "Please use run_all_covar_rawfreq.sh to run the benchmark. Do not use this script directly." >&2
    exit 1
fi

ip_type=$1
topomin=$2
topomax=$3
# dir containing aln/ and pdb/ dirs with benchmark data
aln_pdb_root=$4

if [ $ip_type != 'covar' ] && [ $ip_type != 'rawfreq' ]; then
    echo "argument 1 must be one of 'covar' or 'rawfreq'." >&2
    exit 1
fi

aln_path=${aln_pdb_root}/aln
pdb_path=${aln_pdb_root}/pdb

if [ ! -s nndef.py ]; then
    echo "$0 : ERROR: no nndef.py found in current working dir. Exiting." >&2
    exit 2
fi

resolved_path_to_nndef=$(readlink -f nndef.py)
out_root=$(dirname $resolved_path_to_nndef)

#these may not exist yet

con_path=${out_root}/con
results_path=${out_root}/results_min${topomin}_max${topomax}

# paths to executables/scripts
cov21stats=${PWD}/bin/cov21stats
contactbench=${PWD}/contactbench/contactbench

# because sometimes the system default python isn't the python you want:
python=python


if [ $ip_type == 'rawfreq' ]; then
    c21_options='-p'
    path_21stats=${aln_pdb_root}/21m
else
    c21_options=''
    path_21stats=${aln_pdb_root}/21c
fi

# this needs to be given 21stats file as argv[1]:
nnpredictor=${PWD}/predictor.py
 
mkdir -p $path_21stats $con_path $results_path

export THEANO_FLAGS='device=cpu,floatX=float32'

for n in $(cat ${aln_pdb_root}/target.lst.sorted)
do
    echo "$0 : Processing $n ..."

    aln_file=${aln_path}/${n}.aln
    pdb_file=${pdb_path}/${n}.pdb

    if [ ! -s $aln_file ]; then 
	echo "$0 : WARNING: alignment file $aln_file not found or is empty. Skipping this protein." >&2
	continue
    fi

    if [ ! -s $pdb_file ]; then 
	echo "$0 : WARNING: pdb file $pdb_file not found or is empty. Skipping this protein." >&2
	continue
    fi

    file_21stats=${path_21stats}/${n}.21stats

    if [ ! -s $file_21stats ]; then 
	$cov21stats $c21_options $aln_file $file_21stats
    else
	echo "$0 : Existing 21stats file found. Proceeding to contact prediction..."
    fi

    # predict contacts from the 21stats file. Once again, use saved results if existent.

    con_file=${con_path}/${n}.con

    if [ ! -s $con_file ]; then
	echo "$0 : Preparing contact file..."
	$python $nnpredictor $file_21stats > $con_file
    else
	echo "$0 : Existing contact file found. Proceeding to bench..."
    fi

    # run contactbench on the predictions and the pdb file

    results_file=${results_path}/${n}.results

    echo "$0 : Running benchmark..."
    $contactbench $pdb_file $con_file $topomin $topomax > $results_file 
    
    echo 'Done!'
done

unset THEANO_FLAGS
