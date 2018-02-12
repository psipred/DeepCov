#!/bin/bash

###############################################
# Helper script to run DeepCov
# v1.0 by David T. Jones and Shaun M. Kandathil
###############################################

# usage message

help_text="""
usage: $0 [-h] [-m model_type] [-r receptive_field] -i input_file [-o output_contact_file]

DeepCov v1.0 by David T. Jones and Shaun M. Kandathil

required arguments:
  -i input_file, --input-file input_file
                        Path to input alignment in PSICOV format.
optional arguments:
  -h, --help            show this help message and exit
  -m model_type, --model model_type     
                        Model type. Should be one of: 'pairfreq', 'covar', or
                        'covshallow'. Defaults to 'covar'.
  -r receptive_field, --receptive-field receptive_field
                        Receptive field size. Defaults to 41. See DeepCov paper
                        for valid values. Ignored if --model-type is 'covshallow'.
  -o output_contact file
                        Filename for output contacts. The default is to create a
                        file in PWD with '.con' appended to the input filename.
"""

if [ -z "$*" ]; then
    echo "$help_text"
    exit 0
fi

python='python3' ### Redefine this var if your python executable isn't called python3

# Do not change these 3:
params='FINAL_fullmap_metapsicov_model.npz'
nndef='nndef.py'
predictor='predictor.py'

deepcov_root=$(dirname $0)

bin=${deepcov_root}/bin
cov21stats=${bin}/cov21stats

# check that the cov21stats executable exists.
if [ ! -s $cov21stats ] || [ ! -x $cov21stats ]; then
    echo "$0 : Unable to find cov21stats executable in ${bin}. Did you run setup.sh?" >&2
    exit 6
fi

# Parameter defaults
rf=41
model_type='covar'
input_file=''
output_file=''

####### Parse commandline options
# TODO: add an option to provide a list of input filenames, and iterate over them with a given model
while [[ $# -ge 1 ]]
do
    key="$1"
    
    case $key in
	-h|--help)
	    echo "$help_text"
	    exit 0
	    ;;
	-r|--receptive-field)
	    rf="$2"
	    shift # extra shift since this arg has a value
	    ;;
	-m|--model)
	    model_type="$2" # this is a string so we won't check it here
	    shift
	    ;;
	-i|--input-file)
	    input_file="$2"
	    shift
	    ;;
	-o|--output-file)
	    output_file="$2"
	    shift
	    ;;
	*)
            # process any other argument not specified above.
	    # cannot recognise when such an argument also has a value e.g. -f <filename>
	    # but we'll treat that value as another (hopefully unrecognised) argument.
	    echo "$0 : Warning: unrecognised option $key will be ignored." >&2
	    ;;
    esac
    shift # this is the "default" shift
done

####### Check inputs

# check input file
if [ -z $input_file ]; then
    echo "$0 :Error: an input filename must be specified with option -i." >&2
    exit 1
elif [ ! -s $input_file ]; then
    echo "$0 : Error: input alignment file $input_file is either empty or does not exist. Exiting." >&2
    exit 3
fi

# check that rf is a (positive) nonzero integer and not empty
if ! [[ $rf =~ ^[0-9]+$ ]]; then 
	echo "$0 : Error: receptive field value must be a positive integer. You supplied '$rf'. Exiting." >&2
	exit 2
fi

aln_bn=$(basename $input_file)
con_file=''
if [ -z $output_file ]; then
    con_file=${aln_bn}.con
else
    con_file=$output_file
fi

filename_21stats=${aln_bn}.21stats
cov21stats_options=''

# echo "$rf $model_type $input_file $con_file"

case $model_type in
    covshallow|covar)
	if [ $model_type == 'covshallow' ]; then 
	    rf=41
	    model_dir=${deepcov_root}/covar-model-shallow
	else
	    model_dir=${deepcov_root}/covar-model-all/covar${rf}
	fi
	;;
    pairfreq)
	model_dir=${deepcov_root}/pairfreq-model-all/pairfreq${rf}
	cov21stats_options='-p'
	;;
    *)
	echo "$0 : Error: unrecognised model type: '$model_type'. Please see usage by running $0 -h." >&2
	exit 4
	;;
esac

model_file=${model_dir}/$nndef
param_file=${model_dir}/$params

# Determine if user-specified receptive_field is valid before going any further
# For valid model_type (checked above), this check should fail iff there isn't a model with the requested RF size.
if [ ! -s $model_file ] || [ ! -s $param_file ]; then
    echo "$0 : Error: No model available for specified receptive field size ${rf}. Please see DeepCov documentation/paper for valid values." >&2
    exit 5
fi

######### Inputs are (likely) okay; proceed with running

# symlink nndef and params for correct pipeline
# create symlink to nndef.py and .npz file for current arch, removing any file in cwd with the same names
# nndef needs to be in $deepcov_root, and params needs to be in $PWD.
nndef_linkname=${deepcov_root}/$nndef

ln -sf ${model_file} $nndef_linkname
ln -sf ${param_file} $params

# get correct version of input tensor from cov21stats. TODO: If cov21stats fails, stop
echo 'Creating input features...'
$cov21stats $cov21stats_options $input_file $filename_21stats

# Actually run DeepCov
echo 'Running DeepCov...'

# For inference we'll run on the CPU as we've found it to be faster
export THEANO_FLAGS='device=cpu,floatX=float32' #, blas.ldflags="-L/lib64 -lopenblas"'
$python $deepcov_root/$predictor $filename_21stats > $con_file
unset THEANO_FLAGS

echo 'Cleaning up...'
# Remove (frequently large) 21stats file
rm $filename_21stats

# cleanup symlinks and one .pyc file. TODO: is there any benefit in keeping the pycache dir?
rm -rf $nndef_linkname $params ${nndef_linkname}c "${deepcov_root}/__pycache__/"

echo "Results in $con_file"
