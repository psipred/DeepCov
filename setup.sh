#!/bin/bash

# Setup DeepCov dependencies

mkdir -p bin

# compile cov21stats. Using only -O1 to *try* and keep things consistent across platforms.
# -O2 and -O3 are okay to use.
echo '***Compiling cov21stats...'
cc -O2 -m64 -o bin/cov21stats src/cov21stats.c -lm

# compile contactbench
echo '***Compiling contactbench...'
cc -O2 -m64 -o contactbench/contactbench contactbench/contactbench.c -lm

echo '***Testing DeepCov...'

## Test the pipeline using deepcov.sh

ref_con='test/example_io/1guuA.con'
new_con='test/1guuA_test.con'

./deepcov.sh -i test/example_io/1guuA.aln -o $new_con

echo '***Checking results...'

# diff -q $ref_con $new_con
# We can't simply `diff` the new and reference prediction due to different FP behaviour on different systems, which produce slightly different answers.
# Use compare_contacts.py to compare the ordering of the predictions instead, since the ordering of the contacts is what matters most (for measures such as top-L/x precision).

nl_ref=$(wc -l $ref_con | awk '{print $1}')
nl_test=$(wc -l $new_con | awk '{print $1}')

if [ $nl_ref == $nl_test ]; then

    python test/compare_contacts.py $ref_con $new_con

    exit_stat=$?
    
    if [ $exit_stat == 1 ]; then
	echo '***DeepCov test was unsuccessful.' >&2
    elif [ $exit_stat == 2 ]; then
	echo '***Contact file comparison script reported exit status 2. Contact author.' >&2
    else
	echo '***Test passed; DeepCov is ready to use.'
	#    rm $new_con # remove the test prediction if the file is OK
    fi
else
    echo "Errors in output; not proceeding to test." >&2
fi
