# DeepCov v1.0
### Fully convolutional neural networks for protein residue-residue contact prediction

David T. Jones and Shaun M. Kandathil

University College London

Requirements:
-------------
- Bash shell
- Working C and C++ compilers (tested with GCC 4.8.5)
- Python 2 (tested on 2.7.5) or 3 (tested on 3.4.5) with development libraries and headers

- The following Python modules (version numbers in brackets were used during development/testing):
  - numpy (1.13.1)
  - Theano (0.9.0)
  - Lasagne (0.2.dev1)

At the time of writing, pip will install Lasagne 0.1 by default, which will not work due to changes in Theano 0.9. You may need to use the 'bleeding-edge' install of Lasagne:

`$ pip install --upgrade https://github.com/Lasagne/Lasagne/archive/master.zip`

On some distributions, the C++ compiler is a separate add-on package and may not be installed by default. For example, on CentOS you will need to `yum install` packages `gcc` AND `gcc-c++`.

To get the Python development headers and libs you may need to install a separate package, the name of which will depend on your package manager. For example, on CentOS this is `python-devel` or `python34-devel`.

Setup and testing:
------------------
Run `setup.sh`.

This will compile and test a C executable, cov21stats. This executable generates covariance or pair frequency data from your input alignment.
The script will also test the DeepCov prediction pipeline on a test input alignment, so make sure all dependencies listed above are in place before running it. By default, the scripts will use whatever `cc` and `python3` point to in your shell. These can be changed in `deepcov.sh` and `setup.sh`.

The testing procedure will compare a newly generated contact prediction (in `test/`) against the reference file found in `test/example_io`. Since different OS/compiler combinations can lead to very slightly different contact scores, only the ranking of the contacts is evaluated when deciding whether the test was successful. To see if there are any differences, please compare the two contact files using a program such as `sdiff`.

Running:
--------
`$ /path/to/deepcov.sh [-h] [-m model_type] [-r receptive_field] -i input_file [-o output_contact_file]`

The optional arguments -m and -r are primarily a means to reproduce results in our paper. For most 'production' purposes, you can leave these set to their defaults (covariance model + receptive field of 41 residues).

The input alignment must be in the PSICOV format. If your alignment is in a different format, we recommend using the [ConKit Python module](https://pypi.python.org/pypi/conkit) to reformat it.

The output is in the CASP contact format.

An example input alignment is provided at `test/example_io/1guuA.aln`. The corresponding DeepCov output contact file is `test/example_io/1guuA.con`.

Tips:
-----
For inferring contacts for single alignments, we find that running DeepCov on a (reasonably recent) CPU is faster than running on a GPU, when considering end-to-end runtime on our benchmark sets. For this reason, DeepCov will run on your CPU by default. If you'd like to change this behaviour, edit `deepcov.sh` and change the value of the `THEANO_FLAGS` variable near the end of the script (see http://deeplearning.net/software/theano/library/config.html for more details on this and other variables). You will also need to install other prerequisites for running on the GPU; please refer to Theano's documentation.

Benchmarking scripts:
------
We've included some additional scripts that should reproduce results from our paper.
For running the benchmarking scripts, you will need a recent install of R in addition to the dependencies listed above.
You will also need the PSICOV150 test set, which comes with its own README and can be downloaded [here.](http://bioinfadmin.cs.ucl.ac.uk/downloads/contact_pred_datasets/)

Once the dataset is in place, edit `run_all_covar_rawfreq.sh` to specify the location of the psicov150 set, and then run it, e.g.

`./run_all_covar_rawfreq.sh covar 6 11`

where 6 and 11 refer to the min and max sequence separation you want to consider, and 'covar' refers to the covariance model.

With these inputs, output will be generated in your DeepCov installation directory, in a file named `all_windowsize_results_MEAN_covar_min6_max11.txt`.

PLEASE NOTE: the benchmarking process does create a number of rather large files. Use with caution if you have limited storage.

Training scripts:
-------
An example training script and a README can be found in `training/`, which includes a link to where training data can be found. 

Citing:
-------
If you find DeepCov useful, please cite our paper in Bioinformatics:

Jones DT and Kandathil SM (2018). High precision in protein contact prediction using fully convolutional neural networks and minimal sequence features. _Bioinformatics_ **34**(19): 3308-3315. [Link](https://doi.org/10.1093/bioinformatics/bty341)
