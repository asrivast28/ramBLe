# Experiments
This file is intended to serve as a guide for reproducing the results presented in our paper.

## Contents
* [Run-time Environment](#run-time-environment)
  * [Hardware](#hardware)
  * [Software](#software) 
* [Initial Setup](#initial-setup)
  * [Cloning](#cloning)
  * [Building](#building)
  * [Data sets](#data-sets)
    * [Discretization](#discretization)
* [Measuring Performance](#measuring-performance)
  * [Sequential Execution](#sequential-execution)
    * [Running _bnlearn_](#running-bnlearn) 
  * [Parallel Execution](#parallel-execution)
  * [Scalability Experiments](#scalability-experiments)


## Run-time Environment
We provide a summary of the hardware and the software required for running the experiments. The run-time environment for the experiments (collected using the file [`collect_environment.sh`](https://github.com/SC-Tech-Program/Author-Kit/blob/master/collect_environment.sh))
is available at [`hive_environment.log`](hive_environment.log)

### Hardware
We used the [Hive cluster at Georgia Tech](https://docs.pace.gatech.edu/hive/gettingStarted/) for our experiments. Each node in the cluster has a 2.7 GHz 24-core Intel Xeon 6226 processor and main memory of 192 GB or more. The nodes are connected via EDR (100 Gbps) InfiniBand. More details on the cluster resources can be found in the [cluster documentation](https://docs.pace.gatech.edu/hive/resources/).

### Software
We conducted our experiments on nodes running Linux [**RHEL** _v7.6_](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/7.6_release_notes/index) operating system. We used the following versions of the compiler and other libraries for experimenting with _ramBLe_. 
* [**gcc** _v9.2.0_](https://gcc.gnu.org/gcc-9/changes.html)
* [**MVAPICH** _v2.3.3_](http://mvapich.cse.ohio-state.edu/static/media/mvapich/mvapich2-2.3.3-userguide.html)
* [**Boost** _v1.70.0_](https://www.boost.org/users/history/version_1_70_0.html)
* [**SCons** _v3.1.2_](https://scons.org/doc/3.1.2/HTML/scons-user.html)
* [**Google Test** v1.10.0](https://github.com/google/googletest/releases/tag/release-1.10.0)  
  _Used for unit testing. Not required if [building without unit tests](https://github.com/asrivast28/ramBLe#unit-tests)._

The purpose of all the libraries is explained in more detail in [`README.md`](README.md#requirements).

## Initial Setup
### Cloning
_ramBLe_ can be downloaded by cloning this Github repo, along with all the submodules, by executing the following for cloning via SSH:
<pre><code>git clone --recurse-submodules git@github.com:asrivast28/ramBLe.git
</code></pre>
or the following for cloning via HTTPS:
<pre><code>git clone --recurse-submodules https://github.com/asrivast28/ramBLe.git
</code></pre>
Similarly, the tree at a particular revision can be checked out along with the corresponding version of all the submodules by executing the following:
<pre><code>git checkout --recurse-submodules &lt;tree-ish&gt;
</code></pre>

### Building
The simplest way to build the code for measuring performance is to execute the following:
<pre><code>scons TIMER=1
</code></pre>
More information on building can be found in [`README.md`](README.md#building)

### Data sets
We used the following three gene expression data sets from two model organisms, _Saccharomyces cerevisiae_ and _Arabidopsis thaliana_, for our experiments.
* A data set created by [Tchourine et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5987223/) from multiple RNA-seq studies of _S. cerevisiae_  
2,577 observations for 5,716 genes; can be downloaded from [Zenodo](https://zenodo.org/record/3355524#.Xpx0t1NKhhE)
* An unpublished data set created from multiple microarray studies of _A. thaliana_  
16,838 observations for 18,380 genes; will be made available soon
* A manually curated subset of the above data set focusing only on the studies of the development process in _A. thaliana_
5,102 observations for 18,373 genes; can also be downloaded from [Zenodo](https://zenodo.org/record/4672797#.YG9TQhNKhQI)

#### Discretization
We discretize the expression levels in the data sets using the methodology suggested by [Fridman et al.](https://www.ncbi.nlm.nih.gov/pubmed/11108481)  
For example, the _S. cerevisiae_ data set can be discretized (using [`discretize.py`](https://github.com/asrivast28/bn-utils/blob/31f4957cbbf4c6cf0451b4139f3b54f9bd4cee90/scripts/discretize.py)) as follows:
<pre><code>common/scripts/discretize.py -f yeast_microarray_expression.tsv -s '\t' -c -v -i -o yeast_microarray_expression_discretized.tsv
</code></pre>

## Measuring Performance
### Sequential Execution
_ramBLe_ can be used for learning a Bayesian network, using any of the [supported algorithms](README.md#algorithms), as described in [`README.md`](README.md#execution)  
For example, in order to measure the performance of our network for learning the network from the _S. cerevisiae_ data set using the _GS_ algorithm, the following can be executed:
<pre><code>./ramble -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

#### Running _bnlearn_
We have also provided an [Rscript](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/Rscript), [`ramble_bnlearn.R`](https://github.com/asrivast28/bn-utils/blob/31f4957cbbf4c6cf0451b4139f3b54f9bd4cee90/scripts/ramble_bnlearn.R), for running [_bnlearn_](https://www.bnlearn.com/) with the same arguments as _ramBLe_.  
For example, the performance of _bnlearn_ in learning the network from the _S. cerevisiae_ data set using the _GS_ algorithm can be measured by executing:
<pre><code>common/scripts/ramble_bnlearn.R -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

### Parallel Execution
The performance of _ramBLe_ when run in parallel using MPI can be measured as follows:
<pre><code>mpirun -np 16 ./ramble -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

### Scalability Experiments
We have provided a utility script, [`ramble_experiments.py`](https://github.com/asrivast28/bn-utils/blob/31f4957cbbf4c6cf0451b4139f3b54f9bd4cee90/scripts/ramble_experiments.py), for running the scalability experiments using _ramBLe_.  
By default, this script runs the scalability experiments for all three data sets while varying the number of processors between 1 and 1024 and records the measured times in a CSV file.
The runs can be customized using different arguments to the script, which can be seen by executing:
<pre><code>common/scripts/ramble_experiments.py -h</code</pre>
