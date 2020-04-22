# Experiments

## Data sets
We used the following three gene expression data sets from two model organisms, _Saccharomyces cerevisiae_ and _Arabidopsis thaliana_, for our experiments.
* A data set created by [Tchourine et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5987223/) from multiple RNA-seq studies of _S. cerevisiae_  
2,577 observations for 5,716 genes; can be downloaded from [here](https://zenodo.org/record/3355524#.Xpx0t1NKhhE)
* An unpublished data set created from multiple microarray studies of _A. thaliana_  
16,838 observations for 18,380 genes; will be made available soon
* A manually curated subset of the above data set focusing only on the studies of the development process in _A. thaliana_
5,102 observations for 18,373 genes; will be made available soon

### Discretization
We discretize the expression levels in the data sets using the methodology suggested by [Fridman et al.](https://www.ncbi.nlm.nih.gov/pubmed/11108481)  
For example, the _S. cerevisiae_ data set can be discretized (using [`discretize.py`](scripts/discretize.py)) as follows:
<pre><code>scripts/discretize.py -f yeast_microarray_expression.tsv -s '\t' -c -v -i -o yeast_microarray_expression_discretized.tsv
</code></pre>

## Measuring Performance

### Building
The simplest way to build the code for measuring performance is to execute the following:
<pre><code>scons TIMER=1
</code></pre>
More information on building can be found in [`README.md`](README.md#building)

### Run-time Environment
The run-time environment for the experiments (collected using the file [`collect_environment.sh`](https://github.com/SC-Tech-Program/Author-Kit/blob/master/collect_environment.sh))
is available at [`hive_environment.log`](scripts/hive_environment.log)

### Sequential Execution
The code can be used for learning the Bayesian network, using any of the [supported algorithms](README.md#algorithms), as described in [`README.md`](README.md#execution)  
For example, in order to measure the performance of our network for learning the network from the _S. cerevisiae_ data set using the _GS_ algorithm, the following can be executed:
<pre><code>./ramble -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

#### Running _bnlearn_
We have also provided an [Rscript](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/Rscript), [`ramble_bnlearn.R`](scripts/ramble_bnlearn.R), for running [_bnlearn_](https://www.bnlearn.com/) with the same arguments as our executable.  
For example, the performance of _bnlearn_ in learning the network from the _S. cerevisiae_ data set using the _GS_ algorithm can be measured by executing:
<pre><code>scripts/ramble_bnlearn.R -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

### Parallel Execution
The performance of our executable when run in parallel using MPI can be measured as follows:
<pre><code>mpirun -np 16 ./ramble -n 5716 -m 2577 -f yeast_microarray_expression_discretized.tsv -s '\t' -c -v -i -a gs -o yeast_network.dot -d
</code></pre>

### Scalability Experiments
We have provided a utility script, [`run_experiments.py`](scripts/run_experiments.py), for running the scalability experiments using our executable.  
By default, this script runs the scalability experiments for all three data sets while varying the number of processors between 1 and 1024 and records the measured times in a CSV file.
The runs can be customized using different arguments to the script, which can be seen by executing:
<pre><code>scripts/run_experiments.py -h
</code</pre>
