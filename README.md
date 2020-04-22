# ramBLe - A Parallel Framework for Bayesian Learning
[![](https://github.com/asrivast28/discover-mb/workflows/Build%20and%20Unit%20Tests/badge.svg)](https://github.com/asrivast28/discover-mb/actions)
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

ramBLe (A Parallel F**ram**ework for **B**ayesian **Le**arning) supports multiple constraint-based algorithms for structure learning from data in parallel.

## Requirements
* **gcc** (with C++14 support) is used for compiling the project.  
_This project has been tested only on Linux platform, using version [9.2.0](https://gcc.gnu.org/gcc-9/changes.html)._
* **[Boost](http://boost.org/)** libraries are used for parsing the command line options, logging, and a few other purposes.  
_Tested with version [1.70.0](https://www.boost.org/users/history/version_1_70_0.html)._
* **[MPI](https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/mpi31-report.htm)** is used for execution in parallel.  
_Tested with [MVAPICH2 version 2.3.3](http://mvapich.cse.ohio-state.edu/static/media/mvapich/mvapich2-2.3.3-userguide.html)._
* **[SCons](http://scons.org/)** is required for building the project.  
_Tested with version [3.1.2](https://scons.org/doc/3.1.2/HTML/scons-user.html)._
* The following repositories are used as submodules:
  * **[mxx](https://gitlab.com/patflick/mxx)** is used as a C++ wrapper for MPI.  
  * **[Graph API](https://github.com/asrivast28/cpp-utils)** is used as a lightweight wrapper around Boost.Graph.  
  * **[C++ Utils](https://github.com/asrivast28/cpp-utils)** are used for logging and timing.  
* **[Google Test](https://github.com/google/googletest)** (optional) framework is used for unit testing in this project.   
If this dependency is not satisfied, then the unit tests are not built. See the relevant section in [Building](#unit-tests) for more information.  
_Tested with version [1.10.0](https://github.com/google/googletest/releases/tag/release-1.10.0)._

## Building
After the dependencies have been installed, the project can be built as:  
<pre><code>scons
</code></pre>  
This will create an executable named `ramble`, which can be used for constraint-based structure learning.  
By default, all the paths from the environment in `CPATH` and `LIBRARY_PATH` variables are used as include paths and library paths.  
Path to external includes and libraries at non-default locations can also be specified as:  
<pre><code>scons LOCALINCLUDES=&lt;comma-delimited list of paths&gt; LOCALLIBS=&lt;comma-delimited list of paths&gt;
</code></pre>

#### Unit Tests
The unit tests are built by default. The following can be executed for building only the executable:  
<pre><code>scons TEST=0
</code></pre>  

#### Debug
For building the debug version of the executable, the following can be executed:
<pre><code>scons DEBUG=1
</code></pre>  
Debug version of the executable is named `ramble_debug`.

#### Logging
By default, logging is disabled in the release build and enabled in the debug build.
In order to change the default behavior, `LOGGING=[0,1]` argument can be passed to `scons`:  
<pre><code>scons LOGGING=1 # Enables logging in the release build
</code></pre>
Please be aware that enabling logging will affect the performance.

#### Timing
Timing of high-level operations can be enabled by passing `TIMER=1` argument to `scons`.

## Execution
Once the project has been built, the executable can be used for learning BN as follows:
<pre><code>./ramble -f test/coronary.csv -n 6 -m 1841 -d -o test/coronary.dot
</code></pre>  
For running in parallel, the following can be executed:
<pre><code> mpirun -np 8 ./ramble -f test/coronary.csv -n 6 -m 1841 -d -o test/coronary.dot
</code></pre>  
Please execute the following for more information on all the options that the executable accepts:
<pre><code>./ramble --help
</code></pre>

## Algorithms
The algorithm for learning BNs can be chosen by specifying the desired algorithm as an option to the executable, using `-a` option. The currently supported algorithms are listed below.

### Blanket Learning
This class of algorithms first finds the Markov blanket (MB) of the variable to get the parents and the children (PC).
* `gs` corresponds to the [Grow-Shrink (GS)](https://papers.nips.cc/paper/1685-bayesian-network-induction-via-local-neighborhoods) algorithm by Margaritis & Thrun.
* `iamb` corresponds to the [Incremental Association MB (IAMB)](https://www.aaai.org/Library/FLAIRS/2003/flairs03-073.php) algorithm by Tsamardinos et al.
* `inter.iamb` corresponds to the [Interleaved Incremental Association MB (InterIAMB)](https://www.aaai.org/Library/FLAIRS/2003/flairs03-073.php) by Tsamardinos et al.

## Publication
_Currently under double-blind review._

## Licensing
Our code is licensed under the Apache License 2.0 (see [`LICENSE`](LICENSE)).
The licensing does not apply to the `ext` folder. Please refer to the individual files in the folder for their licensing terms.
