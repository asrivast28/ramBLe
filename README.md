# Bayesian Networks - Constraint-based Structure Learning
[![](https://github.com/asrivast28/discover-mb/workflows/Build%20and%20Unit%20Tests/badge.svg)](https://github.com/asrivast28/discover-mb/actions)
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

This repository implements multiple constraint-based algorithms for structure learning from data.

## Requirements
* **gcc** (with C++14 support)  
This project has been tested only on Linux platform, using gcc with C++14 support.
* **[Boost](http://boost.org/)**  
Boost libraries are used for parsing the command line options, logging, and a few other purposes.  
* **[SCons](http://scons.org/)**  
SCons, a cross-platform Python based build environment, is required for building the project.  
* The following repositories are used as submodules:
  * **[mxx](https://gitlab.com/patflick/mxx)**  
  mxx is used as a C++ wrapper for MPI.
  * **[C++ Utils](https://github.com/asrivast28/cpp-utils)**  
  Logging functionality from the repository is used.
* **[Google Test](https://github.com/google/googletest)** (optional)   
Google Test framework is used for unit testing in this project.  
If this dependency is not satisfied, then the unit tests are not built. See the relevant section in [Building](#unit-tests) for more information.   

## Building
After all the dependencies have been installed, the project can be built as:  
<pre><code>scons
</code></pre>  
This will create an executable named `csl`, which can be used for constraint-based structure learning.  
Path to the Boost libraries can be specified as:  
<pre><code>scons BOOSTLIBPATH=&lt;Boost library path, default is /usr/lib/x86_64-lib-gnu&gt;
</code></pre>

#### Unit Tests
The unit tests are built by default. The following can be executed for building only the executable:  
<pre><code>scons TEST=0
</code></pre>  

#### Debug
For building the debug version of the executable, the following can be executed:
<pre><code>scons DEBUG=1
</code></pre>  
Debug version of the executable is named `csl_debug`.

#### Logging
By default, logging is disabled in the release build and enabled in the debug build.
In order to change the default behavior, `LOGGING=[0,1]` argument can be passed to `scons`:  
<pre><code>scons LOGGING=1 # Enables logging in the release build
</code></pre>
Please be aware that enabling logging will affect the performance.

## Execution
Once the project has been built, the executable can be used for MB discovery as follows:
<pre><code>./csl -f test/coronary.csv -n 6 -m 1841 -t Smoking -b ## Expected Output: M. Work,P. Work,Pressure,Proteins,
</code></pre>  
Please execute the following for more information on all the options that the executable accepts:
<pre><code>./csl --help
</code></pre>

## Algorithms

### Direct Discovery
This class of algorithms first finds the MB of the variable.
We have added a symmetry correction step to all these algorithms.
* `gs`  
[Grow-Shrink (GS)](https://papers.nips.cc/paper/1685-bayesian-network-induction-via-local-neighborhoods) algorithm by Margaritis & Thrun.
* `iamb`  
[Incremental Association MB (IAMB)](https://www.aaai.org/Library/FLAIRS/2003/flairs03-073.php) algorithm by Tsamardinos et al.
* `inter.iamb`  
[Interleaved Incremental Association MB (InterIAMB)](https://www.aaai.org/Library/FLAIRS/2003/flairs03-073.php) by Tsamardinos et al.

### Topological Discovery
This class of algorithms first finds the parent-child (PC) sets of nodes.
* `mmpc`  
[Max-Min PC (MMPC)](https://link.springer.com/article/10.1007/s10994-006-6889-7) algorithm by Tsamardinos et al. and corrected by Pena et al.
* `hiton`  
[HITON-PC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1480117/) algorithm by Aliferis et al. and corrected by Pena et al.
* `si.hiton.pc`  
[Semi-interleaved HITON-PC](http://www.jmlr.org/papers/v11/aliferis10a.html) algorithm by Aliferis et al.
* `getpc`  
[Get PC](https://www.sciencedirect.com/science/article/pii/S0888613X06000600) algorithm by Pena et al.

## Licensing
Our code is licensed under the Apache License 2.0 (see LICENSE).
The licensing does not apply to the `ext` folder. Please refer to the individual files for their licensing terms.
