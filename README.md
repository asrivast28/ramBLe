# Bayesian Networks - Discovering Markov Blankets

This repository implements multiple constraint-based algorithms for discovering Markov Blankets (MB) of a single variable.

## Requirements
* **gcc** (with C++14 support)  
This project has been tested only on Linux platform, using gcc with C++14 support.
* **[Boost](http://boost.org/)**  
Boost libraries are used for parsing the command line options, logging, and a few other purposes.  
* **[SCons](http://scons.org/)**  
SCons, a cross-platform Python based build environment, is required for building the project.  
* The following two repositories are used as submodules:  
  * **[SABNAtk](https://gitlab.com/SCoRe-Group/SABNAtk)**  
  SABNAtk library is used for executing counting queries over the datasets.  
  * **[C++ Utils](https://github.com/asrivast28/cpp-utils)**  
  Logging functionality from the repository is used.
* **[Google Test](https://github.com/google/googletest)** (optional)   
Google Test framework is used for unit testing in this project.  
This dependency can be avoided by not building the unit tests. See the relevant section in [Building](#unit-tests) for more information.   

## Building
After all the dependencies have been installed, the project can be built as:  
<pre><code>scons
</code></pre>  
This will create an executable named `mb`, which can be used for MB discovery.  
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
Debug version of the executable is named `mb_debug`.

#### Logging
By default, logging is disabled in the release build and enabled in the debug build.
In order to change the default behavior, `LOGGING=[0,1]` argument can be passed to `scons`:  
<pre><code>scons LOGGING=1 # Enables logging in the release build
</code></pre>
Please be aware that enabling logging will affect the performance.

## Execution
Once the project has been built, the executable can be used for MB discovery as follows:
<pre><code>./mb -f test/coronary.csv -n 6 -m 1841 -t Smoking ## Expected Output: M. Work,P. Work,Pressure,Proteins,
</code></pre>  
Please execute the following for more information on all the options that the executable accepts:
<pre><code>./mb --help 
</code></pre>
