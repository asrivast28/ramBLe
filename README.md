# Bayesian Networks - Discovering Markov Blankets

This repository implements multiple constraint-based algorithms for discovering Markov Blankets (MB) of a single variable.

## Requirements
* **gcc** (with C++14 support)  
This project has been tested only on Linux platform, using gcc with C++14 support.
* **Boost**  
Boost libraries are used for parsing the command line options, logging, and a few other purposes.  
More information about Boost libraries can be found at [http://boost.org/](http://boost.org/).
* **SCons**  
SCons, a cross-platform Python based build environment, is required for building the project.  
More information about the tool can be found at [http://scons.org/](http://scons.org/).
* **SABNAtk**  
SABNAtk library is used for executing counting queries over the datasets.  
More information about the library can be found at [https://gitlab.com/SCoRe-Group/SABNAtk](https://gitlab.com/SCoRe-Group/SABNAtk).
* **C++ Utils**  
Logging functionality from the [cpp-utils](https://github.com/asrivast28/cpp-utils) repository is used.
* **Google Test**   
Google Test framework is used for unit testing in this project.  
This dependency can be avoided by not building the unit tests. See "Building" section for more information.   
More information about the framework can be found at [https://github.com/google/googletest](https://github.com/google/googletest).

## Building
After all the dependencies have been installed, the project can be built as:  
<pre><code>scons
</code></pre>  
This will create an executable named `mb`, which can be used for MB discovery. 
The unit tests are built by default. The following can be executed for building only the executable  
<pre><code>scons TEST=0
</code></pre>  

Path to the Boost libraries can be specified as:  
<pre><code>scons BOOSTLIBPATH=&lt;Boost library path, default is /usr/lib/x86_64-lib-gnu&gt;
</code></pre>

For building the debug version of the executable, following can be executed:
<pre><code>scons DEBUG=1
</code></pre>  
Debug version of the executable is named `mb_debug`.

## Execution
Once the project has been built, the executable can be used for MB discovery as follows.
<pre><code>./mb -f test/coronary.csv -n 6 -m 1841 -t Smoking ## Expected Output: M. Work,P. Work,Pressure,Proteins,
</code></pre>  
Please execute the following for more information on all the options that the executable accepts.
<pre><code>./mb --help 
</code></pre>
