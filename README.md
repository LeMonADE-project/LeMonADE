# LeMonADE
The abbreviation LeMonADE stands for 
"**L**attice-based **e**xtensible **Mon**te-Carlo **A**lgorithm and **D**evelopment **E**nvironment".

The aim of the LeMonaDE-project is an open source implementation of the 
bond-fluctuation-model [1] [BFM1] [2] [BFM2] for simulating polymeric systems using generic
template metaprogramming in C++. 

[BFM1]: http://dx.doi.org/10.1021/ma00187a030  "I. Carmesin, K. Kremer; Macromolecules 21, 2819-2823 (1988)"
 
[BFM2]: http://dx.doi.org/10.1063/1.459901 "H. P. Deutsch, K. Binder; J. Chem. Phys. 94, 2294-2304 (1990)"


## Installation

* Clone `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install cmake (minimum version 2.6.2)
* Just do for standard compilation (library+examples):

````sh
    # generates the application in build-directory
    ./configure -DINSTALLDIR_LEMONADE=/path/to/install/LeMonADE/ -DBUILDDIR=/path/to/build/ 
    make
    make install #only if you want to install the software after build
````

 or
 
````sh
    # generates the lib and examples
    mkdir build
    cd build
    cmake -DINSTALLDIR_LEMONADE=/path/to/install/LeMonADE/ ..
    make
    make install #only if you want to install the software after build
````

* The options -DINSTALLDIR_LEMONADE and -DBUILDDIR for the configure script are 
  optional. The build directory defaults to ./build, the install directory defaults
  to /usr/local/ . The same goes for the option -DINSTALLDIR_LEMONADE when invoking 
  cmake (the second compilation and installation way).
* When installing using 'make install', the library is installed to
  /given/installation/path/lib/LeMonADE/, or if no path was specified to
  /usr/local/lib/LeMonADE/. Similar for the header files, which go to
  /given/installation/path/include/LeMonADE/ or /usr/local/include/LeMonADE/.
* For uninstalling simply remove the files created in the installation step.
* If you do not specify the install directory, and the default /usr/local/ is used,
  you need root access when installing. In this case

````sh
    # installs the library and headers
    sudo make install #only if you want to install the software after build
````
* If you also want to compile and run the tests, there is another option, which can be
  passed to either the configure script (first build method), or to cmake (second build method)
  This option is -DLEMONADE_TESTS=ON . You need internet access, because the process will
  download the googletest library
* Another option that can be passed is -DCMAKE_BUILD_TYPE=Release/Debug. The default value 
  is Release, which uses compiler flags for optimization. If the option "Debug" is chosen,
  no compiler optimizations are used, compiler warnings are enabled by -Wall, and 
  debugging information is compiled into the binary with -g. Other flags then "Release",
  or "Debug" are not allowed and will lead to termination and an error message by cmake.

## Getting Started

This can be found in the documentation.


## Build the documentation

* Install doxygen 
* Just do for documentation, only:

````sh
    # generates the documentation
    ./configure -DINSTALLDIR_LEMONADE=/path/to/install/LeMonADE/ -DBUILDDIR=/path/to/build/ 
    make docs
````


```sh
    # generates the docs
    mkdir build
    cd build
    cmake ..
    make docs
```
## References

This library has been used in the following publications (without guarantee for completeness)
* C. Jentzsch, J.-U. Sommer; J. Chem. Phys. 10, 104908 (2014)    
* H. Rabbel, M. Werner, J.-U. Sommer; Macromolecules 48, 4724-4732 (2015)
* M. Lang, M. Werner, R. Dockhorn, T. Kreer; Macromolecules 49, 5190-5201 (2016)
* M. Wengenmayr, R. Dockhorn, J.-U. Sommer; Macromolecules 49, 9215-9227 (2016)

## License

See the LICENSE in the root directory.
