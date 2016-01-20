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
    ./configure -DINSTALLDIR_LEMONADE=/path/to/install/LeMonADE/
    make
````

 or
 
````sh
    # generates the lib and examples
    mkdir build
    cd build
    cmake ..
    make
````

## Getting Started

This can be find in the documentation.


## Build the documentation

* Install doxygen 
* Just do for documentation, only:

```sh
    # generates the docs
    mkdir build
    cd build
    cmake ..
    make docs
```
    
## License

See the LICENSE in the root directory.
