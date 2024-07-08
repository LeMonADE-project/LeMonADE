[![Build Status(CircleCI)](https://circleci.com/gh/LeMonADE-project/LeMonADE.svg?style=svg)](https://circleci.com/gh/LeMonADE-project/LeMonADE)
[![DOI](https://zenodo.org/badge/47562779.svg)](https://zenodo.org/badge/latestdoi/47562779)
# LeMonADE 
The abbreviation LeMonADE stands for 
"**L**attice-based **e**xtensible **Mon**te-Carlo **A**lgorithm and **D**evelopment **E**nvironment".

The aim of the LeMonADE-project is an open source implementation of the 
[bond fluctuation model](https://en.wikipedia.org/wiki/Bond_fluctuation_model) ([BFM1], [BFM2]) for simulating polymeric systems using generic
template metaprogramming in `C++`. 

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

* The options `-DINSTALLDIR_LEMONADE` and `-DBUILDDIR` for the configure script are 
  optional. The build directory defaults to `./build`, the install directory defaults
  to `/usr/local/` . The same goes for the option `-DINSTALLDIR_LEMONADE` when invoking 
  cmake (the second compilation and installation way).
* When installing using `make install`, the library is installed to
  `/given/installation/path/lib/LeMonADE/`, or if no path was specified to
  `/usr/local/lib/LeMonADE/`. Similar for the header files, which go to
  `/given/installation/path/include/LeMonADE/` or `/usr/local/include/LeMonADE/`.
* For uninstalling simply remove the files created in the installation step.
* If you do not specify the install directory, and the default `/usr/local/` is used,
  you need root access when installing. In this case

````sh
    # installs the library and headers
    sudo make install #only if you want to install the software after build
````
* If you also want to compile and run the tests, there is another option, which can be
  passed to either the configure script (first build method), or to cmake (second build method)
  This option is `-DLEMONADE_TESTS=ON` . You need internet access, because the process will
  download the googletest library. After compilation you find the executable `LeMonADE-tests`
  in your `BUILDDIR/tests`.  
  `./LeMonADE-tests --gtest_list_tests` prints all tests.  
  `./LeMonADE-tests` runs all test modules.  
  `./LeMonADE-tests --gtest_filter=TestFeatureLinearForce*` runs a specific test module.  
* Another option that can be passed is `-DCMAKE_BUILD_TYPE=Release/Debug`. The default value 
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

or

```sh
    # generates the docs
    mkdir build
    cd build
    cmake ..
    make docs
```


## Contributions

* Any contributions/suggestions making the program become better are welcome
* Every contribution should be made by following the guidelines made by the [gitflow workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)


## Authors

Find the information about active developers, former contributors, and people who contributed in the [AUTHORS](AUTHORS.md) file.

## References

* *"Structural Characterization of Model Gels under Preparation Conditions and at Swelling Equilibrium"*; R. Scholz, M. Lang; [Macromolecules 57, 2539-2555 (2024)](https://doi.org/10.1021/acs.macromol.3c02607)
* *"On the Swelling of Polymer Network Strands"*; M. Lang, R. Scholz; [Macromol. Rapid Commun 45, 2400025 (2024)](https://doi.org/10.1002/marc.202400025)
* *"Machine learning of an implicit solvent for dynamic Monte Carlo simulations"*; A. Checkervarty, J.-U. Sommer, M. Werner; [J. Chem. Phys. 158, 124904 (2023)](https://doi.org/10.1063/5.0116618)
* *"Theory of chain walking catalysis: From disordered dendrimers to dendritic bottle-brushes"*; R. Dockhorn, J.-U. Sommer; [J. Chem. Phys. 157, 044902 (2022)](https://doi.org/10.1063/5.0098263)
* *"On the Reference Size of Chains in a Network and the Shear Modulus of Unentangled Networks Made of Real Chains"*; M. Lang, T. M&uuml;ller; [Macromolecules 55, 8950-8959 (2022)](https://doi.org/10.1021/acs.macromol.2c00660)
* *"Elasticity of Tendomer Gels"*; T. M&uuml;ller, J.-U. Sommer, M. Lang; [Macromolecules 55, 7540-7555 (2022)](https://doi.org/10.1021/acs.macromol.2c01229)
* *"Amphiphilic Model Networks Based on PEG and PCL Tetra-arm Star Polymers with Complementary Reactivity"*; C. Bunk, L. Löser, N. Fribiczer, H. Komber, L. Jakisch, R. Scholz, B. Voit, S. Seiffert, K. Saalwächter, M. Lang, F. Böhme; [Macromolecules 55, 6573-6589 (2022)](https://doi.org/10.1021/acs.macromol.2c00693)
* *"Swelling and Residual Bond Orientations of Polymer Model Gels: The Entanglement-Free Limit"*; M. Lang, R. Scholz, L. Löser, C. Bunk, N. Fribiczer, S. Seiffert, F. Böhme, K. Saalwächter; [Macromolecules 55, 5997-6014 (2022)](https://doi.org/10.1021/acs.macromol.2c00589)
* *"Multimolecular Structure Formation with Linear Dendritic Copolymers"*; M. Wengenmayr, R. Dockhorn, J.-U. Sommer; [Macromolecules 54, 6937–6946 (2021)](https://doi.org/10.1021/acs.macromol.1c00226)
* *"Swelling of Tendomer Gels"*; T. M&uuml;ller, J.-U. Sommer, M. Lang; [Macromolecules 54, 4601-4614 (2021)](https://doi.org/10.1021/acs.macromol.1c00258)
* *"Analysis of the Gel Point of Polymer Model Networks by Computer Simulations"*; M. Lang, T. M&uuml;ller; [Macromolecules 53, 498-512 (2020)](https://doi.org/10.1021/acs.macromol.9b02217)
* *"Polyolefins Formed by Chain Walking Catalysis—A Matter of Branching Density Only?"*; R. Dockhorn, L. Pl&uuml;schke, M. Geisler, J. Zessin, P. Lindner, R. Mundil, J. Merna, J.-U. Sommer, A. Lederer; [J. Am. Chem. Soc. 141, 15586–15596 (2019)](https://pubs.acs.org/doi/10.1021/jacs.9b06785)
* *"Tendomers – force sensitive bis-rotaxanes with jump-like deformation behavior"*; T. M&uuml;ller, J.-U. Sommer, M. Lang; [Soft Matter, 15, 3671-3679  (2019)](https://doi.org/10.1039/C9SM00292H)
* *"Dendrimers in Solution of Linear Polymers: Crowding Effects"*; M. Wengenmayr, R. Dockhorn, J.-U. Sommer; [Macromolecules 52, 2616–2626 (2019)](https://doi.org/10.1021/acs.macromol.9b00010)
* *"Testing the physics of knots with a Feringa nanoengine"*; M. Lang, C. Schuster, R. Dockhorn, M. Wengenmayr, J.-U. Sommer; [Phys. Rev. E 98, 052501 (2018)](https://doi.org/10.1103/PhysRevE.98.052501)
* *"Formation and stabilization of pores in bilayer membranes by peptide-like amphiphilic polymers"*;  A. Checkervarty, M. Werner and J.-U. Sommer; [Soft Matter 14, 2526-2534 (2018)](https://doi.org/10.1039/C7SM02404E)
* *"Swelling Behavior of Single-Chain Polymer Nanoparticles: Theory and Simulation"*; H. Rabbel, P. Breier, J.-U. Sommer; [Macromolecules 50, 7410–7418 (2017)](https://doi.org/10.1021/acs.macromol.7b01379)
* *"A Highly Parallelizable Bond Fluctuation Model on the Body-Centered Cubic Lattice"*; C. Jentzsch, R. Dockhorn, J.-U. Sommer; [Parallel Processing and Applied Mathematics. Lecture Notes in Computer Science, vol 9574. pp 301-311 (2016)](http://dx.doi.org/10.1007/978-3-319-32152-3_28)
* *"Multicore Unimolecular Structure Formation in Single Dendritic–Linear Copolymers under Selective Solvent Conditions"*; M. Wengenmayr, R. Dockhorn, J.-U. Sommer; [Macromolecules 49, 9215-9227 (2016)](http://dx.doi.org/10.1021/acs.macromol.6b01712)
* *"Interactions of Amphiphilic Triblock Copolymers with Lipid Membranes: Modes of Interaction and Effect on Permeability Examined by Generic Monte Carlo Simulations"*; H. Rabbel, M. Werner, J.-U. Sommer; [Macromolecules 48, 4724-4732 (2015)](http://dx.doi.org/10.1021/acs.macromol.5b00720)
* *"Polymer brushes in explicit poor solvents studied using a new variant of the bond fluctuation model"*; C. Jentzsch, J.-U. Sommer; [J. Chem. Phys. 10, 104908 (2014)](http://dx.doi.org/10.1063/1.4895555)    


## License

See the [LICENSE](LICENSE) in the root directory.

## Changelog

Major changes are mentioned in the [CHANGELOG](CHANGELOG.md) file.
