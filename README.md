# MCTrans++ -- A 0D modelling tool for centrifugal mirrors

MCTrans++ is a C++ reimplementation of the original MCTrans Mathematica tool developed by Prof. A. Hassam of
University of Maryland at College Park to model the Maryland Centrifugal Experiment.

MCTrans++ uses a variety of simplified models for the behaviour of a centrifugally-confined mirorr plasma to provide a capability 
for rapid scoping of experimental operating points. This tool is not intended to replace detailed ab initio modelling -- benchmarking studies 
will be listed in this README to provide some sense of where this tool is known to be representative of the detailed physics. This tool
cannot provide quantitatively accurate results for any given experiment! If you are very lucky these results may be good to 50%!

Detailed analyses that lead to the formulae implemented in the code are provided in the document in the notes/ subdirectory of this distribution.

## Getting Started

You will need to download this codebase and compile it in order to run MCTrans++. 

### Prerequisites

To compile and use MCTrans++ you will need a system with the following

 - A C++17 compliant C++ compiler.
 - The Boost C++ Template Library
 - The TOML11 library
 - The SUNDIALS library

Precise dependencies have not been exhaustively tested. 

### Building MCTrans++

All the build options are set in the file `Makefile.local`, which you need to provide for your system.
An example is provided in `Makefile.local.example` -- copy this file to `Makefile.local` and make any edits needed. 
This file is in GNU-compatable Makefile format, and you can set and override all the compilation options here.
For example, if you are not using the default compiler (g++), then you can add a line to `Makefile.local` that reads `CXX = /path/to/my/c++/compiler`.
 
If you're happy with this, let's proceed!

 1. Install the Boost library, either using your system package manager or manually by downloading from [here](https://www.boost.org). If this is a system-wide install, 
 proceed to step 2. If you downloaded the Boost libraries, add a line to `Makefile.local` which sets `BOOST_DIR = /path/to/boost`.
 2. Clone the [TOML11](http://github.com/toruniina/toml11) library into a directory of your choice. If you clone it into the default location of MCTrans/toml11, proceed to step 3. As with Boost, set `TOML11_DIR = /path/to/toml11` in `Makefile.local`.
 4. Install [SUNDIALS](https://computing.llnl.gov/projects/sundials) and edit Makefile.local to set `SUNDIALS_DIR` to the location you have installed the Sundials library in. If you are only using SUNDIALS for MCTrans++, a quick intro to installing SUNDIALS is inclued below.
 5. Set any other options, e.g. setting the variable `DEBUG` to any value will build a version that you can use to develop MCTrans++ and that includes debug information.
 6. Run `make MCTrans++`. 

To compile the manual you will need a copy of pdflatex, bibtex, and revtex4. A pre-compiled copy of the manual is also distributed in this repository.

#### Installing SUNDIALS 

If this is your first use of SUNDIALS, a quick guide to installing the base libraries follows here.

Pick where you want the sundials sources / build tree / compiled libraries to go. We will call these directories
SUNDIALS_SOURCE, SUNDIALS_BUILD, and SUNDIALS_INSTALL in the following. One suggestion would be
```
SUNDIALS_SOURCE = ~/sundials/source
SUNDIALS_BUILD  = ~/sundials/build
SUNDIALS_INSTALL = ~/sundials/
```

With these directories picked, we can download and compile SUNDIALS.

 1. Download the SUNDIALS source from [here](https://computing.llnl.gov/projects/sundials) or [here](https://github.com/LLNL/sundials) into `SUNDIALS_SOURCE`
 2. Move to `SUNDIALS_BUILD`. Configure the SUNDIALS build with 
 ```
 cmake $SUNDIALS_SOURCE -DCMAKE_INSTALL_PREFIX=$SUNDIALS_INSTALL -DEXAMPLES_INSTALL=off
 ```
	   If this gives you any errors (lack of C compiler, etc), refer to the SUNDIALS documentation.
 3. Compile SUNDIALS with `make -j install`. 
 4. You now have sundials installed into the `SUNDIALS_INSTALL` directory. This is the path you should set `SUNDIALS_DIR` to in your MCTrans `Makefile.local`


### Example Configurations

Example configurations live in the `examples/` subdirectory. For each example there is a corresponding example output file named `.report`. 
There is also a `check_examples.sh` script that will check all the examples still work.

## Built With

* [Boost](http://boost.org) - C++ Template library that radically extends the STL
* [TOML11](http://github.com/toruniina/toml11) - For parsing configuration files written in [TOML](https://github.com/toml-lang/toml)
* [Sundials](https://computing.llnl.gov/projects/sundials) - Suite of libraries from Lawrence Livermore National Laboratory for numerical solution of Nonlinear Algebraic Equations, ODEs and DAEs

## Known Issues

General issues include the fact that the steady-state equations may have multiple roots.

Specific known issues are listed here.

 - Debug builds currently fail on OSX due to the use of `feenableexcept()`.

## Contributing

Contributions to this project are welcome. However, as we are currently in active development, please email the authors if you wish to get involved.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/ianabel/MCTrans/tags). 

## Authors

* **Ian Abel** - *C++ Version* - [Ian Abel at UMD](https://ireap.umd.edu/faculty/abel)
* **Adil Hassam** - *Original MCTrans Code*

For full copyright attribution, see the [COPYRIGHT](COPYRIGHT) file.
For a summary of contributors, see the [contributors](http://github.com/ianabel/MCTrans/contributors) page.

## License

This project is licensed under the 3-Clause BSD Licence - see the [LICENSE.md](LICENSE.md) file for details

