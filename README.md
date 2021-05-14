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

Precise dependencies have not been exhaustively tested. 

### Installing

 1. Install the Boost library, either using your system package manager or manually. Ensure this is in the global C++ include path.
 2. Clone the [TOML11](http://github.com/toruniina/toml11) library into a directory of your choice (if you clone it into MCTrans/toml11 skip 
 the next step)
 3. Make sure the MCTrans/toml11 symlink points to your toml11 installation.
 4. Run `make MCTrans++`. 

To compile the manual you will need a copy of pdflatex, bibtex, and revtex4. A pre-compiled copy of the manual is also distributed in this repository.

### Example Configurations

Example configurations live in the `examples/` subdirectory. For each example there is a corresponding example output file named `.report`. 
There is also a `check_examples.sh` script that will check all the examples still work.

## Built With

* [Boost](http://boost.org) - C++ Template library that radically extends the STL
* [TOML11](http://github.com/toruniina/toml11) - For parsing configuration files written in [TOML](https://github.com/toml-lang/toml)

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

