# README #

How to get the DVM code up and running

### What is this repository for? ###

The DVM code is a research code initially created by Aaron Morris and continued by Peter Clarke. The code is a discrete velocity method for solving the Boltzmann equation.

* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

The DVM code requires several libraries/programs before it can be compiled and used:

* autotools
* boost - libboost regex-dev
* grvy
    * Provided in the repository. Open buildgrvy.sh and edit the environment variables and pathnames to match your system then run the file.

After downloading and installing the libraries, open your bashrc file and add the following lines:

* export COMPILER=<compiler>
* export COMPILER_VERSION=<version>
* export BOOST_VERSION=<version>
* export GRVY_VERSION=0.29.1
* export BOOST_DIR=<dir path>
* export GRVY_DIR=<dir path>/grvy-$GRVY_VERSION-$COMPILER-$COMPILER_VERSION-boost-$BOOST_VERSION

To configure the code navigate into trunk and enter the following commands:

* ./bootstrap
* ./configure --prefix=<install directory> <options>

Suggested options for configure include builds for:

* *Optimized:* FCFLAGS='-O3'
* *Simple testing:* --enable-debug=yes FCFLAGS='-g -O0 -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow'
* *Full testing:* --enable-debug=yes FCFLAGS='-g -O0 -fbacktrace -fbounds-check -Wall -Wextra -pedantic -fimplicit-none -ffpe-trap=invalid,zero,overflow,underflow'
* *Profiling:* --enable-debug=yes FCFLAGS='-pg -O3'

After configuring is complete run:

* make
* make install

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Peter Clarke - peter.clarke29@gmail.com