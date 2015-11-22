This program simulates an evolving viral population where the epidemiology
is described by SIRS dynamics and the molecular evolution is given by
codons mutating according to the Kimura two-parameter model.
This is described in more detail in [here](http://dx.doi.org/10.1101/027995)

Install
-------
This program uses the R-standalone library which is included in the R source files
Note that it is not installed by default

    > cd R-<version>
    > ./configure

This step is system-dependent but problems can be inconsistent libaries
This can be resolved by explicitly defining the architecture, i.e

    > ./configure r_arch=x86_64 CC='gcc -arch x86_64' CXX='g++ -arch x86_64'
          F77='gfortran -arch x86_64' FC='gfortran x86_64'
          OBJC='gcc -arch x86_64' --with-x=no

Once the configuration step is performed, makefiles will be generated in
all the sub-directories. For the standalone library, we don't need the top-level files

    > cd src/nmath/standalone
    > make

If the install occurs properly this should generate 'libRmath.a'
The first lines in the Makefile should be modified to match the location
of the R library.

Once the library installed properly, simEpiCodon can be compiled by typing 'make'
in the directory containing the source code

Usage
-----
To run the simulation type

    > ./epiCodon -o <outfile> [options]

