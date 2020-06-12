MPSOLVE 3.2.1
=============

[![Build Status](https://travis-ci.org/robol/MPSolve.svg?branch=master)](https://travis-ci.org/robol/MPSolve)

 MPSolve is a C package to solve polynomials and secular equations. It released under the terms
 of the GNU General public license as it specified in the COPYING file inside the source directory.

 Additional information can be found on the [official website](https://numpi.dm.unipi.it/software/mpsolve)

### How to install MPSolve

If you have downloaded an official MPSolve tarball you can 
install MPSolve simply by typing these commands in a shell:

    ./configure 
    make 
    [sudo] make install

whilst, if you have checked out the git repository directly
you have to use the script `autogen.sh` first, and then the
usual `configure`, `make`, `make install` sequence. 

The last command is optional and install mpsolve system-wide. You can simply
use the mpsolve executable built in its directory by launching 
`./src/mpsolve/mpsolve` from the source directory.

Full install is needed to use MPSolve as a library in other C, 
FORTRAN, Matlab, ... software without further tweaking. 

The `examples/` folder contains a mix of example source files that use
MPSolve and bindings for other programming languages such as Python,
Octave, Matlab (TM), ...

### Using the mpsolve binary

The mpsolve binary is thought as a simple way to solve polynomials 
and/or secular equation using a text file as input. A generic 
input file for mpsolve is composed by a preamble and 
a body. The comments are identified by lines starting with '!'.

In the preamble will be specified all the options for solving 
and some general information about the polynomial. 
The body will contain the coefficients.

Every option in the preamble has to be specified as Key; or Key=value;
This is an example of a valid input for mpsolve that specifies 
the polynomial x^5 - 1

    ! File: nroots5.pol 
    Degree=5;
    Monomial;
    Real;
    Integer;
    
    -1
     0
     0
     0
     0
     1
    
    ! EOF

 The same polynomial can be specified by using sparse notation:

    ! File nroots5sparse.pol
    Degree=5;
    Monomial;
    Real;
    Integer;
    Sparse;
    
    5  1  ! Highest degree coefficient
    0  -1 ! Coefficient of degree 0
    
    ! EOF

 If the `Real` option is not specified real and imaginary part of the coefficients
 will be needed as input. Rational input is also used in here:

    ! File: random-poly.pol
    Degree=3;
    Monomial;
    Rational;

    45/9 7/4
    3/23 293/34234
    234/2369234 2348234/324
    324 234324/23

