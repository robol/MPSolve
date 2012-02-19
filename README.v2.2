                        MPSolve, Version 2.2 
             by Dario Andrea Bini and Giuseppe Fiorentino
               (bini@dm.unipi.it   fiorent@dm.unipi.it)
                             May 2001

  Copyright 2001, Dipartimento di Matematica, Universita' di Pisa, Italy

  1. This software is offered without warranty of any kind, either expressed
     or implied.  The authors would appreciate, however, any reports of bugs
     or other difficulties that may be encountered.
  2. If modifications or enhancements are made to this software by others, 
     Dipartimento di Matematica reserves the right to obtain this enhanced 
     software at no cost and with no restrictions on its usage.
  3. The authors and Dipartimento di Matematica should be acknowledged in 
     any published paper based on computations using this software.  
     Please send a copy of such papers to the authors.
  4. Free use of this software is allowed for scientific purposes only.
     No commercial use is allowed unless explicitly authorized by the
     authors by means of a written agreement.


PACKAGE CONTENT

The MPSolve package contains the following:

      Directories and scripts
Data         The directory with the test polynomials
Results      The directory with some output results
Doc          This directory contains documentation in LaTeX, HTML,
             Postscript and PDF format

Makefile     The main makefile, you may need to customize this
             to fit your operating system
Maketest     A script file to run unisolve on all test polynomials


      C code
uni_solv.c   The unisolve main
mps_*.c      The MPSolve library
mps.h        The MPSolve library header file

usr_*.c      Sample user defined polynomials (rename as mps_user.c)

rur*         RUR related source code


      Extended arithmetic libraries
mt.*, mpc.*, link.*, tools.*, mptemp.*, gmptools.*


PRE-INSTALLATION REQUIREMENTS

MPSolve requires a GNU compatible "make" and C compiler (gcc) which
are standard in any linux distribution and are available on most other
platforms.  It also requires that the GMP multiprecision package has
been correctly installed and that both the header file "gmp.h" and the
library file "libgmp.a" can be accessed by GCC. These files can also
be copied to the Gmp subdirectory once the GMP package has been
successfully configured and compiled (see the GMP files README and
INSTALL for a quick install of the library).


NOTES ON GMP VERSIONS

MPSolve does not include the GMP package which can be downloaded from
the URL http://www.swox.com/gmp/

The package has been developed and fully tested with many GMP versions
from 2.0.2 up to the latest 3.1.1 (at the time of writing); however,
the former contains some known bugs which may lead to slow and badly
formatted outputs. These errors can be fixed by applying some patches
which can be downloaded from the above URL.
We suggest to always download and use the latest available version.


INSTALLATION

In order to create the executable file (unisolve), make sure you have a
GNU-compatible 'make' program (it is called 'gmake' on some systems),
then simply type:

  make

in the MPSolve directory. After a successful compilation, all
unneeded object files and libraries can be removed by issuing:

  make clean

Finally, the executable file can be tested with the command:

  make check
 
which invokes unisolve on a simple polynomial.


USAGE

In order to run the program, type

 unisolve options input_file

where 'input_file' is the name of the file holding the input polynomial,
and 'options' is the list of the options such as the goal of the computation,
the requested output precision, the subset of the complex plane
where the roots are sought and many others (see below for the full list).

If the 'input_file' is missing then unisolve will read from the standard
input stream (typically the keyboard).


INPUT FORMAT

The input file should contain:

1- Some optional comment lines starting with "!", for instance

    !  This is the Wilkinson polynomial of degree n
    !  defined by p(x)=(x-1)(x-2)...(x-n)

2- Three characters related to the nature of the polynomial:
    -the first can take the values
       d (dense polynomial)
       s (sparse polynomial)
       u (user polynomial)

    -the second can take the values
       r (real coefficients)
       c (complex coefficients)

    -the third can take the values
       i (integer coefficients)
       q (rational coefficients)
       f (float coefficients)

   For instance, "dri" stands for a dense real polynomial with
   integer coefficients

3- The precision (in decimal digits) of the coefficients, where 0 denotes
   the infinite precision typical of integer and rational coefficients.

4- The degree of the polynomial

5  - If the polynomial is coded as dense, the list of the coefficients
     arranged in order of increasing degree;

   - if the polynomial is coded as sparse, the number of nonzero
     coefficients followed by the list of all nonzero monomials coded
     as degree and coefficient written on different lines.

   - User defined polynomials must be assigned as a C program in the file
     mps_user.c (see the report mpsolve.ps in the Doc directory for more
     details).

   For instance, the polynomial x^4+2x+5 can be stored in dense form
   in the following way

     !  p(x)=x^4+2x+5
     dri
     0
     4
       5
       2
       0
       0
       1

   or, in sparse form as

     !  p(x)=x^4+2x+5
     sri
     0
     4
     3
     0
       5
     1
       2
     4
       1

   In the above examples we have indented the coefficient values for
   a better readability  (extra spaces are not important).
   Indeed, polynomials like x^500+1 can be easily written in sparse
   form as

     ! p(x)=x^500+1
     sri
     0
     500
     2
     0
       1
     500
       1

   Complex values of the coefficients are coded as real part and
   imaginary part written on different lines.
   Rational numbers are coded as numerator and denominator written on
   different lines.

   Complex numbers with rational real/imaginary parts are coded as a
   quadruple:

     numerator of the real part
     denominator of the real part
     numerator of the imaginary part
     denominator of the imaginary part.

   Below we report some other example of polynomials where the coefficients
   are complex/float/rational.

   ! Example 1,  complex/integer
   ! p(x)=x^200+(1+3i)x^70 +1
   sci
   0
   200
   3
   0
     1
     0
   70
     1
     3
   200
     1
     0


   ! Example 2,  real/rational
   ! p(x)=(3/4)x^200 + (1/3)x^70 + 1
   srq
   0
   200
   3
   0
     1
     1
   70
     1
     3
   200
     3
     4


   ! Example 3,  complex/rational
   ! p(x)=x^2 + (3/2+i*5/7)x + 2
   dcq
   0
   2
    2
    1
    0
    1

    3
    2
    5
    7

    1
    1
    0
    1


   ! Example 4,  real/float with 15 digits of precision
   ! p(x)=1.234e200 x^3 + 1.5e-200 x^2 + 1.3e-200
   drf
   15
   3
     1.3e-200
     0.0e0
     1.5e-200
     1.234e200


OPTIONS

unisolve can be used with several options, in the following list
a "*" denotes the default value:

- goal:
           * isolate the roots:             -Gi
           - approximate the roots          -Ga
           - count the roots                -Gc

- set: chooses the set where roots are searched.
           * all the complex plane          -Sa
           - left halfplane                 -Sl
           - right halfplane                -Sr
           - upper halfplane                -Su
           - lower halfplane (down)         -Sd
           - inside the unit disk           -Si
           - outside the unit disk          -So
           - Real line                      -SR
           - Imaginary line                 -SI

- multiplicity: enables or disables multiplicity detection
           - Detect multiplicity of roots   -M+
           * Do not detect multiplicity     -M-

- precision:
           - number of maximum output
             digits of the roots            -o# where # is the number
                                                of digits (default 30)
           - number of input digits
             of the coefficients
             (override the file info)       -i# where # is the number
                                                 of digits
- detection:
           - Detect real roots              -Dr
           - Detect imaginary roots         -Di
           - Detect both                    -Db
           * Do not detect                  -Dn

- output format:
           * compact form: (re, im)         -Oc
           - bare format: re \tab im        -Ob
           - "Gnuplot" low prec.: re im     -Og
           - verbose: Root(i) = re + I im   -Ov
           - full information:              -Of
             (re, im), inclusion radius, 
             status of the approximation.

- debug:  write auxiliary information
          to the standard error stream      -d
          to the standard output stream     -d1

- Parameters: (for expert users)
          max. number of packet iterations (default 1000):
                                            -Lp# 
                          
          max. number of iterations per packet (default 10):
                                            -Li#  
	  where # is the number of iterations. 

An example of usage is

 unisolve -o1000 -Ga -SR  my_poly

that computes at least 1000 digits of the real roots of my_poly, or

  unisolve -o100000 -Gc -Si < my_poly > result

which counts the number of roots of my_poly inside the unit disk using
at most 100000 digits in order to decide if a root is in or out.
The result is written in the file "result".

All files in the directory Results have been obtained with the options:

  -Gi -o1000
  
i.e., with the goal 'isolate' and an output precision of at most 1000
digits.

In the Makefile the user will also find some other options affecting
the way the program is compiled and behaves with respect to the use of
randomization and custom memory managers.

DOCUMENTATION

The main theoretical results behind MPSolve can be found in the paper:
"Design, Analysis, and Implementation of a Multiprecision Polynomial
Rootfinder" by D. A. Bini and G. Fiorentino published in Numerical
Algorithms, Volume 23 (2000) pages 127-173.

The HTML pages in the Doc directory are the main source of documentation
for the MPSolve package; these can be accessed with any browser starting
from the page:

  file://<MPSolve dir>/Doc/index.htm

where <MPSolve dir> is the directory where the MPSolve package has been
installed.

More information on the advanced features of MPSolve and on the polynomial
test suite can be found in the document "mpsolve" available in LaTeX,
PostScript and PDF format in the Doc directory.
