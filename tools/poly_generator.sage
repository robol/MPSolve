#!
# -*- coding: utf-8 -*-
#
# Polynomial generator using sage

import subprocess, random, sys

# We need imaginary unit over Q
I = QuadraticField (-1, 'I').gen ()

random.seed()

class PolynomialGenerator():

    def __init__(self):
        self._degree = 0
        self._x = var('x')
        self._pol = 1
        self._eps = Rational("1 / 100000000")

        self._comment =  "! Polygen generated polynomial\n"
        self._comment += "! This polynomial has the following roots: \n"

    def add_complex_root (self, complex_root, mult = 1):
        self._pol *= (self._x - complex_root)**mult
        self._comment += "! #) %s\n" % str(complex_root)

    def add_complex_conjugate_roots (self, complex_root, mult = 1):
        self._pol *= (self._x - complex_root)**mult
        self._pol *= (self._x - complex_root.conjugate ())**mult
        
    def add_cluster (self, cluster_center, n):
        for i in range(n):
            gamma = Rational(RR(cos(2*i*pi/n))) + Rational(RR(sin(2*i*pi/n))) * I
            gamma = gamma * self._eps
            self.add_complex_root (cluster_center + gamma)

    def get_pol_file (self):
        self._pol = expand (self._pol)

        coeffs = list (map (lambda x : x[0], self._pol.coefficients()))
        imag_coeffs = list (map (lambda x : x.imag (), coeffs))
                       
        real = reduce (lambda x, y : x and y, map (lambda x : x == 0, imag_coeffs))

        pol_file = self._comment + "\n"
        pol_file += "Degree = %d;\n" % (len (coeffs) - 1)

        # Check if the coefficients are real
        if real:
            pol_file += "Real;\n"
        else:
            pol_file += "Complex;\n"
            
        # Determine poly structure
        pol_file += "Rational;\n"
        pol_file += "Monomial;\n"
        pol_file += "\n"

        for c in coeffs:
            r = c.real ()
            i = c.imag ()

            pol_file += str (r) + " "
            if not real:
                pol_file += str (i) + " "

            pol_file += "\n"

        return pol_file
               
            
def usage ():
    print "%s degree [type]" % sys.argv[0]
    print ""
    print " type can be one of: "
    print "  - separate: Well separate integer roots (Wilkinson-style), default choice"
    print "  - clustered: Polynomial with clusters of roots"
    print "  - mixed: Polynomial with both clusters and isolated roots"
    print "  - multiple: Polynomial with multiple roots"
    print "  - all: Polynomial with all the above"
    print ""
    sys.exit (1)

if __name__ == "__main__":

    if (len (sys.argv) < 2 or len (sys.argv) > 3):
        usage ()

    polygen = PolynomialGenerator ()
    n = int(sys.argv[1])

    t = "separate"
    if (len (sys.argv) == 3):
        t = sys.argv[2]

    integer_grid = range (1, 100)

    if t == "separate":
        for i in range(n):
            polygen.add_complex_root (random.choice (integer_grid) + 
                                      random.choice (integer_grid) * I)
    elif t == "clustered":
        degree = 0
        while degree < n:
            step = random.choice (range (1, n - degree + 1))
            polygen.add_cluster (random.choice (integer_grid) + 
                                 random.choice (integer_grid) * I, step)
            degree += step
    elif t == "mixed":
        degree = 0
        while degree < n:
            if random.random() <= 0.7:
                polygen.add_complex_root (random.choice (integer_grid) + 
                                          random.choice (integer_grid) * I)
                degree += 1
            else:
                step = random.choice (range (1, n - degree + 1))
                polygen.add_cluster (random.choice (integer_grid) + 
                                     random.choice (integer_grid) * I, step)
                degree += step

    elif t == "multiple":
        pass

    elif t == "all":
        pass

    else:
        usage ()
            
    print polygen.get_pol_file ()
