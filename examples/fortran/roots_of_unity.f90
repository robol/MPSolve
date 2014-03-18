! Fortran program that computes root of unity
! using MPSolve routines
!
! Author: Leonardo Robol <robol@mail.dm.unipi.it>
!
! Compile this file with gfortran linking the module
! mps, libmps, libgmp, libpthread and libm
!
!
PROGRAM roots_of_unity

	! Double precision
	INTEGER, PARAMETER :: dp = KIND(0.d0)
	
	! Degree of the polynomial
	INTEGER, PARAMETER :: n = 5
	
	! Auxiliary variables
	INTEGER :: i
	
	! Coefficients and roots
	COMPLEX(dp), DIMENSION(n + 1) :: coeff
	COMPLEX(dp), DIMENSION(n) :: roots

	! Set coefficients of x^n - 1
	coeff = 0
	coeff(1) = -1
	coeff(n + 1) = 1

	! Call MPSolve routine
	CALL MPS_ROOTS(n, coeff, roots)
	
	! Write out roots
	WRITE(*,*) "Printing computed roots of unity for n = ", n
	DO i = 1,n
		WRITE(*,*) roots(i)
	END DO

END PROGRAM
