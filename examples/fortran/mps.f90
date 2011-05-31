MODULE MPS
!
! Module used to interface Fortran 95 routines with MPSolve
!
! The function to use to call the polynomial solver is
! MPS_ROOTS(n, coeff, roots)
!
! PARAMETERS:
!  INTEGER :: n							Degree of the polynomial
! 
!  COMPLEX(2), DIMENSION(n+1) coeff		Array with the complex coefficients of the
!										polynomial. 
!
!  COMPLEX(2), DIMENSION(n)   roots		Preallocated array where the roots will be
!										saved. 
!
!
CONTAINS
	SUBROUTINE MPS_ROOTS(n, coeff, roots)
		IMPLICIT NONE
		INTEGER, PARAMETER :: dp = KIND(0.d0)
		INTEGER :: n, i
		COMPLEX(dp), DIMENSION(n)   :: roots
		COMPLEX(dp), DIMENSION(n+1) :: coeff

		! Call C bridge with MPSolve
		CALL MPS_ROOTS_IMPL (coeff, roots, n)
	END SUBROUTINE
END MODULE
