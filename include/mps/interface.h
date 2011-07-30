/* 
 * File:   interface.h
 * Author: leonardo
 *
 * Created on 23 aprile 2011, 10.35
 */

#ifndef MPS_INTERFACE_H
#define	MPS_INTERFACE_H

#ifdef	__cplusplus
extern "C" {
#endif

/* Octave module workardound */
#ifdef __UNDEF_CPLUSPLUS
#undef __cplusplus
#endif

#include <gmp.h>
#include <mps/types.h>
#include <mps/mt.h>
#include <mps/mpc.h>
#include <mps/tools.h>
#include <stdlib.h>
#include <stdio.h>


/*
 *    ====== ROUTINES EXPOSED TO THE INTERFACE ======
 */

/* functions in mps_defaults.c */
void mps_set_default_values(mps_status* s);

/* Functions in mps_main.c */
void mps_mpsolve(mps_status* s);
void mps_standard_mpsolve(mps_status* s);

/* functions in mps_interface.c */
mps_status* mps_status_new();
void mps_status_free(mps_status* s);
int mps_status_set_poly_d(mps_status* s, cplx_t* coeff, long unsigned int n);
int mps_status_set_poly_i(mps_status* s, int* coeff, long unsigned int n);
int mps_status_get_roots_d(mps_status* s, cplx_t* roots, double* radius);
int mps_status_set_poly_u(mps_status* s, int n, mps_fnewton_ptr fnewton,
		mps_dnewton_ptr dnewton, mps_mnewton_ptr mnewton);
void mps_status_allocate_poly_inplace(mps_status* s, int n);
void mps_status_select_algorithm(mps_status* s, mps_algorithm algorithm);
void mps_status_set_degree(mps_status* s, int n);
int mps_status_get_roots_d(mps_status* s, cplx_t* roots, double* radius);
int mps_status_get_roots_m(mps_status* s, mpc_t* roots, rdpe_t* radius);

#ifdef	__cplusplus
}
#endif

#ifdef __UNDEF_CPLUSPLUS
}
#endif

#endif	/* MPS_INTERFACE_H */

