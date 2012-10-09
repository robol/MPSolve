/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **               Version 2.9, September 2012              **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2012, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>

mps_approximation *
mps_approximation_new (mps_context * s)
{
  mps_approximation * appr = mps_new (mps_approximation);
  mpc_init2 (appr->mvalue, s->mpwp);
  appr->again = true;
  appr->approximated = false;
  return appr;
}

void
mps_approximation_free (mps_context * s, mps_approximation * appr)
{
  mpc_clear (appr->mvalue);
  free (appr);
}
