#ifndef MPS_ROOT_H_
#define MPS_ROOT_H_

#ifdef __cplusplus
extern "C" {
#endif

  struct mps_approximation {

    cplx_t fvalue;
    cdpe_t dvalue;
    mpc_t  mvalue;

    double frad;
    rdpe_t drad;

    mps_boolean approximated;
    mps_boolean again;
    long int wp;

  };


#ifdef __cplusplus
}
#endif

#endif /* #ifndef MPS_ROOT_H */
