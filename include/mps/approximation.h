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

    mps_root_status status;

    /**
     * @brief Attributes that have been set on
     * the roots.
     */
    mps_root_attrs  attrs;

    /**
     * @brief Inclusion status of the root
     * in the target set specified in the field
     * <code>input_config->search_set</code>.
     */
    mps_root_inclusion inclusion;

  };


#ifdef __cplusplus
}
#endif

#endif /* #ifndef MPS_ROOT_H */
