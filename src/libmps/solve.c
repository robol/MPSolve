/***********************************************************
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.2, May 2001                  **
 **                                                        **
 **                      Written by                        **
 **       Dario Andrea Bini and Giuseppe Fiorentino        **
 **       (bini@dm.unipi.it)  (fiorent@dm.unipi.it)        **
 **                                                        **
 ** (C) 2001, Dipartimento di Matematica, FRISCO LTR 21024 **
 ***********************************************************/

#include <math.h>
#include <mps/gmptools.h>
#include <float.h>
#include <mps/core.h>
#include <mps/threading.h>

/**
 * @brief Set <code>again[i]</code> to <code>true</code> or to <code>false</code> 
 * according to the values of <code>status[i,*]</code> and <code>goal</code>.
 *
 * More precisely:
 *
 * - If goal is "count": .true. only for statu='**u' but not for 
 *         'f*u', 'a*u', 'o*u'                            (1) 
 *      - multipl. on: true also for 'c**'
 *      - Real on    : true also for '*u*' excluded (1)
 *      - Imag. on   : true also for '*v*' excluded (1)
 *
 * Goal isolate: true only for statu='**u' or statu='c**' but
 *          not for 'f**', 'a**', 'o**', 'i*i', 'i*o'     (2) 
 *      - multipl. on: true also for 'c**'
 *      - Real on    : true also for '*u*' excluded (2)
 *      - Imag. on   : true also for '*v*' excluded (2)
 * Goal approximate: the same as isolate.
 */
void
mps_update (mps_status * s)
{
  int i;

  for (i = 0; i < s->n; i++)
    s->again[i] = false;
  switch (s->goal[0])
    {

    case 'c':                  /*  count */
      for (i = 0; i < s->n; i++)
        {
          if (s->status[i][2] == 'u')
            if (s->status[i][0] != 'f' && s->status[i][0] != 'a'
                && s->status[i][0] != 'o')
              s->again[i] = true;
          if (s->goal[2] == 'm' && s->status[i][0] == 'c' && s->status[i][0]
              != 'o')
            s->again[i] = true;

          switch (s->goal[3])
            {

            case 'r':          /* real option */
              if (s->status[i][1] == 'w' && (s->status[i][2] != 'u'
                                             || (s->status[i][0] != 'f'
                                                 && s->status[i][0] != 'a'
                                                 && s->status[i][0] != 'o')))
                s->again[i] = true;
              break;

            case 'i':          /* imaginary option */
              if (s->status[i][1] == 'w' && (s->status[i][2] != 'u'
                                             || (s->status[i][0] != 'f'
                                                 && s->status[i][0] != 'a'
                                                 && s->status[i][0] != 'o')))
                s->again[i] = true;
              break;

            case 'b':          /* both imaginary and real options */
              if (s->status[i][1] == 'w' && (s->status[i][2] != 'u'
                                             || (s->status[i][0] != 'f'
                                                 && s->status[i][0] != 'a'
                                                 && s->status[i][0] != 'o')))
                s->again[i] = true;
              break;
            }
        }

      break;

    case 'i':                  /* isolate */
      for (i = 0; i < s->n; i++)
        {
          if (s->status[i][2] == 'u' || (s->status[i][0] == 'c'
                                         && s->status[i][2] == 'i'))
            if (s->status[i][0] != 'f' && s->status[i][0] != 'a'
                && s->status[i][0] != 'o' && (s->status[i][0] != 'i'
                                              || s->status[i][2] != 'i'))
              s->again[i] = true;
          if (s->goal[2] == 'm' && s->status[i][0] == 'c' && s->status[i][2]
              != 'o')
            s->again[i] = true;

          switch (s->goal[3])
            {

            case 'r':          /* real option */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;
              break;

            case 'i':          /* imaginary option */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;     /* DARIO RIVEDERE */
              break;

            case 'b':          /* both imaginary and real options */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;
              break;
            }
        }

      break;

    case 'a':                  /* approximate (the same as isolate) */
      for (i = 0; i < s->n; i++)
        {
          if (s->status[i][2] == 'u' || (s->status[i][0] == 'c'
                                         && s->status[i][2] == 'i'))
            if (s->status[i][0] != 'f' && s->status[i][0] != 'a'
                && s->status[i][0] != 'o')
              s->again[i] = true;

          if (s->goal[2] == 'm' && s->status[i][0] == 'c' && s->status[i][2]
              != 'o')
            s->again[i] = true;

          switch (s->goal[3])
            {

            case 'r':          /* real option */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;
              break;

            case 'i':          /* imaginary option */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;
              break;

            case 'b':          /* both imaginary and real options */
              if (s->status[i][1] == 'w' && (s->status[i][0] != 'f'
                                             && s->status[i][0] != 'a'
                                             && s->status[i][0] != 'o'))
                s->again[i] = true;
              break;
            }
        }
      break;
    }
}

/**
 * @brief Compute super center and super radius
 * 
 * This routines the super radius of the <code>i</code>-th cluster,
 * i.e. the radius of the inclusion disc for the whole cluster
 *
 * @param s The <code>mps_status</code> associated with the current computation.
 * @param i Index of the cluster;
 * @param sc Center of the cluster;
 * @param sr Double that will be set to the super radius of the cluster;
 */
void
mps_fsrad (mps_status * s, int i, cplx_t sc, double *sr)
{
  cplx_t ctmp;
  double sum;
  int j, l;

  sum = 0.0;
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      sum += s->frad[l];
    }
  cplx_set (sc, cplx_zero);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      cplx_mul_d (ctmp, s->froot[l], s->frad[l]);
      cplx_add_eq (sc, ctmp);
    }
  cplx_div_eq_d (sc, sum);
  *sr = 0.0;
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      cplx_sub (ctmp, sc, s->froot[l]);
      *sr = MAX (*sr, s->frad[l] + cplx_mod (ctmp));
    }
}

/**
 * @brief <code>dpe</code> version of <code>fsrad()</code>
 */
void
mps_dsrad (mps_status * s, int i, cdpe_t sc, rdpe_t sr)
{
  cdpe_t ctmp;
  rdpe_t sum, rtmp;
  int j, l;

  rdpe_set (sum, rdpe_zero);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      rdpe_add_eq (sum, s->drad[l]);
    }
  cdpe_set (sc, cdpe_zero);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      cdpe_mul_e (ctmp, s->droot[l], s->drad[l]);
      cdpe_add_eq (sc, ctmp);
    }
  cdpe_div_eq_e (sc, sum);
  rdpe_set (sr, rdpe_zero);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      cdpe_sub (ctmp, sc, s->droot[l]);
      cdpe_mod (rtmp, ctmp);
      rdpe_add_eq (rtmp, s->drad[l]);
      if (rdpe_lt (sr, rtmp))
        rdpe_set (sr, rtmp);
    }
}

/**
 * @brief Multiprecision versione of <code>fsrad()</code>
 */
void
mps_msrad (mps_status * s, int i, mpc_t sc, rdpe_t sr)
{
  int j, l;
  rdpe_t rtmp;
  cdpe_t cdtmp;
  tmpf_t ftmp, sum;
  tmpc_t ctmp;

  tmpc_init2 (ctmp, s->mpwp);
  tmpf_init2 (ftmp, s->mpwp);
  tmpf_init2 (sum, s->mpwp);

  mpf_set_ui (sum, 0);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      mpf_set_rdpe (ftmp, s->drad[l]);
      mpf_add (sum, sum, ftmp);
    }

  mpc_set_ui (sc, 0, 0);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      mpf_set_rdpe (ftmp, s->drad[l]);
      mpc_mul_f (ctmp, s->mroot[l], ftmp);
      mpc_add_eq (sc, ctmp);
    }

  mpc_div_eq_f (sc, sum);
  rdpe_set (sr, rdpe_zero);
  for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
    {
      l = s->clust[s->punt[i] + j];
      mpc_sub (ctmp, sc, s->mroot[l]);
      mpc_get_cdpe (cdtmp, ctmp);
      cdpe_mod (rtmp, cdtmp);
      rdpe_add_eq (rtmp, s->drad[l]);
      if (rdpe_lt (sr, rtmp))
	  rdpe_set (sr, rtmp);
      else
	{
	  MPS_DEBUG_RDPE (s, sr, "sr");
	  MPS_DEBUG_RDPE (s, rtmp, "rtmp");
	}
    }

  tmpf_clear (sum);
  tmpf_clear (ftmp);
  tmpc_clear (ctmp);
}

/**
 * @brief Modify the vector 'status' according to the goal, and
 * to the location of the roots.
 *
 * @param s The mps_status associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * The subroutine is used also for marking the new cluster
 * that have been detected between two consecutive packets
 * of Aberth's iteration.
 * 
 * -# The subroutine changes into 'C' the components of
 * status[:1) corresponding to old clusters, keeping
 * status[:,1)='c' for the new formed clusters.
 * In this way applying restart selects new starting
 * approximations only for the new detected clusters.
 *
 * -# For the components for which 
 * status[:,1]!='C', 'f', 'x' performs the following
 * analysis:
 * If the cluster has mult=1 mark it with status[:1)='i'
 * if is also approximated mark it with status(:1)='a'
 * Check if c*u and i*u (i.e., uncertain set) can 
 * be made certain according to goal[1]
 *
 * -# Perform the same with options, that is,
 * If multiplicity is on then check if a cluster
 * corresponds to a  multiple root
 * If detect real then detect real roots
 * if detect imaginary then detect imaginary roots
 * If detect both then detect both imaginary and
 * real roots
 */
void
mps_fmodify (mps_status * s, mps_boolean track_new_cluster)
{
  int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, nf, j2, l2;
  double sr, tmpr, afri, sep1;
  cplx_t sc;
  rdpe_t rtmp;
  mps_boolean tcr, tcr1;

  /* ==1== Change into 'C' the components of status for old clusters */
  nf = 2 * s->n;                /* Isolation factor */

  /* Mark old cluster with 'C' if requested */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->status[i][0] == 'c')
	  s->status[i][0] = 'C';
    }

  i_old = 0;
  i_new = 0;
  for (i = 1; i <= s->nclust && i_new < s->n; i++)
    {                           /*  loop1: DO i=1, nclust */
      if (s->oldpunt[i_old + 1] == s->punt[i_new + 1])
        {
          i_old++;
          i_new++;
          continue;
        }
      else
        {
          for (j = i_new + 1; j < s->nclust; j++)
            {                   /* loop2: DO j=i_new+1, nclust */
              if (s->oldpunt[i_old + 1] != s->punt[j + 1])
                continue;
              else
                {
                  nnewclust = j - i_new + 1;    /* scan each new cluster */
                  for (k = 0; k < nnewclust; k++)
                    {           /* loop3: DO k=1, nnewclust */
                      i1 = i_new + k;
                                                /*********************************
						 scan the entries of each new cluster
						 set status[l][0]='i' if the cluster has multip=1
						 and mark with 'c'
						 the ones which are different from 'i'
						 **********************************/
                      if (s->punt[i1 + 1] - s->punt[i1] == 1
                          && s->status[s->clust[s->punt[i1]]][0] != 'x'
                          && s->status[s->clust[s->punt[i1]]][0] != 'f')
                        s->status[s->clust[s->punt[i1]]][0] = 'i';
                      for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++)
                        {       /* loop4: */
                          ip1 = s->clust[s->punt[i1] + l];
                          if (s->status[ip1][0] != 'i' && s->status[ip1][0]
                              != 'x' && s->status[ip1][0] != 'f'
                              && s->status[ip1][0] != 'a'
                              && s->status[ip1][0] != 'o')
                            s->status[ip1][0] = 'c';
                        }
                    }
                  i_new = j + 1;
                  i_old++;
                  break;
                }               /* else */
            }                   /* for */
        }                       /* else */
    }                           /* for */

        /*=2== Scan all the clusters */
  for (i = 0; i < s->nclust; i++)
    {                           /* scan : DO i=1,s->nclust */
      /* check isolation/approximation */
      if (s->punt[i + 1] - s->punt[i] == 1
          && s->status[s->clust[s->punt[i]]][0] != 'x'
          && s->status[s->clust[s->punt[i]]][0] != 'f')
        {
          s->status[s->clust[s->punt[i]]][0] = 'i';
          tmpr = cplx_mod (s->froot[s->clust[s->punt[i]]]);
          tmpr = s->frad[s->clust[s->punt[i]]] / tmpr;
          tmpr = log (tmpr);
          if (tmpr < -s->output_config->prec * LOG2)
            s->status[s->clust[s->punt[i]]][0] = 'a';
        }
      /* Scan inside the cluster */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan_in */
          l = s->clust[s->punt[i] + j];
                        /**************************************************
			 first check for inside/outside unit disk in the case where
			 there are very large and/or very small roots (statu='x')
			 and for counting only
			 *************************************************/
          afri = cplx_mod (s->froot[i]);
	  
	  /* Check if the root, even if clustered, is approximated */
	     rdpe_set_d (rtmp, s->frad[l] / cplx_mod (s->froot[l]));     
	     if (rdpe_le (rtmp, s->eps_out))     
	       s->status[l][0] = 'o';     

          if (s->status[l][0] == 'x' && s->goal[0] == 'c')
            {
              if ((s->goal[1] == 'i' && afri < 1) || (s->goal[1] == 'o'
                                                      && afri > 1))
                s->status[l][2] = 'i';
              if ((s->goal[1] == 'i' && afri > 1) || (s->goal[1] == 'o'
                                                      && afri < 1))
                s->status[l][2] = 'o';
            }
          /* Now check the standard cases */
          if (s->status[l][0] == 'x')
            continue;
          if ((s->status[l][0] == 'c' || s->status[l][0] == 'i'
               || s->status[l][0] == 'C' || s->status[l][0] == 'o')
              && s->status[l][2] == 'u')
            {
              /* Check if the approximation is inside/outside the set */
              switch (s->goal[1])
                {

                case 'a':      /* all */
                  s->status[l][2] = 'i';
                  break;

                case 'i':      /* inside unit circle */
                  if (!mps_ftouchunit (s, nf, l))
                    {
                      if (cplx_mod (s->froot[l]) < 1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'o':      /* outside unit circle */
                  if (!mps_ftouchunit (s, nf, l))
                    {
                      if (cplx_mod (s->froot[l]) > 1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'l':      /* left half plane  */
                  if (!mps_ftouchimag (s, nf, l))
                    {
                      if (cplx_Re (s->froot[l]) < 0)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'r':      /* right half plane */
                  if (!mps_ftouchimag (s, nf, l))
                    {
                      if (cplx_Re (s->froot[l]) > 0)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'u':      /* upper half plane */
                  if (!mps_ftouchreal (s, nf, l))
                    {
                      if (cplx_Im (s->froot[l]) > 0)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'd':      /* lower half plane */
                  if (!mps_ftouchreal (s, nf, l))
                    {
                      if (cplx_Im (s->froot[l]) < 0)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'R':      /* Real line  NEW */
                  if (s->status[l][1] != 'w')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_ftouchreal (s, 1, l))
                        {
                          if (s->data_type[1] == 'r')
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'R';
                            }
                          else
                            {
                              /* fsrad(i, sc, &sr); DARIO */
                              sr = s->frad[l];
                              if (log (sr) < s->sep - s->n * s->lmax_coeff)
                                {
                                  s->status[l][2] = 'i';
                                  s->status[l][1] = 'R';
                                }
                              else
                                {
                                  s->status[l][2] = 'u';
                                  s->status[l][1] = 'w';
                                }
                            }
                        }
                      else
                        {       /* do not touch real */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'r';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_ftouchreal (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_ftouchreal (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /*  mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano R */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_fsrad (s, i, sc, &sr);
                              sep1 = s->sep;
                              if (s->data_type[1] == 'c')
                                sep1 = s->sep - s->n * s->lmax_coeff;
                              if (log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'R';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano R */

                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][2] = 'o';
                              s->status[l2][1] = 'r';
                            }
                        }
                    }
                  break;

                case 'I':      /* Imaginary line  NEW */
                  if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_ftouchimag (s, nf, l))
                        {
                          /* fsrad(i, sc, &sr); DARIO */
                          sr = s->frad[l];
                          if (log (sr) < s->sep - s->n * s->lmax_coeff)
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'I';
                            }
                          else
                            {
                              s->status[l][2] = 'u';
                              s->status[l][1] = 'w';
                            }
                        }
                      else
                        {       /* do not touch imag */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_ftouchimag (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_ftouchimag (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /* mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano I */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_fsrad (s, i, sc, &sr);
                              sep1 = s->sep;
                              sep1 = s->sep - s->n * s->lmax_coeff;
                              if (log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'I';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano I */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                    }
                  break;

                case 'S':      /* Set provided by the user */
                  mps_error (s, 1, "User set not implemented yet");
                  break;

                default:
                  mps_error (s, 1, "Mistake in goal");
                  break;
                }
            }
        }
      /* If some cluster still contains an uncertain disk then set
       * all the disks uncertain */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {
          l = s->clust[s->punt[i] + j];
          if (s->status[l][2] == 'u')
            {
              for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                {
                  l1 = s->clust[s->punt[i] + j1];
                  s->status[l1][2] = 'u';
                }
              break;
            }
        }
    scan:;                     /* scan the next component */
    }

        /*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < s->nclust; i++)
    {                           /* scan1 */
      if (s->punt[i + 1] - s->punt[i] == 1)
        continue;
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan1_in */
          l = s->clust[s->punt[i] + j];
          if (s->status[l][0] == 'x' || s->status[l][0] == 'f'
              || s->status[l][0] == 'm')
            goto scan1;
          if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /* NEW */
                                    s->status[l][0] == 'C'))
            {                   /* multiplicity on */
              if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /* float coeff. */
                mps_error (s, 1,
                           "Fatal: Float coefficients - impossible to detect multiplicity");

              /* compute super center and super radius */
              mps_fsrad (s, i, sc, &sr);

              if (log (sr) < s->sep)
                for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                  {
                    l1 = s->clust[s->punt[i] + j1];     /* NEW j-> j1 */
                    s->status[l1][0] = 'm';
                  }
              goto scan1;
            }
        }
    scan1:;                    /* scan next component */
    }

  /* Option Real check */
  if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan2 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan2_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][1] != 'w' && s->status[l][1] != 'i')
                continue;
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan2;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_ftouchreal (s, nf, l))
                    {
                      if (s->data_type[1] == 'r')
                        s->status[l][1] = 'R';
                      else
                        {
                          /* fsrad(i, sc, &sr); DARIO */
                          sr = s->frad[l];
                          if (log (sr) < s->sep - s->n * s->lmax_coeff)
                            s->status[l][1] = 'R';
                          else
                            s->status[l][1] = 'w';
                        }
                    }
                  else
                    /* do not touch real */
                    s->status[l][1] = 'r';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_ftouchreal (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_ftouchreal (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan2_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano R */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_fsrad (s, i, sc, &sr);
                          sep1 = s->sep;
                          if (s->data_type[1] == 'c')
                            sep1 = s->sep - s->n * s->lmax_coeff;
                          if (log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'R';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano R */
                    s->status[l][1] = 'r';
                }
            scan2_in:;
            }
        scan2:;
        }
    }
  /* Option Imaginary check */
  if (s->goal[3] == 'i' || s->goal[3] == 'b')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan3 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan3_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan3;
              if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                continue;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_ftouchimag (s, nf, l))
                    {
                      /* fsrad(i, sc, &sr); DARIO */
                      sr = s->frad[l];
                      if (log (sr) < s->sep - s->n * s->lmax_coeff)
                        s->status[l][1] = 'I';
                      else
                        s->status[l][1] = 'w';
                    }
                  else
                    /* do not touch imag */
                    s->status[l][1] = 'i';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_ftouchimag (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_ftouchimag (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan3_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano I */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_fsrad (s, i, sc, &sr);
                          sep1 = s->sep;
                          sep1 = s->sep - s->n * s->lmax_coeff;
                          if (log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'I';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano I */
                    s->status[l][1] = 'i';
                }
            scan3_in:;
            }
        scan3:;
        }
    }
}

/**
 * @brief The DPE version of <code>mps_fmodify()</code>.
 *
 * @param s The mps_status associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * @see mps_fmodify()
 */
void
mps_dmodify (mps_status * s, mps_boolean track_new_cluster)
{
  int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, j2, l2, nf;
  double rtmp, sep1;
  rdpe_t sr, tmpr, tmpr2;
  cdpe_t sc;
  mps_boolean tcr, tcr1;

  /* ==1==  Change into 'C' the components of status for old clusters */
  nf = 2 * s->n;                /* Isolation factor */

  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->status[i][0] == 'c')
	  s->status[i][0] = 'C';
    }

  i_old = 0;
  i_new = 0;
  for (i = 1; i <= s->nclust && i_new < s->n; i++)
    {                           /*  loop1: DO i=1, s->nclust */
      if (s->oldpunt[i_old + 1] == s->punt[i_new + 1])
        {
          i_old++;
          i_new++;
          continue;
        }
      else
        {
          for (j = i_new + 1; j < s->nclust; j++)
            {                   /* loop2:  */
              if (s->oldpunt[i_old + 1] != s->punt[j + 1])
                continue;
              else
                {
                  nnewclust = j - i_new + 1;    /*  scan each new cluster */
                  for (k = 0; k < nnewclust; k++)
                    {           /* loop3:  */
                      i1 = i_new + k;

		      /* Check if the root, even if clustered, is approximated */
		      rdpe_set (tmpr, s->drad[i1]);
		      cdpe_mod (tmpr2, s->droot[i1]);
		      rdpe_div_eq (tmpr, tmpr2);
		      if (rdpe_le (tmpr, s->eps_out)) 
			s->status[i1][0] = 'o'; 

                      /* scan the entries of each new cluster
                       * set status[l][0]='i' if the cluster has multip=1
                       * and mark with 'c' those which are different from 'i'
                       */
                      if (s->punt[i1 + 1] - s->punt[i1] == 1
                          && s->status[s->clust[s->punt[i1]]][0] != 'x'
                          && s->status[s->clust[s->punt[i1]]][0] != 'f')
                        s->status[s->clust[s->punt[i1]]][0] = 'i';
                      for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++)
                        {       /* loop4: */
                          ip1 = s->clust[s->punt[i1] + l];
                          if (s->status[ip1][0] != 'i' && s->status[ip1][0]
                              != 'x' && s->status[ip1][0] != 'f'
                              && s->status[ip1][0] != 'a'
                              && s->status[ip1][0] != 'o')
                            s->status[ip1][0] = 'c';
                        }
                    }
                  i_new = j + 1;
                  i_old++;
                  break;
                }
            }
        }
    }

        /*=2== Scan all the clusters */
  for (i = 0; i < s->nclust; i++)
    {                           /* scan: */
      /* check isolation/approximation */
      if (s->punt[i + 1] - s->punt[i] == 1
          && s->status[s->clust[s->punt[i]]][0] != 'x'
          && s->status[s->clust[s->punt[i]]][0] != 'f')
        {
          s->status[s->clust[s->punt[i]]][0] = 'i';
          cdpe_mod (tmpr, s->droot[s->clust[s->punt[i]]]);
          rdpe_div (tmpr, s->drad[s->clust[s->punt[i]]], tmpr);
          rtmp = rdpe_log (tmpr);
          if (rtmp < -s->output_config->prec * LOG2)
            s->status[s->clust[s->punt[i]]][0] = 'a';
        }
      /* Scan inside the cluster */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan_in: */
          l = s->clust[s->punt[i] + j];

          /* Now check the standard cases */
          if (s->status[l][0] == 'x')
            continue;
          if ((s->status[l][0] == 'c' || s->status[l][0] == 'i'
               || s->status[l][0] == 'C' || s->status[l][0] == 'o')
              && s->status[l][2] == 'u')
            {
              /* Check if the approximation is inside/outside the set */
              switch (s->goal[1])
                {

                case 'a':      /* all */
                  s->status[l][2] = 'i';
                  break;

                case 'i':      /* inside unit circle */
                  if (!mps_dtouchunit (s, nf, l))
                    {
                      cdpe_mod (tmpr, s->droot[l]);
                      if (rdpe_lt (tmpr, rdpe_one))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'o':      /* outside unit circle */
                  if (!mps_dtouchunit (s, nf, l))
                    {
                      cdpe_mod (tmpr, s->droot[l]);
                      if (rdpe_gt (tmpr, rdpe_one))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'l':      /* left half plane  */
                  if (!mps_dtouchimag (s, nf, l))
                    {
                      rdpe_set (tmpr, cdpe_Re (s->droot[l]));
                      if (rdpe_lt (tmpr, rdpe_zero))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'r':      /* right half plane */
                  if (!mps_dtouchimag (s, nf, l))
                    {
                      rdpe_set (tmpr, cdpe_Re (s->droot[l]));
                      if (rdpe_gt (tmpr, rdpe_zero))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'u':      /* upper half plane */
                  if (!mps_dtouchreal (s, nf, l))
                    {
                      rdpe_set (tmpr, cdpe_Im (s->droot[l]));
                      if (rdpe_gt (tmpr, rdpe_zero))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'd':      /* lower half plane */
                  if (!mps_dtouchreal (s, nf, l))
                    {
                      rdpe_set (tmpr, cdpe_Im (s->droot[l]));
                      if (rdpe_lt (tmpr, rdpe_zero))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'R':      /* Real line  NEW */
                  if (s->status[l][1] != 'w')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_dtouchreal (s, 1, l))
                        {
                          if (s->data_type[1] == 'r')
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'R';
                            }
                          else
                            {
                              /* dsrad(i, sc, sr); DARIO */
                              rdpe_set (sr, s->drad[l]);
                              if (rdpe_log (sr) < s->sep - s->n
                                  * s->lmax_coeff)
                                {
                                  s->status[l][2] = 'i';
                                  s->status[l][1] = 'R';
                                }
                              else
                                {
                                  s->status[l][2] = 'u';
                                  s->status[l][1] = 'w';
                                }
                            }
                        }
                      else
                        {       /* do not touch real */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'r';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_dtouchreal (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_dtouchreal (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /*  mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano R */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_dsrad (s, i, sc, sr);
                              sep1 = s->sep;
                              if (s->data_type[1] == 'c')
                                sep1 = s->sep - s->n * s->lmax_coeff;
                              if (rdpe_log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'R';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano R */

                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][2] = 'o';
                              s->status[l2][1] = 'r';
                            }
                        }
                    }
                  break;

                case 'I':      /* Imaginary line  NEW */
                  if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_dtouchimag (s, nf, l))
                        {
                          /* dsrad(i, sc, sr); */
                          rdpe_set (sr, s->drad[l]);
                          if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'I';
                            }
                          else
                            {
                              s->status[l][2] = 'u';
                              s->status[l][1] = 'w';
                            }
                        }
                      else
                        {       /* do not touch imag */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_dtouchimag (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_dtouchimag (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /* mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano I */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_dsrad (s, i, sc, sr);
                              sep1 = s->sep;
                              sep1 = s->sep - s->n * s->lmax_coeff;
                              if (rdpe_log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'I';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano I */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                    }
                  break;

                case 'S':      /* Set provided by the user */
                  mps_error (s, 1, "Custom region not implemented yet");
                  break;
                default:
                  mps_error (s, 1, "mistake in goal");
                  break;
                }
            }
        }
      /* If some cluster still contains an uncertain disk then set
       * all the disks uncertain
       */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {
          l = s->clust[s->punt[i] + j];
          if (s->status[l][2] == 'u')
            {
              for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                {
                  l1 = s->clust[s->punt[i] + j1];
                  s->status[l1][2] = 'u';
                }
              break;
            }
        }
    scan:;
    }

        /*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < s->nclust; i++)
    {                           /* scan1 */
      if (s->punt[i + 1] - s->punt[i] == 1)
        continue;
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan1_in */
          l = s->clust[s->punt[i] + j];
          if (s->status[l][0] == 'x' || s->status[l][0] == 'f'
              || s->status[l][0] == 'm')
            goto scan1;
          if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /* NEW */
                                    s->status[l][0] == 'C'))
            {                   /* multiplicity on */

              if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /* float coeff. */
                mps_error (s, 1,
                           "Fatal: Float coefficients - impossible to detect multiplicity");

              /* compute super center and super radius */
              mps_dsrad (s, i, sc, sr);

              if (rdpe_log (sr) < s->sep)
                for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                  {
                    l1 = s->clust[s->punt[i] + j1];     /* NEW j-> j1 */
                    s->status[l1][0] = 'm';
                  }
              goto scan1;
            }
        }
    scan1:;                    /* scan next component */
    }

  /* Option Real check */
  if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan2 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan2_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][1] != 'w' && s->status[l][1] != 'i')
                continue;
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan2;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_dtouchreal (s, nf, l))
                    {
                      if (s->data_type[1] == 'r')
                        s->status[l][1] = 'R';
                      else
                        {
                          /* dsrad(i, sc, sr); DARIO */
                          rdpe_set (sr, s->drad[l]);
                          if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                            s->status[l][1] = 'R';
                          else
                            s->status[l][1] = 'w';
                        }
                    }
                  else
                    /* do not touch real */
                    s->status[l][1] = 'r';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_dtouchreal (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_dtouchreal (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan2_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano R */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_dsrad (s, i, sc, sr);
                          sep1 = s->sep;
                          if (s->data_type[1] == 'c')
                            sep1 = s->sep - s->n * s->lmax_coeff;
                          if (rdpe_log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'R';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano R */
                    s->status[l][1] = 'r';
                }
            scan2_in:;
            }
        scan2:;
        }
    }

  /* Option Imaginary check */
  if (s->goal[3] == 'i' || s->goal[3] == 'b')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan3 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan3_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan3;
              if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                continue;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_dtouchimag (s, nf, l))
                    {
                      /* dsrad(i, sc, sr); DARIO */
                      rdpe_set (sr, s->drad[l]);
                      if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                        s->status[l][1] = 'I';
                      else
                        s->status[l][1] = 'w';
                    }
                  else
                    /* do not touch imag */
                    s->status[l][1] = 'i';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_dtouchimag (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_dtouchimag (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan3_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano I */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_dsrad (s, i, sc, sr);
                          sep1 = s->sep;
                          sep1 = s->sep - s->n * s->lmax_coeff;
                          if (rdpe_log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'I';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano I */
                    s->status[l][1] = 'i';
                }
            scan3_in:;
            }
        scan3:;
        }
    }
}

/**
 * @brief The multiprecision version of the routine
 * <code>mps_fmodify()</code>. 
 *
 * @param s The mps_status associated to the current computation.
 * @param track_new_cluster true if old clusters should be marked
 * with 'C' instead of 'c', so they are recognizable (for shifting).
 *
 * @see mps_fmodify()
 */
void
mps_mmodify (mps_status * s, mps_boolean track_new_cluster)
{
  int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, nf, j2, l2;
  double rtmp, sep1;
  rdpe_t sr, tmpr, tmpr2;
  cdpe_t tmpc;
  mps_boolean tcr, tcr1;
  tmpf_t tmpf;
  tmpc_t sc;

  tmpc_init2 (sc, s->mpwp);
  tmpf_init2 (tmpf, s->mpwp);

  /* ==1==  Change into 'C' the components of status for old clusters */
  nf = 2 * s->n;                /* Isolation factor */

  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->status[i][0] == 'c')
	  s->status[i][0] = 'C';
    }

  i_old = 0;
  i_new = 0;
  for (i = 1; i <= s->nclust && i_new < s->n; i++)
    {                           /* loop1: */
      if (s->oldpunt[i_old + 1] == s->punt[i_new + 1])
        {
          i_old++;
          i_new++;
          continue;
        }
      else
        {
          for (j = i_new + 1; j < s->nclust; j++)
            {                   /* loop2: */
              if (s->oldpunt[i_old + 1] != s->punt[j + 1])
                continue;
              else
                {
                  nnewclust = j - i_new + 1;    /* scan each new cluster */
                  for (k = 0; k < nnewclust; k++)
                    {           /* loop3: */
                      i1 = i_new + k;

		      /* Check if the root, even if clustered, is approximated */
		      rdpe_set (tmpr, s->drad[i1]);
		      cdpe_mod (tmpr2, s->droot[i1]);
		      rdpe_div_eq (tmpr, tmpr2);
		      if (rdpe_le (tmpr, s->eps_out)) 
			s->status[i1][0] = 'o'; 

                                                /*****************************************
						 scan the entries of each new cluster set
						 status[l][0]='i' if the cluster has multip=1 and
						 mark with 'c' the ones which are different from 'i'
						 *****************************************/
                      if (s->punt[i1 + 1] - s->punt[i1] == 1
                          && s->status[s->clust[s->punt[i1]]][0] != 'x'
                          && s->status[s->clust[s->punt[i1]]][0] != 'f')
                        s->status[s->clust[s->punt[i1]]][0] = 'i';
                      for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++)
                        {       /* loop4: */
                          ip1 = s->clust[s->punt[i1] + l];
                          if (s->status[ip1][0] != 'i' && s->status[ip1][0]
                              != 'x' && s->status[ip1][0] != 'f'
                              && s->status[ip1][0] != 'a'
                              && s->status[ip1][0] != 'o')
                            s->status[ip1][0] = 'c';
                        }
                    }
                  i_new = j + 1;
                  i_old++;
                  break;
                }
            }
        }
    }

        /*=2== Scan all the clusters */
  for (i = 0; i < s->nclust; i++)
    {                           /*  scan : DO i=1,s->nclust */
      /* check isolation/approximation */
      if (s->punt[i + 1] - s->punt[i] == 1
          && s->status[s->clust[s->punt[i]]][0] != 'x'
          && s->status[s->clust[s->punt[i]]][0] != 'f')
        {
          s->status[s->clust[s->punt[i]]][0] = 'i';
          mpc_get_cdpe (tmpc, s->mroot[s->clust[s->punt[i]]]);
          cdpe_mod (tmpr, tmpc);
          rdpe_div (tmpr, s->drad[s->clust[s->punt[i]]], tmpr);
          rtmp = rdpe_log (tmpr);
          if (rtmp < -s->output_config->prec * LOG2)
            {
              s->status[s->clust[s->punt[i]]][0] = 'a';
            }
        }
      /* Scan inside the cluster */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan_in: */
          l = s->clust[s->punt[i] + j];

          /* Now check the standard cases */
          if (s->status[l][0] == 'x')
            continue;
          if ((s->status[l][0] == 'c' || s->status[l][0] == 'i'
               || s->status[l][0] == 'a' || s->status[l][0] == 'C'
               || s->status[l][0] == 'o') && s->status[l][2] == 'u')
            {
              /* Check if the approximation is inside/outside the set */
              switch (s->goal[1])
                {

                case 'a':      /* all */
                  s->status[l][2] = 'i';
                  break;

                case 'i':      /* inside unit circle */
                  if (!mps_mtouchunit (s, nf, l))
                    {
                      mpc_get_cdpe (tmpc, s->mroot[l]);
                      cdpe_mod (tmpr, tmpc);
                      if (rdpe_lt (tmpr, rdpe_one))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'o':      /* outside unit circle */
                  if (!mps_mtouchunit (s, nf, l))
                    {
                      mpc_get_cdpe (tmpc, s->mroot[l]);
                      cdpe_mod (tmpr, tmpc);
                      if (rdpe_gt (tmpr, rdpe_one))
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'l':      /* left half plane  */
                  if (!mps_mtouchimag (s, nf, l))
                    {
                      /* mpc_get_re(tmpf, s->mroot[l]); */
                      mpf_set (tmpf, mpc_Re (s->mroot[l]));
                      if (mpf_sgn (tmpf) == -1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'r':      /* right half plane */
                  if (!mps_mtouchimag (s, nf, l))
                    {
                      /* mpc_get_re(tmpf, s->mroot[l]); */
                      mpf_set (tmpf, mpc_Re (s->mroot[l]));
                      if (mpf_sgn (tmpf) == 1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'u':      /* upper half plane */
                  if (!mps_mtouchreal (s, nf, l))
                    {
                      /* mpc_get_im(tmpf, s->mroot[l]); */
                      mpf_set (tmpf, mpc_Im (s->mroot[l]));
                      if (mpf_sgn (tmpf) == 1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'd':      /* lower half plane */
                  if (!mps_mtouchreal (s, nf, l))
                    {
                      /* mpc_get_im(tmpf, s->mroot[l]); */
                      mpf_set (tmpf, mpc_Im (s->mroot[l]));
                      if (mpf_sgn (tmpf) == -1)
                        s->status[l][2] = 'i';
                      else
                        s->status[l][2] = 'o';
                    }
                  break;

                case 'R':      /* Real line  NEW */
                  if (s->status[l][1] != 'w')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_mtouchreal (s, 1, l))
                        {
                          if (s->data_type[1] == 'r')
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'R';
                            }
                          else
                            {
                              rdpe_set (sr, s->drad[l]);
                              /* msrad(i, sc, sr);**#DARIO */
                              if (rdpe_log (sr) < s->sep - s->n
                                  * s->lmax_coeff)
                                {
                                  s->status[l][2] = 'i';
                                  s->status[l][1] = 'R';
                                }
                              else
                                {
                                  s->status[l][2] = 'u';
                                  s->status[l][1] = 'w';
                                }
                            }
                        }
                      else
                        {       /* do not touch real */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'r';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_mtouchreal (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_mtouchreal (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /*  mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano R */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_msrad (s, i, sc, sr);
                              sep1 = s->sep;
                              if (s->data_type[1] == 'c')
                                sep1 = s->sep - s->n * s->lmax_coeff;
                              if (rdpe_log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'R';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano R */

                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][2] = 'o';
                              s->status[l2][1] = 'r';
                            }
                        }
                    }
                  break;

                case 'I':
                  if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                    continue;
                  if (s->punt[i + 1] - s->punt[i] == 1)
                    {           /* one disk */
                      if (mps_mtouchimag (s, nf, l))
                        {
                          /* msrad(i, sc, sr); *//*#DARIO */
                          rdpe_set (sr, s->drad[l]);
                          if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                            {
                              s->status[l][2] = 'i';
                              s->status[l][1] = 'I';
                            }
                          else
                            {
                              s->status[l][2] = 'u';
                              s->status[l][1] = 'w';
                            }
                        }
                      else
                        {       /* do not touch imag */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                      continue;
                    }
                  else
                    {           /* cluster */
                      tcr = mps_mtouchimag (s, nf, l);
                      for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                        {
                          l1 = s->clust[s->punt[i] + j1];
                          tcr1 = mps_mtouchimag (s, nf, l1);
                          if ((tcr && tcr1) || (!tcr && !tcr1))
                            continue;
                          else
                            {   /* mixed situation */
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                              goto scan;
                            }
                        }
                      if (tcr)
                        {       /* tutti i dischi intersecano I */
                          if (s->data_type[2] == 'f' || s->data_type[2]
                              == 'b')
                            {
                              for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                   j2++)
                                {
                                  l2 = s->clust[s->punt[i] + j2];
                                  s->status[l2][2] = 'u';
                                  s->status[l2][1] = 'w';
                                }
                            }
                          else
                            {   /* integer/rational polynomial */
                              mps_msrad (s, i, sc, sr);
                              sep1 = s->sep;
                              sep1 = s->sep - s->n * s->lmax_coeff;
                              if (rdpe_log (sr) < sep1)
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'i';
                                      s->status[l2][1] = 'I';
                                    }
                                }
                              else
                                {
                                  for (j2 = 0; j2 < s->punt[i + 1]
                                       - s->punt[i]; j2++)
                                    {
                                      l2 = s->clust[s->punt[i] + j2];
                                      s->status[l2][2] = 'u';
                                      s->status[l2][1] = 'w';
                                    }
                                }
                            }
                        }
                      else
                        {       /* tutti i dischi non intersecano I */
                          s->status[l][2] = 'o';
                          s->status[l][1] = 'i';
                        }
                    }
                  break;

                case 'S':      /* Set provided by the user */
                  mps_error (s, 1, "Custom region not implemented yet");
                  break;
                default:
                  mps_error (s, 1, "mistake in goal");
                  break;
                }
            }
        }
      /* If some cluster still contains an uncertain disk then set
       * all the disks uncertain
       */
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {
          l = s->clust[s->punt[i] + j];
          if (s->status[l][2] == 'u')
            {
              for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                {
                  l1 = s->clust[s->punt[i] + j1];
                  s->status[l1][2] = 'u';
                }
              break;
            }
        }

    scan:;
    }

        /*==3==  now check the options */

  /* Option multiplicity */
  for (i = 0; i < s->nclust; i++)
    {                           /* scan1 */
      if (s->punt[i + 1] - s->punt[i] == 1)
        continue;
      for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
        {                       /* scan1_in */
          l = s->clust[s->punt[i] + j];
          if (s->status[l][0] == 'x' || s->status[l][0] == 'f'
              || s->status[l][0] == 'm')
            goto scan1;
          if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /* NEW */
                                    s->status[l][0] == 'C'))
            {                   /* multiplicity on */

              if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /* float coeff. */
                mps_error (s, 1,
                           "Fatal: Float coefficients - impossible to detect multiplicity");

              /* compute super center and super radius */
              mps_msrad (s, i, sc, sr);

              if (rdpe_log (sr) < s->sep)
                for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                  {
                    l1 = s->clust[s->punt[i] + j1];     /* NEW j-> j1 */
                    s->status[l1][0] = 'm';
                  }
              goto scan1;
            }
        }
    scan1:;                    /* scan next component */
    }

  /* Option Real check */
  if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan2 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan2_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][1] != 'w' && s->status[l][1] != 'i')
                continue;
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan2;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_mtouchreal (s, nf, l))
                    {
                      if (s->data_type[1] == 'r')
                        s->status[l][1] = 'R';
                      else
                        {
                          /* msrad(i, sc, sr); *//*#DARIO */
                          rdpe_set (sr, s->drad[l]);
                          if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                            s->status[l][1] = 'R';
                          else
                            s->status[l][1] = 'w';
                        }
                    }
                  else
                    /* do not touch real */
                    s->status[l][1] = 'r';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_mtouchreal (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_mtouchreal (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan2_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano R */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_msrad (s, i, sc, sr);
                          sep1 = s->sep;
                          if (s->data_type[1] == 'c')
                            sep1 = s->sep - s->n * s->lmax_coeff;
                          if (rdpe_log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'R';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano R */
                    s->status[l][1] = 'r';
                }
            scan2_in:;
            }
        scan2:;
        }
    }
  /* Option Imaginary check */
  if (s->goal[3] == 'i' || s->goal[3] == 'b')
    {
      for (i = 0; i < s->nclust; i++)
        {                       /* scan3 */
          for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++)
            {                   /* scan3_in */
              l = s->clust[s->punt[i] + j];
              if (s->status[l][0] == 'x' || s->status[l][0] == 'f')
                goto scan3;
              if (s->status[l][1] != 'w' && s->status[l][1] != 'r')
                continue;
              if (s->punt[i + 1] - s->punt[i] == 1)
                {               /* one disk */
                  if (mps_mtouchimag (s, nf, l))
                    {
                      rdpe_set (sr, s->drad[l]);
                      /* msrad(i, sc, sr); DARIO */
                      if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff)
                        s->status[l][1] = 'I';
                      else
                        s->status[l][1] = 'w';
                    }
                  else
                    /* do not touch imag */
                    s->status[l][1] = 'i';
                  continue;
                }
              else
                {               /* cluster */
                  tcr = mps_mtouchimag (s, nf, l);
                  for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++)
                    {
                      l1 = s->clust[s->punt[i] + j1];
                      tcr1 = mps_mtouchimag (s, nf, l1);
                      if ((tcr && tcr1) || (!tcr && !tcr1))
                        continue;
                      else
                        {       /* mixed situation */
                          for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                            {
                              l2 = s->clust[s->punt[i] + j2];
                              s->status[l2][1] = 'w';
                            }
                          goto scan3_in;
                        }
                    }
                  if (tcr)
                    {           /* tutti i dischi intersecano I */
                      if (s->data_type[2] == 'f' || s->data_type[2] == 'b')
                        for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++)
                          {
                            l2 = s->clust[s->punt[i] + j2];
                            s->status[l2][1] = 'w';
                          }
                      else
                        {       /* integer/rational polynomial */
                          mps_msrad (s, i, sc, sr);
                          sep1 = s->sep;
                          sep1 = s->sep - s->n * s->lmax_coeff;
                          if (rdpe_log (sr) < sep1)
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'I';
                              }
                          else
                            for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i];
                                 j2++)
                              {
                                l2 = s->clust[s->punt[i] + j2];
                                s->status[l2][1] = 'w';
                              }
                        }
                    }
                  else
                    /* tutti i dischi non intersecano I */
                    s->status[l][1] = 'i';
                }
            scan3_in:;
            }
        scan3:;
        }
    }

  tmpf_clear (tmpf);
  tmpc_clear (sc);
}

/**
 * @brief Check if the roots are computed with the required
 * precision.
 *
 * Set <code>computed</code> to <code>true</code> 
 * if the stop condition is satisfied,
 * otherwise set <code>computed</code> to <code>false</code>.
 *
 * The stop condition is obtained from the vector <code>status</code> 
 * as follows:
 *
 * If the <code>goal</code> is count stop if
 *  -  <code>**u</code> does not exist, except for 
 *     <code>a*u</code>, <code>o*u</code>, <code>f*u</code>
 *  - Mult. and does not exist <code>c**</code>
 *  - Real. and does not exist <code>*u*</code>, except for <code>au*</code>, 
 *    <code>ou*</code>;
 *  - Imag  and does not exist <code>*v*</code>, except for 
 *    <code>av*</code>, <code>ov*</code>;
 *
 * If the <code>goal</code> is isolate or approximate stop if:
 * - <code>**u</code> does not exist, except for  <code>a*u</code>, 
 *   <code>o*u</code>, <code>f*u</code>
 *   and if <code>c*i</code> does not exist;
 * - Mult. and does not exist <code>c*i</code>, <code>o*i</code>
 * - Real. and does not exist <code>*ui</code>, except for <code>aui</code>, 
 *   <code>oui</code>;
 * - Imag  and does not exist <code>*vi</code>, 
 *   except for <code>avi</code>, <code>ovi</code>;
 * 
 * @see status
 */
mps_boolean
mps_check_stop (mps_status * s)
{
  MPS_DEBUG_THIS_CALL;

  int i;
  mps_boolean computed;

  computed = false;
  /* count */
  if (s->goal[0] == 'c')
    {
      for (i = 0; i < s->n; i++)
        {
          if (s->status[i][2] == 'u' && s->status[i][0] != 'f'
              && s->status[i][0] != 'o' && s->status[i][0] != 'a')
            return computed;
          if (s->goal[2] == 'm' && s->status[i][0] == 'c' && s->status[i][2]
              != 'o')
            return computed;
          if (s->goal[3] == 'r' && s->status[i][1] == 'w' && s->status[i][2] != 'o' && s->status[i][0] != 'a' && s->status[i][0] != 'o' && s->status[i][0] != 'm')      /* NEW */
            return computed;
          if (s->goal[3] == 'i' && s->status[i][1] == 'w' && s->status[i][2]
              != 'o' && s->status[i][0] != 'a' && s->status[i][0] != 'o'
              && s->status[i][0] != 'm')
            return computed;
          if (s->goal[3] == 'b' && s->status[i][2] != 'o' && s->status[i][1]
              == 'w' && s->status[i][0] != 'a' && s->status[i][0] != 'o'
              && s->status[i][0] != 'm')
            return computed;
        }
      computed = true;
    }
  /* isolate or approximate */
  if (s->goal[0] == 'i' || s->goal[0] == 'a')
    {
      for (i = 0; i < s->n; i++)
        {
          if (s->status[i][2] == 'u' && s->status[i][0] != 'f'
              && s->status[i][0] != 'o' && s->status[i][0] != 'a')
            return computed;
          if (s->status[i][0] == 'c' && s->status[i][2] != 'o')
            return computed;
          if (s->goal[2] == 'm' && s->status[i][2] != 'o' && (s->status[i][0]
                                                              == 'c'))
            return computed;
          if (s->goal[3] == 'r' && s->status[i][1] == 'w' && s->status[i][2]
              != 'o' && s->status[i][0] != 'a' && s->status[i][0] != 'o'
              && s->status[i][0] != 'm')
            return computed;
          if (s->goal[3] == 'i' && s->status[i][1] == 'w' && s->status[i][2]
              != 'o' && s->status[i][0] != 'm' && s->status[i][0] != 'a'
              && s->status[i][0] != 'o' && s->status[i][0] != 'm')
            return computed;
          if (s->goal[3] == 'b' && s->status[i][2] != 'o' && s->status[i][0]
              != 'm' && s->status[i][1] == 'w' && s->status[i][0] != 'a'
              && s->status[i][0] != 'o')
            return computed;
        }
      computed = true;
    }

  return computed;
}

/**
 * @brief Actually solve the polynomial
 *
 * This routine performs the following computations:
 * -# Select starting approximations and check if some of them need dpe.
 *    Initialize the vector again which is true if the corresponding
 *    approximation is out of the root neighbourhood.
 * -# Performs max_pack packets of Aberth iterations on all the
 *    components out of the root neighbourhood belonging to the set S
 *    and on which it is possible to iterate with float.
 *    More precisely, each packet performs max_it iterations on all the
 *    components where again is true. At each iteration check if the
 *    current approximation is in the root neighbourhood; in this case
 *    set 'again' to false.  
 * -# If at the end of the general packet all the approximations
 *    are inside the root neighbourhood, i.e., 'again' is false in all
 *    the components then return.
 *    else, perform cluster analysis, select new starting approximations
 *    update the vector 'statu', update the vector 'again' that selects 
 *    the components on which to iterate, according to the goal, and 
 *    repeat until the max number of allowed packets is reached. 
 *    In the latter case output FAILURE.
 *
 * The local variable <code>again</code> controls the iteration: i.e., 
 *   <code>again[i]=true</code> means iterate on the <code>i</code>-th 
 *   component
 *
 * @param s The <code>mps_status</code> associated with the current computation.
 * @param d_after_f this variable is <code>true</code> if dpe
 * are needed after the floating point pass. 
 */
void
mps_fsolve (mps_status * s, mps_boolean * d_after_f)
{
  mps_boolean excep;
  int it_pack, iter, nit, oldnclust, i, j;
  rdpe_t eps_out;
  mps_monomial_poly *p = s->monomial_poly;

  /* == 1 ==  Initialize variables */
  it_pack = 0;
  s->nclust = 0;
  for (i = 0; i < s->n; i++)
    {
      s->again[i] = true;
      cplx_set (s->froot[i], cplx_zero);
      s->frad[i] = DBL_MAX;
    }

  /* choose starting approximations */
  if (s->DOLOG)
    fprintf (s->logstr, "FSOLVE: call fstart");

  /* If there is a custom starting point function use it,
   * otherwise use the default one */
  if (s->fstart_usr)
    (*s->fstart_usr) (s, s->n, 0, 0.0, 0.0, eps_out);
  else
    mps_fstart (s, s->n, 0, 0.0, 0.0, eps_out, p->fap);

        /***************
	 this part of code performs shift in the gravity center of the roots
	 In order to use it, uncomment the part below and comment the
	 instruction above. Dangerous for overflow.
	 ************/
  /*
     {
     cplx_t ft;
     cplx_mul_d(ft, fpc[n], -n);
     cplx_div(ft, fpc[n-1], ft);
     fshift(n, 0, 100, ft, eps_out);
     }
   *//* till here */

  if (s->DOLOG)
    mps_dump (s, s->logstr);

  /* Check if there are too large or too small approximations */
  *d_after_f = false;
  for (i = 0; i < s->n; i++)
    if (s->status[i][0] == 'x')
      {
        s->again[i] = false;
        *d_after_f = true;
      }

  /* == 2 ==  Perform max_pack packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE:  call fpolzer\n");
  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* floop: */

      /* mps_fpolzer(s, &nit, &excep); */
      mps_thread_fpolzer (s, &nit, &excep);
      it_pack += nit;

      if (s->DOLOG)
        fprintf (s->logstr, "Packet %d  iterations= %d\n", iter, nit);

      /* perform cluster analysis, shift, restart, update 'statu', and
       * update 'again'
       */
      if (excep)
        {
          for (i = 0; i <= s->n; i++)
            s->oldpunt[i] = s->punt[i];
          oldnclust = s->nclust;

          if (s->DOLOG)
            fprintf (s->logstr, "   FSOLVE: call fcluster\n");
          /* cluster analysis */
          mps_fcluster (s, 2 * s->n);   /* Isolation factor */
          if (oldnclust == s->nclust)
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: cycle\n");
              continue;
            }
          else
            {
              /* modify the vector status and mark also the old
               * clusters with 'C'
               */
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call modify\n");
              mps_fmodify (s, true);

              if (iter == 0)
                for (i = 0; i < s->n; i++)
                  if (s->status[i][0] == 'C')
                    s->status[i][0] = 'c';

              /* If the polynomial is not given in terms of its coeff. then
               * skip the restart stage */
              if (s->data_type[0] != 'u')
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr, "   FSOLVE: call frestart\n");
                  mps_frestart (s);
                }
              /* reset the status vector */
              for (j = 0; j < s->n; j++)
                {
                  if (s->status[j][0] == 'C')
                    s->status[j][0] = 'c';
                  s->again_old[j] = s->again[j];
                }

              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call update\n");
              mps_update (s);

              /* adjust 'again' This is needed since we are
               * between two packets */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->again[i] = false;
              if (s->DOLOG)
                fprintf (s->logstr, "   FSOLVE: call checkstop\n");
              /* Check the stop condition */
              if (mps_check_stop (s))
                return;
            }
        }
      else
        break;
    }

  /* The 'floop' has been completed:
   * If the max number of iteration has been reached then output FAILURE */
  if (iter == s->max_pack)
    {
      mps_dump (s, s->logstr);
      mps_error (s, 1,
                 "Float: reached the maximum number of packet iterations");
    }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double. */

  if (s->DOLOG)
    fprintf (s->logstr, "FLOAT: nit= %d\n", it_pack);

  /* Update */
  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE: call fcluster\n");
  for (i = 0; i <= s->n; i++)
    s->oldpunt[i] = s->punt[i];
  oldnclust = s->nclust;
  mps_fcluster (s, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "   FSOLVE: call modify\n");
  mps_fmodify (s, true);

  /* reset the status vector */
  for (j = 0; j < s->n; j++)
    if (s->status[j][0] == 'C')
      s->status[j][0] = 'c';
}

/**
 * @brief This routine applies <code>nit</code> iterations of 
 * Aberth's method.
 * 
 * The method is applied to
 * the <code>i</code>-th component of the approximations for which 
 * <code>again[i]</code>
 * is <code>true</code>. Set <code>again[i]=false</code> if the <code>i</code>-th 
 * approximation is in
 * the root neighbourhood. Stop if <code>again[i]=false</code> for any <code>i</code>.
 *
 * @param s The <code>mps_status</code> associated with the current computation.
 * @param it Index of the component on which the iteration is needed. 
 * @param excep This variable is set to <code>true</code> if after <code>nit</code> 
 * iterations some approximation is still out of the root neighbourhood.
 */
void
mps_fpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  int i, iter, nzeros;
  cplx_t corr, abcorr;
  double rad1, modcorr;
  mps_monomial_poly * p = s->monomial_poly;

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    return;

  /* Start Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "FPOLZER: starts aberth it\n");

  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter : DO iter=1,nit */

      if (s->DOLOG)
        {
          fprintf (s->logstr, "FPOLZER: iteration %d\n", iter);
          mps_dump (s, s->logstr);
        }

      for (i = 0; i < s->n; i++)
        {                       /* do_index */
          if (s->again[i])
            {
              (*it)++;
              rad1 = s->frad[i];
              if (s->data_type[0] != 'u')
                {
                  mps_fnewton (s, s->n, s->froot[i], &s->frad[i], corr,
                               p->fpc, p->fap, &s->again[i]);
                  if (iter == 0 && !s->again[i] && s->frad[i] > rad1 && rad1
                      != 0)
                    s->frad[i] = rad1;
                                        /***************************************
					 The above condition is needed to cope with the case
					 where at the first iteration the starting point
					 is already in the root neighbourhood and the actually
					 computed radius is too big since the value of the first
					 derivative is too small.
					 In this case the previous radius bound, obtained by
					 means of Rouche' is more reliable and strict
					 **************************************/
                }
              else if (s->fnewton_usr != NULL)
                {
                  (*s->fnewton_usr) (s, s->froot[i], &s->frad[i], corr,
                                     &s->again[i], NULL);
                }
              else
                {
                  mps_fnewton_usr (s, s->froot[i], &s->frad[i], corr,
                                   &s->again[i]);
                }

              if (s->again[i] ||
                  /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
                  s->data_type[0] == 'u' || iter != 0 || s->frad[i] != rad1)
                {
                  mps_faberth (s, i, abcorr);
                  cplx_mul_eq (abcorr, corr);
                  cplx_sub (abcorr, cplx_one, abcorr);
                  cplx_div (abcorr, corr, abcorr);
                  cplx_sub_eq (s->froot[i], abcorr);
                  modcorr = cplx_mod (abcorr);
                  s->frad[i] += modcorr;
                }

              /* check for new approximated roots */
              if (!s->again[i])
                {
                  nzeros++;
                  if (nzeros == s->n)
                    return;
                }

            }
        }
    }
  *excep = true;
}

/**
 * @brief <code>dpe</code> version of <code>fpolzer()</code>.
 */
void
mps_dpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  int iter, i, nzeros;
  rdpe_t rad1, rtmp;
  cdpe_t corr, abcorr;
  mps_monomial_poly * p = s->monomial_poly;

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    return;

  /* Start Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "DPOLZER: starts aberth\n");
  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter: */

      for (i = 0; i < s->n; i++)
        {                       /* do_index: */

          if (s->again[i])
            {
              (*it)++;
              rdpe_set (rad1, s->drad[i]);
              if (s->data_type[0] != 'u')
                {
                  mps_dnewton (s, s->n, s->droot[i], s->drad[i], corr, p->dpc,
                               p->dap, &s->again[i]);
                  if (iter == 0 && !s->again[i] && rdpe_gt (s->drad[i], rad1)
                      && rdpe_ne (rad1, rdpe_zero))
                    rdpe_set (s->drad[i], rad1);
                }
              else if (s->dnewton_usr != NULL)
                {
                  (*s->dnewton_usr) (s, s->droot[i], s->drad[i], corr,
                                     &s->again[i], NULL);
                }
              else
                {
                  mps_dnewton_usr (s, s->droot[i], s->drad[i], corr,
                                   &s->again[i]);
                }

                                /************************************************
				 The above condition is needed to manage with the case where
				 at the first iteration the starting point is already in the
				 root neighbourhood and the actually computed radius is too
				 big since the value of the first derivative is too small.
				 In this case the previous radius bound, obtained by means of
				 Rouche' is more reliable and strict
				 **********************************************/

              if (s->again[i] ||
                  /* the correction is performed only if iter!=1 or rad(i)!=rad1 */
                  s->data_type[0] == 'u' || iter != 0
                  || rdpe_ne (s->drad[i], rad1))
                {
                  mps_daberth (s, i, abcorr);
                  cdpe_mul_eq (abcorr, corr);
                  cdpe_sub (abcorr, cdpe_one, abcorr);
                  cdpe_div (abcorr, corr, abcorr);
                  cdpe_sub_eq (s->droot[i], abcorr);
                  cdpe_mod (rtmp, abcorr);
                  rdpe_add_eq (s->drad[i], rtmp);
                }

              /* check for new approximated roots */
              if (!s->again[i])
                {
                  nzeros++;
                  if (nzeros == s->n)
                    return;
                }

            }
        }
    }
  *excep = true;
}

/**
 * @brief <code>dpe</code> version of <code>fsolve()</code>.
 */
void
mps_dsolve (mps_status * s, mps_boolean d_after_f)
{
  int it_pack, iter, nit, oldnclust, i, j;
  mps_boolean excep;
  rdpe_t dummy;
  mps_monomial_poly * p = s->monomial_poly;

  if (s->DOLOG)
    {
      fprintf (s->logstr, "   DSOLVE: d_after_f= %d\n", d_after_f);
    }

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (d_after_f)
    for (i = 0; i < s->n; i++)
      if (s->status[i][0] == 'x')
        {
          s->again[i] = true;
          rdpe_set_d (s->drad[i], DBL_MAX);
        }
      else
        s->again[i] = false;
  else
    {
      s->nclust = 0;
      for (i = 0; i < s->n; i++)
        {
          s->again[i] = true;
          rdpe_set (s->drad[i], RDPE_MAX);
          cdpe_set (s->droot[i], cdpe_zero);
        }
    }

  /* Choose starting approximations */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: call dstart con again=\n");

  rdpe_set (dummy, rdpe_zero);

  /* Use a custom routine if it is set, otherwise use the default one */
  if (s->dstart_usr)
    (*s->dstart_usr) (s, s->n, 0, dummy, dummy, dummy);
  else
    mps_dstart (s, s->n, 0, dummy, dummy, dummy, p->dap);

  /* Now adjust the status vector */
  if (d_after_f)
    for (i = 0; i < s->n; i++)
      if (s->status[i][0] == 'x')
        s->status[i][0] = 'c';

  if (s->DOLOG)
    mps_dump (s, s->logstr);

  /* == 2 == Perform s->max_pack  packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: call dpolzero\n");

  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* dloop : DO iter=1,s->max_pack */

      /* mps_dpolzer(s, &nit, &excep); */
      mps_thread_dpolzer (s, &nit, &excep);
      it_pack += nit;

      MPS_DEBUG (s, "DPE packet completed in %d iterations", nit);

      if (s->DOLOG)
        fprintf (s->logstr, "Packet %d iterations= %d\n", iter, nit);

      if (excep)
        {
          for (i = 0; i <= s->n; i++)
            s->oldpunt[i] = s->punt[i];
          oldnclust = s->nclust;

          /* cluster analysis */
          if (s->DOLOG)
            fprintf (s->logstr, "   DSOLVE: call dcluster\n");
          mps_dcluster (s, 2 * s->n);   /* Isolation factor */
          if (oldnclust == s->nclust)
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE:  CYCLE\n");
              continue;
            }
          else
            {
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call dmodify\n");
              mps_dmodify (s, true);

              if (iter == 0 && !d_after_f)
                for (i = 0; i < s->n; i++)
                  if (s->status[i][0] == 'C')
                    s->status[i][0] = 'c';

              /* If the polynomial is not given in terms of its
               * coeff. then skip the restart stage */
              if (s->data_type[0] != 'u')
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr, "   DSOLVE: call drestart\n");
                  mps_drestart (s);
                }
              /* reset the status vector */
              for (j = 0; j < s->n; j++)
                if (s->status[j][0] == 'C')
                  s->status[j][0] = 'c';
              for (j = 0; j < s->n; j++)
                s->again_old[j] = s->again[j];

              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call update\n");
              mps_update (s);
              /* adjust 'again'
               * This is needed since we are between two packets
               */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->again[i] = false;
              if (s->DOLOG)
                fprintf (s->logstr, "   DSOLVE: call checkstop\n");
              if (mps_check_stop (s))
                return;
            }
        }
      else
        break;
    }
  if (iter == s->max_pack)
    {
      mps_dump (s, s->logstr);
      mps_error (s, 1,
                 "DPE: reached the maximum number of packet iterations");
    }
  /* Otherwise exit since all the approximations are
   * in the root neighbourhood, except for the ones that cannot be
   * represented as double.
   */

  if (s->DOLOG)
    fprintf (s->logstr, "DPE: nit=%d\n", nit);

  /* Update */
  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: now update: call dcluster\n");
  for (i = 0; i <= s->n; i++)
    s->oldpunt[i] = s->punt[i];
  oldnclust = s->nclust;
  mps_dcluster (s, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "   DSOLVE: now call dmodify\n");
  mps_dmodify (s, true);

  /* reset the status vector */
  for (j = 0; j < s->n; j++)
    if (s->status[j][0] == 'C')
      s->status[j][0] = 'c';
}

/**
 * @brief Multiprecision version of <code>fsolve()</code>.
 */
void
mps_msolve (mps_status * s)
{
  int iter, nit, oldnclust, i, j, it_pack;
  mps_boolean excep;
  int nzc;
  mps_monomial_poly * p;

  /* == 1 == Initialize variables */
  it_pack = 0;

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call restart\n");
  if (s->data_type[0] != 'u')
    mps_mrestart (s);
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call update1\n");
  mps_update (s);
  if (s->DOLOG)
    {
      fprintf (s->logstr, "  MSOLVE: again after update = ");
      for (i = 0; i < s->n; i++)
        fprintf (s->logstr, "%d", s->again[i]);
      fprintf (s->logstr, "\n");
    }

  for (i = 0; i < s->n; i++)
    if (s->again[i])
      s->rootwp[i] = s->mpwp;

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: call checkstop\n");
  if (mps_check_stop (s))
    {
      mps_mmodify (s, true);

      /* reset the status vector */
      for (j = 0; j < s->n; j++)
        if (s->status[j][0] == 'C')
          s->status[j][0] = 'c';

      return;
    }
  nzc = 0;
  if (s->goal[1] == 'a' && s->goal[2] == 'n' && s->goal[3] == 'n')
    for (i = 0; i < s->n; i++)
      if (s->status[i][0] == 'i' || s->status[i][0] == 'a'
          || s->status[i][0] == 'o')
        nzc++;
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: nzc=%d\n", nzc);

  if (nzc == s->n)
    {
      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: call mmodify and return\n");
      mps_mmodify (s, true);

      /* reset the status vector */
      for (j = 0; j < s->n; j++)
        if (s->status[j][0] == 'C')
          s->status[j][0] = 'c';
      return;
    }

  /* Perform s->max_pack  packets of Aberth's iterations */
  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE: Perform packets of Aberth\n");

  for (iter = 0; iter < s->max_pack; iter++)
    {                           /* mloop : DO iter=1,s->max_pack */
      if (s->DOLOG)
        {
          fprintf (s->logstr, "  MSOLVE: packet= %d\n", iter);
          fprintf (s->logstr, "  MSOLVE: again before mpolzer =");
          for (i = 0; i < s->n; i++)
            fprintf (s->logstr, "%d", s->again[i]);
          fprintf (s->logstr, "\n");
          fprintf (s->logstr, "  MSOLVE: call mpolzer\n");
        }
      /* mps_mpolzer(s, &nit, &excep); */
      mps_thread_mpolzer (s, &nit, &excep);

      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: Packet %d: iterations= %d\n", iter,
                 nit);

      it_pack += nit;
      nzc = 0;
      if (s->goal[1] == 'a' && s->goal[2] == 'n' && s->goal[3] == 'n')  /* DARIO APRILE 98 */
        for (i = 0; i < s->n; i++)
          if (s->status[i][0] == 'i' || s->status[i][0] == 'a'
              || s->status[i][0] == 'o')
            nzc++;
      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
      if (nzc == s->n)
        {
          if (s->DOLOG)
            fprintf (s->logstr, "  MSOLVE: call mmodify and return\n");
          mps_mmodify (s, true);

          /* reset the status vector */
          for (j = 0; j < s->n; j++)
            if (s->status[j][0] == 'C')
              s->status[j][0] = 'c';
          if (s->DOLOG)
            {
              fprintf (s->logstr, "           s->status=");
              for (j = 0; j < s->n; j++)
                fprintf (s->logstr, "%3.3s", s->status[j]);
              fprintf (s->logstr, "\n");
            }
          return;
        }

      if (s->DOLOG)
        fprintf (s->logstr, "  MSOLVE: isolated %d roots excep=%d\n", nzc,
                 excep);

      if (excep)
        {
          for (i = 0; i <= s->n; i++)
            s->oldpunt[i] = s->punt[i];
          oldnclust = s->nclust;

          /* cluster analysis */
          if (s->DOLOG)
            fprintf (s->logstr, "  MSOLVE: call mcluster\n");

          mps_mcluster (s, 2 * s->n);   /* Isolation factor */

          s->newtis_old = s->newtis;
          if (s->newtis == 0)
            mps_mnewtis (s);
          if (s->DOLOG)
            fprintf (s->logstr,
                     "  MSOLVE: newtis_old=%d, newtis=%d, oldncl=%d, s->nclust=%d\n",
                     s->newtis_old, s->newtis, oldnclust, s->nclust);

          if (oldnclust == s->nclust && !(s->newtis == 1 && s->newtis_old
                                          == 0))
            /*#             if(&& iter != 0) AGO99 */
            {
              /*#D !newtis */
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: CYCLE\n");
              continue;
            }
          else
            {
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call modify\n");
              mps_mmodify (s, true);
              if (iter == 0)
                /* if first packet: reset the status vector */
                for (j = 0; j < s->n; j++)
                  if (s->status[j][0] == 'C')
                    s->status[j][0] = 'c';
              if (s->DOLOG)
                {
                  fprintf (s->logstr, "  MSOLVE:  status=");
                  for (j = 0; j < s->n; j++)
                    fprintf (s->logstr, "%3.3s", s->status[j]);
                  fprintf (s->logstr, "\n");
                }
              /* If the polynomial is not given in terms of its coeff. then
               * skip the restart stage */
              if (s->data_type[0] != 'u')
                {
                  /* choose new starting approximations only for new clusters */
                  if (s->DOLOG)
                    fprintf (s->logstr,
                             "  MSOLVE: call mrestart for new clusters\n");
                  mps_mrestart (s);
                }
              /* reset the s->status vector */
              for (j = 0; j < s->n; j++)
                {
                  if (s->status[j][0] == 'C')
                    s->status[j][0] = 'c';
                  s->again_old[j] = s->again[j];
                }
              /* update 'again' */
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call update2 : ");
              mps_update (s);
              if (s->DOLOG)
                {
                  fprintf (s->logstr, "  MSOLVE: again = ");
                  for (j = 0; j < s->n; j++)
                    fprintf (s->logstr, "%d", s->again[j]);
                  fprintf (s->logstr, "\n");
                }
              /* adjust 'again'  This is needed since we are between two packets */
              for (i = 0; i < s->n; i++)
                if (!s->again_old[i])
                  s->again[i] = false;
              if (s->DOLOG)
                {
                  fprintf (s->logstr, "  MSOLVE: adjusted again = ");
                  for (j = 0; j < s->n; j++)
                    fprintf (s->logstr, "%d", s->again[j]);
                  fprintf (s->logstr, "\n");
                }
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: call checkstop\n");
              if (mps_check_stop (s))
                {
                  mps_mmodify (s, true);

                  /* reset the s->status vector */
                  for (j = 0; j < s->n; j++)
                    if (s->status[j][0] == 'C')
                      s->status[j][0] = 'c';

                  return;
                }

              nzc = 0;
              if (s->goal[1] == 'a' && s->goal[2] == 'n' && s->goal[3] == 'n')
                for (i = 0; i < s->n; i++)
                  if (s->status[i][0] == 'i' || s->status[i][0] == 'a'
                      || s->status[i][0] == 'o')
                    nzc++;
              if (s->DOLOG)
                fprintf (s->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
              if (nzc == s->n)
                {
                  if (s->DOLOG)
                    fprintf (s->logstr, "  MSOLVE: call mmodify and return");
                  mps_mmodify (s, true);

                  /* reset the s->status vector */
                  for (j = 0; j < s->n; j++)
                    if (s->status[j][0] == 'C')
                      s->status[j][0] = 'c';

                  return;
                }
            }
        }
      else
        break;
    }

  if (iter == s->max_pack)
    {
      mps_error (s, 1, "MP: reached the maximum number of packet iteration");
    }

  if (s->DOLOG)
    {
      fprintf (s->logstr, "  MSOLVE: MP: nit= %d\n", nit);
      fprintf (s->logstr, "  MSOLVE: call mcluster\n");
    }

  mps_mcluster (s, 2 * s->n);   /* Isolation factor */

  if (s->DOLOG)
    fprintf (s->logstr, "  MSOLVE:  call mmodify\n");
  mps_mmodify (s, true);

  for (j = 0; j < s->n; j++)
    if (s->status[j][0] == 'C')
      s->status[j][0] = 'c';

  if (s->DOLOG)
    {
      fprintf (s->logstr, "  MSOLVE: s->status=");
      for (j = 0; j < s->n; j++)
        fprintf (s->logstr, "%3.3s", s->status[j]);
      fprintf (s->logstr, "\n");
      fprintf (s->logstr, "  MSOLVE: call update3 : ");
    }
  mps_update (s);

  if (s->DOLOG)
    {
      for (j = 0; j < s->n; j++)
        fprintf (s->logstr, "%d", s->again[j]);
      fprintf (s->logstr, "\n");
    }
}

/**
 * @brief Multiprecision versione of <code>fpolzer()</code>.
 */
void
mps_mpolzer (mps_status * s, int *it, mps_boolean * excep)
{
  int nzeros, i, j, iter, l;
  tmpc_t corr, abcorr;
  rdpe_t eps, rad1, rtmp;
  cdpe_t ctmp;
  mps_monomial_poly * p = s->monomial_poly;

  tmpc_init2 (abcorr, s->mpwp);
  tmpc_init2 (corr, s->mpwp);

  rdpe_mul_d (eps, s->mp_epsilon, (double) 4 * s->n);

  /* initialize the iteration counter */
  *it = 0;
  *excep = false;

  /* count the number of approximations in the root neighbourhood */
  nzeros = 0;
  for (i = 0; i < s->n; i++)
    if (!s->again[i])
      nzeros++;
  if (nzeros == s->n)
    goto endfun;

  /* Start Aberth's iterations */
  for (iter = 0; iter < s->max_it; iter++)
    {                           /* do_iter: */
      for (j = 0; j < s->nclust; j++)
        {                       /* do_clust: */
          for (i = 0; i < s->punt[j + 1] - s->punt[j]; i++)
            {                   /* do_indice: */
              l = s->clust[s->punt[j] + i];
              if (s->again[l])
                {
                  (*it)++;
                  if (s->data_type[0] != 'u')
                    {
                      /* sparse/dense polynomial */
                      rdpe_set (rad1, s->drad[l]);
                      mps_mnewton (s, s->n, s->mroot[l], s->drad[l], corr,
                                   p->mfpc, p->mfppc, p->dap, p->spar,
                                   &s->again[l], 0);
                      if (iter == 0 && !s->again[l] && rdpe_gt (s->drad[l],
                                                                rad1)
                          && rdpe_ne (rad1, rdpe_zero))
                        rdpe_set (s->drad[l], rad1);

                                                /************************************************
						 The above condition is needed to cope with the case
						 where at the first iteration the starting point is
						 already in the root neighbourhood and the actually
						 computed radius is too big since the value of the
						 first derivative is too small.
						 In this case the previous radius bound, obtained by
						 means of Rouche' is more reliable and strict
						 ***********************************************/
                    }
                  else /* user's polynomial */ if (s->mnewton_usr != NULL)
                    {
                      (*s->mnewton_usr) (s, s->mroot[l], s->drad[l], corr,
                                         &s->again[l], NULL);
                    }
                  else
                    {
                      mps_mnewton_usr (s, s->mroot[l], s->drad[l], corr,
                                       &s->again[l]);
                    }

                  if (s->again[l] ||
                      /* the correction is performed only if iter!=1 or rad[l]!=rad1 */
                      s->data_type[0] == 'u' || iter != 0
                      || rdpe_ne (s->drad[l], rad1))
                    {
                      mps_maberth_s (s, l, j, abcorr);
                      mpc_mul_eq (abcorr, corr);
                      mpc_neg_eq (abcorr);
                      mpc_add_eq_ui (abcorr, 1, 0);
                      mpc_div (abcorr, corr, abcorr);
                      mpc_sub_eq (s->mroot[l], abcorr);
                      mpc_get_cdpe (ctmp, abcorr);
                      cdpe_mod (rtmp, ctmp);
                      rdpe_add_eq (s->drad[l], rtmp);
                    }

                  /* check for new approximated roots */
                  if (!s->again[l])
                    {
                      nzeros++;
                      if (nzeros == s->n)
                        goto endfun;
                    }

                }
            }
        }
      if (nzeros == s->n)
        goto endfun;
    }
  *excep = true;

endfun:                        /* free local MP variables */
  tmpc_clear (corr);
  tmpc_clear (abcorr);
}
