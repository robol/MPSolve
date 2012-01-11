/************************************************************
 **                                                        **
 **             __  __ ___  ___      _                     **
 **            |  \/  | _ \/ __| ___| |_ _____             **
 **            | |\/| |  _/\__ \/ _ \ \ V / -_)            **
 **            |_|  |_|_|  |___/\___/_|\_/\___|            **
 **                                                        **
 **       Multiprecision Polynomial Solver (MPSolve)       **
 **                 Version 2.9, April 2011                **
 **                                                        **
 **                      Written by                        **
 **                                                        **
 **     Dario Andrea Bini       <bini@dm.unipi.it>         **
 **     Giuseppe Fiorentino     <fiorent@dm.unipi.it>      **
 **     Leonardo Robol          <robol@mail.dm.unipi.it>   **
 **                                                        **
 **           (C) 2011, Dipartimento di Matematica         **
 ***********************************************************/

#include <mps/mps.h>
#include <math.h>

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
  int l, i;
  rdpe_t rtmp;
  double eps_out = rdpe_get_d (s->eps_out);

  /* Set isolation factor */
  /* nf = 2 * s->n; */

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->root_status[i] == MPS_ROOT_STATUS_CLUSTERED)
	  s->root_status[i] = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
	{
	  /* Check if the root is already approximated; if that's not the
	   * case set it at least as isolated. */
	  if (s->root_status[l] != MPS_ROOT_STATUS_APPROXIMATED)
	    {
	      s->root_status[root->k] = MPS_ROOT_STATUS_ISOLATED;
	      /* s->root_inclusion[root->k] = MPS_ROOT_INCLUSION_IN; */

	      /* Check if we need to mark this root as approximated */
	      if (s->frad[l] < cplx_mod (s->froot[l]) * eps_out)
		s->root_status[root->k] = MPS_ROOT_STATUS_APPROXIMATED;
	    }

	  /* Grab the next cluster and continue scanning */
	  c_item = c_item->next;
	  continue;
	}

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
	{
	  l = root->k;

	  /* If track_new_cluster is false then we may directly set here the
	   * approximation status of the roots. */
	  if (!track_new_cluster)
	    {	  
	      s->root_status[l] = MPS_ROOT_STATUS_CLUSTERED;
	      rdpe_set_d (rtmp, s->frad[l] / cplx_mod (s->froot[l]));
	      if (rdpe_le (rtmp, s->eps_out))
		s->root_status[l] = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
	    }

	  root = root->next;
	}

      /* TODO: Implement checking of the zone where the roots are. */

      c_item = c_item->next;
    }

  mps_fupdate_inclusions (s);
}

/* void */
/* mps_fmodify2 (mps_status * s, mps_boolean track_new_cluster) */
/* { */
/*   int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, nf, j2, l2; */
/*   double sr, tmpr, afri, sep1; */
/*   cplx_t sc; */
/*   rdpe_t rtmp; */
/*   mps_boolean tcr, tcr1; */

/*   /\* ==1== Change into 'C' the components of status for old clusters *\/ */
/*   nf = 2 * s->n;                /\* Isolation factor *\/ */

/*   /\* Mark old cluster with 'C' if requested *\/ */
/*   if (track_new_cluster) */
/*     { */
/*       for (i = 0; i < s->n; i++) */
/* 	if (s->status[i][0] == 'c') */
/* 	  s->status[i][0] = 'C'; */
/*     } */

/*   i_old = 0; */
/*   i_new = 0; */
/*   for (i = 1; i <= s->nclust && i_new < s->n; i++) */
/*     {                           /\*  loop1: DO i=1, nclust *\/ */
/*       if (s->oldpunt[i_old + 1] == s->punt[i_new + 1]) */
/*         { */
/*           i_old++; */
/*           i_new++; */
/*           continue; */
/*         } */
/*       else */
/*         { */
/*           for (j = i_new + 1; j < s->nclust; j++) */
/*             {                   /\* loop2: DO j=i_new+1, nclust *\/ */
/*               if (s->oldpunt[i_old + 1] != s->punt[j + 1]) */
/*                 continue; */
/*               else */
/*                 { */
/*                   nnewclust = j - i_new + 1;    /\* scan each new cluster *\/ */
/*                   for (k = 0; k < nnewclust; k++) */
/*                     {           /\* loop3: DO k=1, nnewclust *\/ */
/*                       i1 = i_new + k; */
/*                                                 /\********************************* */
/* 						 scan the entries of each new cluster */
/* 						 set status[l][0]='i' if the cluster has multip=1 */
/* 						 and mark with 'c' */
/* 						 the ones which are different from 'i' */
/* 						 **********************************\/ */
/*                       if (s->punt[i1 + 1] - s->punt[i1] == 1 */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'x' */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'f') */
/*                         s->status[s->clust[s->punt[i1]]][0] = 'i'; */
/*                       for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++) */
/*                         {       /\* loop4: *\/ */
/*                           ip1 = s->clust[s->punt[i1] + l]; */
/*                           if (s->status[ip1][0] != 'i' && s->status[ip1][0] */
/*                               != 'x' && s->status[ip1][0] != 'f' */
/*                               && s->status[ip1][0] != 'a' */
/*                               && s->status[ip1][0] != 'o') */
/*                             s->status[ip1][0] = 'c'; */
/*                         } */
/*                     } */
/*                   i_new = j + 1; */
/*                   i_old++; */
/*                   break; */
/*                 }               /\* else *\/ */
/*             }                   /\* for *\/ */
/*         }                       /\* else *\/ */
/*     }                           /\* for *\/ */

/*         /\*=2== Scan all the clusters *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\* scan : DO i=1,s->nclust *\/ */
/*       /\* check isolation/approximation *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1 */
/*           && s->status[s->clust[s->punt[i]]][0] != 'x' */
/*           && s->status[s->clust[s->punt[i]]][0] != 'f') */
/*         { */
/*           s->status[s->clust[s->punt[i]]][0] = 'i'; */
/*           tmpr = cplx_mod (s->froot[s->clust[s->punt[i]]]); */
/*           tmpr = s->frad[s->clust[s->punt[i]]] / tmpr; */
/*           tmpr = log (tmpr); */
/*           if (tmpr < -s->output_config->prec * LOG2) */
/*             s->status[s->clust[s->punt[i]]][0] = 'a'; */
/*         } */
/*       /\* Scan inside the cluster *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan_in *\/ */
/*           l = s->clust[s->punt[i] + j]; */
/*                         /\************************************************** */
/* 			 first check for inside/outside unit disk in the case where */
/* 			 there are very large and/or very small roots (statu='x') */
/* 			 and for counting only */
/* 			 *************************************************\/ */
/*           afri = cplx_mod (s->froot[i]); */
	  
/* 	  /\* Check if the root, even if clustered, is approximated *\/ */
/* 	  if (s->algorithm == MPS_ALGORITHM_SECULAR_GA) */
/* 	    { */
/* 	      rdpe_set_d (rtmp, s->frad[l] / cplx_mod (s->froot[l])); */
/* 	      if (rdpe_le (rtmp, s->eps_out))      */
/* 		s->status[l][0] = 'o';      */
/* 	    } */

/*           if (s->status[l][0] == 'x' && s->goal[0] == 'c') */
/*             { */
/*               if ((s->goal[1] == 'i' && afri < 1) || (s->goal[1] == 'o' */
/*                                                       && afri > 1)) */
/*                 s->status[l][2] = 'i'; */
/*               if ((s->goal[1] == 'i' && afri > 1) || (s->goal[1] == 'o' */
/*                                                       && afri < 1)) */
/*                 s->status[l][2] = 'o'; */
/*             } */
/*           /\* Now check the standard cases *\/ */
/*           if (s->status[l][0] == 'x') */
/*             continue; */
/*           if ((s->status[l][0] == 'c' || s->status[l][0] == 'i' */
/*                || s->status[l][0] == 'C' || s->status[l][0] == 'o') */
/*               && s->status[l][2] == 'u') */
/*             { */
/*               /\* Check if the approximation is inside/outside the set *\/ */
/*               switch (s->goal[1]) */
/*                 { */

/*                 case 'a':      /\* all *\/ */
/*                   s->status[l][2] = 'i'; */
/*                   break; */

/*                 case 'i':      /\* inside unit circle *\/ */
/*                   if (!mps_ftouchunit (s, nf, l)) */
/*                     { */
/*                       if (cplx_mod (s->froot[l]) < 1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'o':      /\* outside unit circle *\/ */
/*                   if (!mps_ftouchunit (s, nf, l)) */
/*                     { */
/*                       if (cplx_mod (s->froot[l]) > 1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'l':      /\* left half plane  *\/ */
/*                   if (!mps_ftouchimag (s, nf, l)) */
/*                     { */
/*                       if (cplx_Re (s->froot[l]) < 0) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'r':      /\* right half plane *\/ */
/*                   if (!mps_ftouchimag (s, nf, l)) */
/*                     { */
/*                       if (cplx_Re (s->froot[l]) > 0) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'u':      /\* upper half plane *\/ */
/*                   if (!mps_ftouchreal (s, nf, l)) */
/*                     { */
/*                       if (cplx_Im (s->froot[l]) > 0) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'd':      /\* lower half plane *\/ */
/*                   if (!mps_ftouchreal (s, nf, l)) */
/*                     { */
/*                       if (cplx_Im (s->froot[l]) < 0) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'R':      /\* Real line  NEW *\/ */
/*                   if (s->status[l][1] != 'w') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_ftouchreal (s, 1, l)) */
/*                         { */
/*                           if (s->data_type[1] == 'r') */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'R'; */
/*                             } */
/*                           else */
/*                             { */
/*                               /\* fsrad(i, sc, &sr); DARIO *\/ */
/*                               sr = s->frad[l]; */
/*                               if (log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                                 { */
/*                                   s->status[l][2] = 'i'; */
/*                                   s->status[l][1] = 'R'; */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   s->status[l][2] = 'u'; */
/*                                   s->status[l][1] = 'w'; */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch real *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'r'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_ftouchreal (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_ftouchreal (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\*  mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano R *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_fsrad (s, i, sc, &sr); */
/*                               sep1 = s->sep; */
/*                               if (s->data_type[1] == 'c') */
/*                                 sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'R'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano R *\/ */

/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][2] = 'o'; */
/*                               s->status[l2][1] = 'r'; */
/*                             } */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'I':      /\* Imaginary line  NEW *\/ */
/*                   if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_ftouchimag (s, nf, l)) */
/*                         { */
/*                           /\* fsrad(i, sc, &sr); DARIO *\/ */
/*                           sr = s->frad[l]; */
/*                           if (log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'I'; */
/*                             } */
/*                           else */
/*                             { */
/*                               s->status[l][2] = 'u'; */
/*                               s->status[l][1] = 'w'; */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch imag *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_ftouchimag (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_ftouchimag (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\* mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano I *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_fsrad (s, i, sc, &sr); */
/*                               sep1 = s->sep; */
/*                               sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'I'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano I *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'S':      /\* Set provided by the user *\/ */
/*                   mps_error (s, 1, "User set not implemented yet"); */
/*                   break; */

/*                 default: */
/*                   mps_error (s, 1, "Mistake in goal"); */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*       /\* If some cluster still contains an uncertain disk then set */
/*        * all the disks uncertain *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         { */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][2] == 'u') */
/*             { */
/*               for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                 { */
/*                   l1 = s->clust[s->punt[i] + j1]; */
/*                   s->status[l1][2] = 'u'; */
/*                 } */
/*               break; */
/*             } */
/*         } */
/*     scan:;                     /\* scan the next component *\/ */
/*     } */

/*         /\*==3==  now check the options *\/ */

/*   /\* Option multiplicity *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\* scan1 *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1) */
/*         continue; */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan1_in *\/ */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][0] == 'x' || s->status[l][0] == 'f' */
/*               || s->status[l][0] == 'm') */
/*             goto scan1; */
/*           if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /\* NEW *\/ */
/*                                     s->status[l][0] == 'C')) */
/*             {                   /\* multiplicity on *\/ */
/*               if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /\* float coeff. *\/ */
/*                 mps_error (s, 1, */
/*                            "Fatal: Float coefficients - impossible to detect multiplicity"); */

/*               /\* compute super center and super radius *\/ */
/*               mps_fsrad (s, i, sc, &sr); */

/*               if (log (sr) < s->sep) */
/*                 for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                   { */
/*                     l1 = s->clust[s->punt[i] + j1];     /\* NEW j-> j1 *\/ */
/*                     s->status[l1][0] = 'm'; */
/*                   } */
/*               goto scan1; */
/*             } */
/*         } */
/*     scan1:;                    /\* scan next component *\/ */
/*     } */

/*   /\* Option Real check *\/ */
/*   if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan2 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan2_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'i') */
/*                 continue; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan2; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_ftouchreal (s, nf, l)) */
/*                     { */
/*                       if (s->data_type[1] == 'r') */
/*                         s->status[l][1] = 'R'; */
/*                       else */
/*                         { */
/*                           /\* fsrad(i, sc, &sr); DARIO *\/ */
/*                           sr = s->frad[l]; */
/*                           if (log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             s->status[l][1] = 'R'; */
/*                           else */
/*                             s->status[l][1] = 'w'; */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* do not touch real *\/ */
/*                     s->status[l][1] = 'r'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_ftouchreal (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_ftouchreal (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan2_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano R *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_fsrad (s, i, sc, &sr); */
/*                           sep1 = s->sep; */
/*                           if (s->data_type[1] == 'c') */
/*                             sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'R'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano R *\/ */
/*                     s->status[l][1] = 'r'; */
/*                 } */
/*             scan2_in:; */
/*             } */
/*         scan2:; */
/*         } */
/*     } */
/*   /\* Option Imaginary check *\/ */
/*   if (s->goal[3] == 'i' || s->goal[3] == 'b') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan3 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan3_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan3; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                 continue; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_ftouchimag (s, nf, l)) */
/*                     { */
/*                       /\* fsrad(i, sc, &sr); DARIO *\/ */
/*                       sr = s->frad[l]; */
/*                       if (log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                         s->status[l][1] = 'I'; */
/*                       else */
/*                         s->status[l][1] = 'w'; */
/*                     } */
/*                   else */
/*                     /\* do not touch imag *\/ */
/*                     s->status[l][1] = 'i'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_ftouchimag (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_ftouchimag (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan3_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano I *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_fsrad (s, i, sc, &sr); */
/*                           sep1 = s->sep; */
/*                           sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'I'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano I *\/ */
/*                     s->status[l][1] = 'i'; */
/*                 } */
/*             scan3_in:; */
/*             } */
/*         scan3:; */
/*         } */
/*     } */
/* } */



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
  int l, i;
  rdpe_t tmpr, tmpr2;

  /* Set isolation factor */
  /* nf = 2 * s->n; */

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->root_status[i] == MPS_ROOT_STATUS_CLUSTERED)
	  s->root_status[i] = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
	{
	  /* Check if the root is already approximated; if that's not the
	   * case set it at least as isolated. */
	  if (s->root_status[l] != MPS_ROOT_STATUS_APPROXIMATED)
	    s->root_status[root->k] = MPS_ROOT_STATUS_ISOLATED;

	  /* Grab the next cluster and continue scanning */
	  c_item = c_item->next;
	  continue;
	}

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
	{
	  l = root->k;

	  /* If track_new_cluster is false then we may directly set here the
	   * approximation status of the roots. */
	  if (!track_new_cluster)
	    {
	      s->root_status[l] = MPS_ROOT_STATUS_CLUSTERED;
	      rdpe_set (tmpr, s->drad[l]);
	      cdpe_mod (tmpr2, s->droot[l]);
	      rdpe_div_eq (tmpr, tmpr2);
	      if (rdpe_le (tmpr, s->eps_out)) 
		s->root_status[l] = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
	    }

	  root = root->next;
	}

      /* TODO: Implement checking of the zone where the roots are. */
      c_item = c_item->next;
    }

  mps_dupdate_inclusions (s);
}


/* void */
/* mps_dmodify (mps_status * s, mps_boolean track_new_cluster) */
/* { */
/*   int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, j2, l2, nf; */
/*   double rtmp, sep1; */
/*   rdpe_t sr, tmpr, tmpr2; */
/*   cdpe_t sc; */
/*   mps_boolean tcr, tcr1; */

/*   /\* ==1==  Change into 'C' the components of status for old clusters *\/ */
/*   nf = 2 * s->n;                /\* Isolation factor *\/ */

/*   if (track_new_cluster) */
/*     { */
/*       for (i = 0; i < s->n; i++) */
/* 	if (s->status[i][0] == 'c') */
/* 	  s->status[i][0] = 'C'; */
/*     } */

/*   i_old = 0; */
/*   i_new = 0; */
/*   for (i = 1; i <= s->nclust && i_new < s->n; i++) */
/*     {                           /\*  loop1: DO i=1, s->nclust *\/ */
/*       if (s->oldpunt[i_old + 1] == s->punt[i_new + 1]) */
/*         { */
/*           i_old++; */
/*           i_new++; */
/*           continue; */
/*         } */
/*       else */
/*         { */
/*           for (j = i_new + 1; j < s->nclust; j++) */
/*             {                   /\* loop2:  *\/ */
/*               if (s->oldpunt[i_old + 1] != s->punt[j + 1]) */
/*                 continue; */
/*               else */
/*                 { */
/*                   nnewclust = j - i_new + 1;    /\*  scan each new cluster *\/ */
/*                   for (k = 0; k < nnewclust; k++) */
/*                     {           /\* loop3:  *\/ */
/*                       i1 = i_new + k; */

/* 		      /\* Check if the root, even if clustered, is approximated *\/ */
/* 		      if (s->algorithm == MPS_ALGORITHM_SECULAR_GA) */
/* 			{ */
/* 			  rdpe_set (tmpr, s->drad[i1]); */
/* 			  cdpe_mod (tmpr2, s->droot[i1]); */
/* 			  rdpe_div_eq (tmpr, tmpr2); */
/* 			  if (rdpe_le (tmpr, s->eps_out))  */
/* 			    s->status[i1][0] = 'o';  */
/* 			} */
			  
/*                       /\* scan the entries of each new cluster */
/*                        * set status[l][0]='i' if the cluster has multip=1 */
/*                        * and mark with 'c' those which are different from 'i' */
/*                        *\/ */
/*                       if (s->punt[i1 + 1] - s->punt[i1] == 1 */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'x' */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'f') */
/*                         s->status[s->clust[s->punt[i1]]][0] = 'i'; */
/*                       for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++) */
/*                         {       /\* loop4: *\/ */
/*                           ip1 = s->clust[s->punt[i1] + l]; */
/*                           if (s->status[ip1][0] != 'i' && s->status[ip1][0] */
/*                               != 'x' && s->status[ip1][0] != 'f' */
/*                               && s->status[ip1][0] != 'a' */
/*                               && s->status[ip1][0] != 'o') */
/*                             s->status[ip1][0] = 'c'; */
/*                         } */
/*                     } */
/*                   i_new = j + 1; */
/*                   i_old++; */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*     } */

/*         /\*=2== Scan all the clusters *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\* scan: *\/ */
/*       /\* check isolation/approximation *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1 */
/*           && s->status[s->clust[s->punt[i]]][0] != 'x' */
/*           && s->status[s->clust[s->punt[i]]][0] != 'f') */
/*         { */
/*           s->status[s->clust[s->punt[i]]][0] = 'i'; */
/*           cdpe_mod (tmpr, s->droot[s->clust[s->punt[i]]]); */
/*           rdpe_div (tmpr, s->drad[s->clust[s->punt[i]]], tmpr); */
/*           rtmp = rdpe_log (tmpr); */
/*           if (rtmp < -s->output_config->prec * LOG2) */
/*             s->status[s->clust[s->punt[i]]][0] = 'a'; */
/*         } */
/*       /\* Scan inside the cluster *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan_in: *\/ */
/*           l = s->clust[s->punt[i] + j]; */

/*           /\* Now check the standard cases *\/ */
/*           if (s->status[l][0] == 'x') */
/*             continue; */
/*           if ((s->status[l][0] == 'c' || s->status[l][0] == 'i' */
/*                || s->status[l][0] == 'C' || s->status[l][0] == 'o') */
/*               && s->status[l][2] == 'u') */
/*             { */
/*               /\* Check if the approximation is inside/outside the set *\/ */
/*               switch (s->goal[1]) */
/*                 { */

/*                 case 'a':      /\* all *\/ */
/*                   s->status[l][2] = 'i'; */
/*                   break; */

/*                 case 'i':      /\* inside unit circle *\/ */
/*                   if (!mps_dtouchunit (s, nf, l)) */
/*                     { */
/*                       cdpe_mod (tmpr, s->droot[l]); */
/*                       if (rdpe_lt (tmpr, rdpe_one)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'o':      /\* outside unit circle *\/ */
/*                   if (!mps_dtouchunit (s, nf, l)) */
/*                     { */
/*                       cdpe_mod (tmpr, s->droot[l]); */
/*                       if (rdpe_gt (tmpr, rdpe_one)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'l':      /\* left half plane  *\/ */
/*                   if (!mps_dtouchimag (s, nf, l)) */
/*                     { */
/*                       rdpe_set (tmpr, cdpe_Re (s->droot[l])); */
/*                       if (rdpe_lt (tmpr, rdpe_zero)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'r':      /\* right half plane *\/ */
/*                   if (!mps_dtouchimag (s, nf, l)) */
/*                     { */
/*                       rdpe_set (tmpr, cdpe_Re (s->droot[l])); */
/*                       if (rdpe_gt (tmpr, rdpe_zero)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'u':      /\* upper half plane *\/ */
/*                   if (!mps_dtouchreal (s, nf, l)) */
/*                     { */
/*                       rdpe_set (tmpr, cdpe_Im (s->droot[l])); */
/*                       if (rdpe_gt (tmpr, rdpe_zero)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'd':      /\* lower half plane *\/ */
/*                   if (!mps_dtouchreal (s, nf, l)) */
/*                     { */
/*                       rdpe_set (tmpr, cdpe_Im (s->droot[l])); */
/*                       if (rdpe_lt (tmpr, rdpe_zero)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'R':      /\* Real line  NEW *\/ */
/*                   if (s->status[l][1] != 'w') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_dtouchreal (s, 1, l)) */
/*                         { */
/*                           if (s->data_type[1] == 'r') */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'R'; */
/*                             } */
/*                           else */
/*                             { */
/*                               /\* dsrad(i, sc, sr); DARIO *\/ */
/*                               rdpe_set (sr, s->drad[l]); */
/*                               if (rdpe_log (sr) < s->sep - s->n */
/*                                   * s->lmax_coeff) */
/*                                 { */
/*                                   s->status[l][2] = 'i'; */
/*                                   s->status[l][1] = 'R'; */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   s->status[l][2] = 'u'; */
/*                                   s->status[l][1] = 'w'; */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch real *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'r'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_dtouchreal (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_dtouchreal (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\*  mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano R *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_dsrad (s, i, sc, sr); */
/*                               sep1 = s->sep; */
/*                               if (s->data_type[1] == 'c') */
/*                                 sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (rdpe_log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'R'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano R *\/ */

/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][2] = 'o'; */
/*                               s->status[l2][1] = 'r'; */
/*                             } */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'I':      /\* Imaginary line  NEW *\/ */
/*                   if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_dtouchimag (s, nf, l)) */
/*                         { */
/*                           /\* dsrad(i, sc, sr); *\/ */
/*                           rdpe_set (sr, s->drad[l]); */
/*                           if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'I'; */
/*                             } */
/*                           else */
/*                             { */
/*                               s->status[l][2] = 'u'; */
/*                               s->status[l][1] = 'w'; */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch imag *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_dtouchimag (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_dtouchimag (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\* mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano I *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_dsrad (s, i, sc, sr); */
/*                               sep1 = s->sep; */
/*                               sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (rdpe_log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'I'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano I *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'S':      /\* Set provided by the user *\/ */
/*                   mps_error (s, 1, "Custom region not implemented yet"); */
/*                   break; */
/*                 default: */
/*                   mps_error (s, 1, "mistake in goal"); */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*       /\* If some cluster still contains an uncertain disk then set */
/*        * all the disks uncertain */
/*        *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         { */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][2] == 'u') */
/*             { */
/*               for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                 { */
/*                   l1 = s->clust[s->punt[i] + j1]; */
/*                   s->status[l1][2] = 'u'; */
/*                 } */
/*               break; */
/*             } */
/*         } */
/*     scan:; */
/*     } */

/*         /\*==3==  now check the options *\/ */

/*   /\* Option multiplicity *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\* scan1 *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1) */
/*         continue; */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan1_in *\/ */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][0] == 'x' || s->status[l][0] == 'f' */
/*               || s->status[l][0] == 'm') */
/*             goto scan1; */
/*           if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /\* NEW *\/ */
/*                                     s->status[l][0] == 'C')) */
/*             {                   /\* multiplicity on *\/ */

/*               if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /\* float coeff. *\/ */
/*                 mps_error (s, 1, */
/*                            "Fatal: Float coefficients - impossible to detect multiplicity"); */

/*               /\* compute super center and super radius *\/ */
/*               mps_dsrad (s, i, sc, sr); */

/*               if (rdpe_log (sr) < s->sep) */
/*                 for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                   { */
/*                     l1 = s->clust[s->punt[i] + j1];     /\* NEW j-> j1 *\/ */
/*                     s->status[l1][0] = 'm'; */
/*                   } */
/*               goto scan1; */
/*             } */
/*         } */
/*     scan1:;                    /\* scan next component *\/ */
/*     } */

/*   /\* Option Real check *\/ */
/*   if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan2 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan2_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'i') */
/*                 continue; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan2; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_dtouchreal (s, nf, l)) */
/*                     { */
/*                       if (s->data_type[1] == 'r') */
/*                         s->status[l][1] = 'R'; */
/*                       else */
/*                         { */
/*                           /\* dsrad(i, sc, sr); DARIO *\/ */
/*                           rdpe_set (sr, s->drad[l]); */
/*                           if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             s->status[l][1] = 'R'; */
/*                           else */
/*                             s->status[l][1] = 'w'; */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* do not touch real *\/ */
/*                     s->status[l][1] = 'r'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_dtouchreal (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_dtouchreal (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan2_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano R *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_dsrad (s, i, sc, sr); */
/*                           sep1 = s->sep; */
/*                           if (s->data_type[1] == 'c') */
/*                             sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (rdpe_log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'R'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano R *\/ */
/*                     s->status[l][1] = 'r'; */
/*                 } */
/*             scan2_in:; */
/*             } */
/*         scan2:; */
/*         } */
/*     } */

/*   /\* Option Imaginary check *\/ */
/*   if (s->goal[3] == 'i' || s->goal[3] == 'b') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan3 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan3_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan3; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                 continue; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_dtouchimag (s, nf, l)) */
/*                     { */
/*                       /\* dsrad(i, sc, sr); DARIO *\/ */
/*                       rdpe_set (sr, s->drad[l]); */
/*                       if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                         s->status[l][1] = 'I'; */
/*                       else */
/*                         s->status[l][1] = 'w'; */
/*                     } */
/*                   else */
/*                     /\* do not touch imag *\/ */
/*                     s->status[l][1] = 'i'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_dtouchimag (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_dtouchimag (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan3_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano I *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_dsrad (s, i, sc, sr); */
/*                           sep1 = s->sep; */
/*                           sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (rdpe_log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'I'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano I *\/ */
/*                     s->status[l][1] = 'i'; */
/*                 } */
/*             scan3_in:; */
/*             } */
/*         scan3:; */
/*         } */
/*     } */
/* } */


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
  int l, i;
  rdpe_t tmpr, tmpr2;
  cdpe_t cdtmp;

  /* If tracking of new cluster is enabled we need to mark old clusters
   * with 'C' to distinguish them from the new ones. */
  if (track_new_cluster)
    {
      for (i = 0; i < s->n; i++)
	if (s->root_status[i] == MPS_ROOT_STATUS_CLUSTERED)
	  s->root_status[i] = MPS_ROOT_STATUS_NEW_CLUSTERED;
    }

  /* Iterate over the cluster to update the status of the roots */
  mps_cluster_item * c_item;
  mps_cluster * cluster;

  c_item = s->clusterization->first;
  while (c_item != NULL)
    {
      mps_root * root;
      cluster = c_item->cluster;
      
      /* Pick the first root in the cluster */
      root = cluster->first;
      l = root->k;

      /* Check if this is an isolated cluster */
      if (cluster->n == 1)
	{
	  /* Check if the root is already approximated; if that's not the
	   * case set it at least as isolated. */
	  if (s->root_status[l] != MPS_ROOT_STATUS_APPROXIMATED)
	    s->root_status[root->k] = MPS_ROOT_STATUS_ISOLATED;

	  /* Grab the next cluster and continue scanning */
	  c_item = c_item->next;
	  continue;
	}

      /* If it's not the case scan the roots in the cluster and set them
       * to 'c'. */
      while (root != NULL)
	{
	  l = root->k;

	  /* If track_new_cluster is false then we may directly set here the
	   * approximation status of the roots. */
	  if (!track_new_cluster)
	    {
	      s->root_status[l] = MPS_ROOT_STATUS_CLUSTERED;
	      rdpe_set (tmpr, s->drad[l]);
	      mpc_get_cdpe (cdtmp, s->mroot[l]);
	      cdpe_mod (tmpr2, cdtmp);
	      rdpe_div_eq (tmpr, tmpr2);
	      if (rdpe_le (tmpr, s->eps_out)) 
		s->root_status[l] = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER; 
	    }

	  root = root->next;
	}

      /* TODO: Implement checking of the zone where the roots are. */
      c_item = c_item->next;
    }

  mps_mupdate_inclusions (s);
}


/* void */
/* mps_mmodify (mps_status * s, mps_boolean track_new_cluster) */
/* { */
/*   int i, j, l, k, nnewclust, i_new, i_old, ip1, i1, l1, j1, nf, j2, l2; */
/*   double rtmp, sep1; */
/*   rdpe_t sr, tmpr, tmpr2; */
/*   cdpe_t tmpc; */
/*   mps_boolean tcr, tcr1; */
/*   mpf_t tmpf; */
/*   mpc_t sc; */

/*   mpc_init2 (sc, s->mpwp); */
/*   mpf_init2 (tmpf, s->mpwp); */

/*   /\* ==1==  Change into 'C' the components of status for old clusters *\/ */
/*   nf = 2 * s->n;                /\* Isolation factor *\/ */

/*   if (track_new_cluster) */
/*     { */
/*       for (i = 0; i < s->n; i++) */
/* 	if (s->status[i][0] == 'c') */
/* 	  s->status[i][0] = 'C'; */
/*     } */

/*   i_old = 0; */
/*   i_new = 0; */
/*   for (i = 1; i <= s->nclust && i_new < s->n; i++) */
/*     {                           /\* loop1: *\/ */

/*       if (s->oldpunt[i_old + 1] ==  */
/* 	  s->punt[i_new + 1]) */
/*         { */
/*           i_old++; */
/*           i_new++; */
/*           continue; */
/*         } */
/*       else */
/*         { */
/*           for (j = i_new + 1; j < s->nclust; j++) */
/*             {                   /\* loop2: *\/ */
/*               if (s->oldpunt[i_old + 1] != s->punt[j + 1]) */
/*                 continue; */
/*               else */
/*                 { */
/*                   nnewclust = j - i_new + 1;    /\* scan each new cluster *\/ */
/*                   for (k = 0; k < nnewclust; k++) */
/*                     {           /\* loop3: *\/ */
/*                       i1 = i_new + k; */

/* 		      /\* Check if the root, even if clustered, is approximated *\/ */
/* 		      if (s->algorithm == MPS_ALGORITHM_SECULAR_GA) */
/* 			{ */
/* 			  rdpe_set (tmpr, s->drad[i1]); */
/* 			  cdpe_mod (tmpr2, s->droot[i1]); */
/* 			  rdpe_div_eq (tmpr, tmpr2); */
/* 			  if (rdpe_le (tmpr, s->eps_out))  */
/* 			    s->status[i1][0] = 'o';  */
/* 			} */

/*                                                 /\***************************************** */
/* 						 scan the entries of each new cluster set */
/* 						 status[l][0]='i' if the cluster has multip=1 and */
/* 						 mark with 'c' the ones which are different from 'i' */
/* 						 *****************************************\/ */
/*                       if (s->punt[i1 + 1] - s->punt[i1] == 1 */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'x' */
/*                           && s->status[s->clust[s->punt[i1]]][0] != 'f') */
/*                         s->status[s->clust[s->punt[i1]]][0] = 'i'; */
/*                       for (l = 0; l < s->punt[i1 + 1] - s->punt[i1]; l++) */
/*                         {       /\* loop4: *\/ */
/*                           ip1 = s->clust[s->punt[i1] + l]; */
/*                           if (s->status[ip1][0] != 'i' && s->status[ip1][0] */
/*                               != 'x' && s->status[ip1][0] != 'f' */
/*                               && s->status[ip1][0] != 'a' */
/*                               && s->status[ip1][0] != 'o') */
/*                             s->status[ip1][0] = 'c'; */
/*                         } */
/*                     } */
/*                   i_new = j + 1; */
/*                   i_old++; */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*     } */

/*         /\*=2== Scan all the clusters *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\*  scan : DO i=1,s->nclust *\/ */
/*       /\* check isolation/approximation *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1 */
/*           && s->status[s->clust[s->punt[i]]][0] != 'x' */
/*           && s->status[s->clust[s->punt[i]]][0] != 'f') */
/*         { */
/*           s->status[s->clust[s->punt[i]]][0] = 'i'; */
/*           mpc_get_cdpe (tmpc, s->mroot[s->clust[s->punt[i]]]); */
/*           cdpe_mod (tmpr, tmpc); */
/*           rdpe_div (tmpr, s->drad[s->clust[s->punt[i]]], tmpr); */
/*           rtmp = rdpe_log (tmpr); */
/*           if (rtmp < -s->output_config->prec * LOG2) */
/*             { */
/*               s->status[s->clust[s->punt[i]]][0] = 'a'; */
/*             } */
/*         } */
/*       /\* Scan inside the cluster *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan_in: *\/ */
/*           l = s->clust[s->punt[i] + j]; */

/*           /\* Now check the standard cases *\/ */
/*           if (s->status[l][0] == 'x') */
/*             continue; */
/*           if ((s->status[l][0] == 'c' || s->status[l][0] == 'i' */
/*                || s->status[l][0] == 'a' || s->status[l][0] == 'C' */
/*                || s->status[l][0] == 'o') && s->status[l][2] == 'u') */
/*             { */
/*               /\* Check if the approximation is inside/outside the set *\/ */
/*               switch (s->goal[1]) */
/*                 { */

/*                 case 'a':      /\* all *\/ */
/*                   s->status[l][2] = 'i'; */
/*                   break; */

/*                 case 'i':      /\* inside unit circle *\/ */
/*                   if (!mps_mtouchunit (s, nf, l)) */
/*                     { */
/*                       mpc_get_cdpe (tmpc, s->mroot[l]); */
/*                       cdpe_mod (tmpr, tmpc); */
/*                       if (rdpe_lt (tmpr, rdpe_one)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'o':      /\* outside unit circle *\/ */
/*                   if (!mps_mtouchunit (s, nf, l)) */
/*                     { */
/*                       mpc_get_cdpe (tmpc, s->mroot[l]); */
/*                       cdpe_mod (tmpr, tmpc); */
/*                       if (rdpe_gt (tmpr, rdpe_one)) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'l':      /\* left half plane  *\/ */
/*                   if (!mps_mtouchimag (s, nf, l)) */
/*                     { */
/*                       /\* mpc_get_re(tmpf, s->mroot[l]); *\/ */
/*                       mpf_set (tmpf, mpc_Re (s->mroot[l])); */
/*                       if (mpf_sgn (tmpf) == -1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'r':      /\* right half plane *\/ */
/*                   if (!mps_mtouchimag (s, nf, l)) */
/*                     { */
/*                       /\* mpc_get_re(tmpf, s->mroot[l]); *\/ */
/*                       mpf_set (tmpf, mpc_Re (s->mroot[l])); */
/*                       if (mpf_sgn (tmpf) == 1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'u':      /\* upper half plane *\/ */
/*                   if (!mps_mtouchreal (s, nf, l)) */
/*                     { */
/*                       /\* mpc_get_im(tmpf, s->mroot[l]); *\/ */
/*                       mpf_set (tmpf, mpc_Im (s->mroot[l])); */
/*                       if (mpf_sgn (tmpf) == 1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'd':      /\* lower half plane *\/ */
/*                   if (!mps_mtouchreal (s, nf, l)) */
/*                     { */
/*                       /\* mpc_get_im(tmpf, s->mroot[l]); *\/ */
/*                       mpf_set (tmpf, mpc_Im (s->mroot[l])); */
/*                       if (mpf_sgn (tmpf) == -1) */
/*                         s->status[l][2] = 'i'; */
/*                       else */
/*                         s->status[l][2] = 'o'; */
/*                     } */
/*                   break; */

/*                 case 'R':      /\* Real line  NEW *\/ */
/*                   if (s->status[l][1] != 'w') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_mtouchreal (s, 1, l)) */
/*                         { */
/*                           if (s->data_type[1] == 'r') */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'R'; */
/*                             } */
/*                           else */
/*                             { */
/*                               rdpe_set (sr, s->drad[l]); */
/*                               /\* msrad(i, sc, sr);**#DARIO *\/ */
/*                               if (rdpe_log (sr) < s->sep - s->n */
/*                                   * s->lmax_coeff) */
/*                                 { */
/*                                   s->status[l][2] = 'i'; */
/*                                   s->status[l][1] = 'R'; */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   s->status[l][2] = 'u'; */
/*                                   s->status[l][1] = 'w'; */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch real *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'r'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_mtouchreal (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_mtouchreal (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\*  mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano R *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_msrad (s, i, sc, sr); */
/*                               sep1 = s->sep; */
/*                               if (s->data_type[1] == 'c') */
/*                                 sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (rdpe_log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'R'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano R *\/ */

/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][2] = 'o'; */
/*                               s->status[l2][1] = 'r'; */
/*                             } */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'I': */
/*                   if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                     continue; */
/*                   if (s->punt[i + 1] - s->punt[i] == 1) */
/*                     {           /\* one disk *\/ */
/*                       if (mps_mtouchimag (s, nf, l)) */
/*                         { */
/*                           /\* msrad(i, sc, sr); *\//\*#DARIO *\/ */
/*                           rdpe_set (sr, s->drad[l]); */
/*                           if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             { */
/*                               s->status[l][2] = 'i'; */
/*                               s->status[l][1] = 'I'; */
/*                             } */
/*                           else */
/*                             { */
/*                               s->status[l][2] = 'u'; */
/*                               s->status[l][1] = 'w'; */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* do not touch imag *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                       continue; */
/*                     } */
/*                   else */
/*                     {           /\* cluster *\/ */
/*                       tcr = mps_mtouchimag (s, nf, l); */
/*                       for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                         { */
/*                           l1 = s->clust[s->punt[i] + j1]; */
/*                           tcr1 = mps_mtouchimag (s, nf, l1); */
/*                           if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                             continue; */
/*                           else */
/*                             {   /\* mixed situation *\/ */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                               goto scan; */
/*                             } */
/*                         } */
/*                       if (tcr) */
/*                         {       /\* tutti i dischi intersecano I *\/ */
/*                           if (s->data_type[2] == 'f' || s->data_type[2] */
/*                               == 'b') */
/*                             { */
/*                               for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                    j2++) */
/*                                 { */
/*                                   l2 = s->clust[s->punt[i] + j2]; */
/*                                   s->status[l2][2] = 'u'; */
/*                                   s->status[l2][1] = 'w'; */
/*                                 } */
/*                             } */
/*                           else */
/*                             {   /\* integer/rational polynomial *\/ */
/*                               mps_msrad (s, i, sc, sr); */
/*                               sep1 = s->sep; */
/*                               sep1 = s->sep - s->n * s->lmax_coeff; */
/*                               if (rdpe_log (sr) < sep1) */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'i'; */
/*                                       s->status[l2][1] = 'I'; */
/*                                     } */
/*                                 } */
/*                               else */
/*                                 { */
/*                                   for (j2 = 0; j2 < s->punt[i + 1] */
/*                                        - s->punt[i]; j2++) */
/*                                     { */
/*                                       l2 = s->clust[s->punt[i] + j2]; */
/*                                       s->status[l2][2] = 'u'; */
/*                                       s->status[l2][1] = 'w'; */
/*                                     } */
/*                                 } */
/*                             } */
/*                         } */
/*                       else */
/*                         {       /\* tutti i dischi non intersecano I *\/ */
/*                           s->status[l][2] = 'o'; */
/*                           s->status[l][1] = 'i'; */
/*                         } */
/*                     } */
/*                   break; */

/*                 case 'S':      /\* Set provided by the user *\/ */
/*                   mps_error (s, 1, "Custom region not implemented yet"); */
/*                   break; */
/*                 default: */
/*                   mps_error (s, 1, "mistake in goal"); */
/*                   break; */
/*                 } */
/*             } */
/*         } */
/*       /\* If some cluster still contains an uncertain disk then set */
/*        * all the disks uncertain */
/*        *\/ */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         { */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][2] == 'u') */
/*             { */
/*               for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                 { */
/*                   l1 = s->clust[s->punt[i] + j1]; */
/*                   s->status[l1][2] = 'u'; */
/*                 } */
/*               break; */
/*             } */
/*         } */

/*     scan:; */
/*     } */

/*         /\*==3==  now check the options *\/ */

/*   /\* Option multiplicity *\/ */
/*   for (i = 0; i < s->nclust; i++) */
/*     {                           /\* scan1 *\/ */
/*       if (s->punt[i + 1] - s->punt[i] == 1) */
/*         continue; */
/*       for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*         {                       /\* scan1_in *\/ */
/*           l = s->clust[s->punt[i] + j]; */
/*           if (s->status[l][0] == 'x' || s->status[l][0] == 'f' */
/*               || s->status[l][0] == 'm') */
/*             goto scan1; */
/*           if (s->goal[2] == 'm' && (s->status[l][0] == 'c' ||   /\* NEW *\/ */
/*                                     s->status[l][0] == 'C')) */
/*             {                   /\* multiplicity on *\/ */

/*               if (s->data_type[2] == 'b' || s->data_type[2] == 'f')     /\* float coeff. *\/ */
/*                 mps_error (s, 1, */
/*                            "Fatal: Float coefficients - impossible to detect multiplicity"); */

/*               /\* compute super center and super radius *\/ */
/*               mps_msrad (s, i, sc, sr); */

/*               if (rdpe_log (sr) < s->sep) */
/*                 for (j1 = 0; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                   { */
/*                     l1 = s->clust[s->punt[i] + j1];     /\* NEW j-> j1 *\/ */
/*                     s->status[l1][0] = 'm'; */
/*                   } */
/*               goto scan1; */
/*             } */
/*         } */
/*     scan1:;                    /\* scan next component *\/ */
/*     } */

/*   /\* Option Real check *\/ */
/*   if ((s->goal[3] == 'r' || s->goal[3] == 'b') && s->goal[1] != 'R') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan2 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan2_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'i') */
/*                 continue; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan2; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_mtouchreal (s, nf, l)) */
/*                     { */
/*                       if (s->data_type[1] == 'r') */
/*                         s->status[l][1] = 'R'; */
/*                       else */
/*                         { */
/*                           /\* msrad(i, sc, sr); *\//\*#DARIO *\/ */
/*                           rdpe_set (sr, s->drad[l]); */
/*                           if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                             s->status[l][1] = 'R'; */
/*                           else */
/*                             s->status[l][1] = 'w'; */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* do not touch real *\/ */
/*                     s->status[l][1] = 'r'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_mtouchreal (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_mtouchreal (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan2_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano R *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_msrad (s, i, sc, sr); */
/*                           sep1 = s->sep; */
/*                           if (s->data_type[1] == 'c') */
/*                             sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (rdpe_log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'R'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano R *\/ */
/*                     s->status[l][1] = 'r'; */
/*                 } */
/*             scan2_in:; */
/*             } */
/*         scan2:; */
/*         } */
/*     } */
/*   /\* Option Imaginary check *\/ */
/*   if (s->goal[3] == 'i' || s->goal[3] == 'b') */
/*     { */
/*       for (i = 0; i < s->nclust; i++) */
/*         {                       /\* scan3 *\/ */
/*           for (j = 0; j < s->punt[i + 1] - s->punt[i]; j++) */
/*             {                   /\* scan3_in *\/ */
/*               l = s->clust[s->punt[i] + j]; */
/*               if (s->status[l][0] == 'x' || s->status[l][0] == 'f') */
/*                 goto scan3; */
/*               if (s->status[l][1] != 'w' && s->status[l][1] != 'r') */
/*                 continue; */
/*               if (s->punt[i + 1] - s->punt[i] == 1) */
/*                 {               /\* one disk *\/ */
/*                   if (mps_mtouchimag (s, nf, l)) */
/*                     { */
/*                       rdpe_set (sr, s->drad[l]); */
/*                       /\* msrad(i, sc, sr); DARIO *\/ */
/*                       if (rdpe_log (sr) < s->sep - s->n * s->lmax_coeff) */
/*                         s->status[l][1] = 'I'; */
/*                       else */
/*                         s->status[l][1] = 'w'; */
/*                     } */
/*                   else */
/*                     /\* do not touch imag *\/ */
/*                     s->status[l][1] = 'i'; */
/*                   continue; */
/*                 } */
/*               else */
/*                 {               /\* cluster *\/ */
/*                   tcr = mps_mtouchimag (s, nf, l); */
/*                   for (j1 = 1; j1 < s->punt[i + 1] - s->punt[i]; j1++) */
/*                     { */
/*                       l1 = s->clust[s->punt[i] + j1]; */
/*                       tcr1 = mps_mtouchimag (s, nf, l1); */
/*                       if ((tcr && tcr1) || (!tcr && !tcr1)) */
/*                         continue; */
/*                       else */
/*                         {       /\* mixed situation *\/ */
/*                           for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                             { */
/*                               l2 = s->clust[s->punt[i] + j2]; */
/*                               s->status[l2][1] = 'w'; */
/*                             } */
/*                           goto scan3_in; */
/*                         } */
/*                     } */
/*                   if (tcr) */
/*                     {           /\* tutti i dischi intersecano I *\/ */
/*                       if (s->data_type[2] == 'f' || s->data_type[2] == 'b') */
/*                         for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; j2++) */
/*                           { */
/*                             l2 = s->clust[s->punt[i] + j2]; */
/*                             s->status[l2][1] = 'w'; */
/*                           } */
/*                       else */
/*                         {       /\* integer/rational polynomial *\/ */
/*                           mps_msrad (s, i, sc, sr); */
/*                           sep1 = s->sep; */
/*                           sep1 = s->sep - s->n * s->lmax_coeff; */
/*                           if (rdpe_log (sr) < sep1) */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'I'; */
/*                               } */
/*                           else */
/*                             for (j2 = 0; j2 < s->punt[i + 1] - s->punt[i]; */
/*                                  j2++) */
/*                               { */
/*                                 l2 = s->clust[s->punt[i] + j2]; */
/*                                 s->status[l2][1] = 'w'; */
/*                               } */
/*                         } */
/*                     } */
/*                   else */
/*                     /\* tutti i dischi non intersecano I *\/ */
/*                     s->status[l][1] = 'i'; */
/*                 } */
/*             scan3_in:; */
/*             } */
/*         scan3:; */
/*         } */
/*     } */

/*   mpf_clear (tmpf); */
/*   mpc_clear (sc); */
/* } */
