/***********************************************************
**       Multiprecision Polynomial Solver (MPSolve)       **
**                Version 2.01, April 1998                **
**                                                        **
**                      Written by                        **
**       Dario Andrea Bini and Giuseppe Fiorentino        **
**       (bini@dm.unipi.it)  (fiorent@di.unipi.it)        **
**                                                        **
** (C) 1998, Dipartimento di Matematica, FRISCO LTR 21024 **
***********************************************************/

/********************************************************
This file contains a sample of a user-defined polynomial.
In this sample the Deflated Chromatic  Polynomials of kind a
are computed. Their degrees are 0, 4, 20, 84, 340, 1356, ... 
(n_0=0, n_{i+1}=4 n_i+4).
Deflated chromatic polynomials of kind a,c and d, are
defined by means of the relation

na=(x-1)^3 a^4 + c^4, nc= c^2 d, d=[(x-2) c^2 +2(x-1)a^2]
with the initial values
a=1, c=(x-2)
The number of steps is i=log(3 dega+4)/log(4)-1, i=log((3degd-1)/2)/log 4
Derivatives:
NA'= 3(x-1)^2 a^4+(x-1)^34 a^3 a' +4c^3 c'
ND'=c^2+2[(x-2)cc'+2(x-1)a^2+2(x-1)^2aa']
NC'=2cc'd+c^2d'
Error bounds:
EA=S1|x-1|^3a^4+S2 c^4 +4(x-1)^3a^3 EA +4 c^3 EC
ED=S4 |x-2| c^2 + 2S5 (x-1)^2 a^2 + 2c(x-2) EC +4(x-1)^2 a EA
EC= c^2 d S3 + 2cd EC +c^2 ED 

S1=(10 sqrt 2 +2)u<16.15u, S2=(4sqrt 2+1)u<6.66u, S3=(4 sqrt 2)u<5.66u, S4=(4 sqrt 2+2)u<7.66u
S5=(6 sqrt 2+2)u<10.49u
The three programs below compute the values of a(x), a'(x) and
c(x), c'(x), d(x), d'(x) by means of the above formulae in the float 
version, dpe version and mp version. 
THESE PROGRAMS OUTPUT THE VALUE corr=a/a'

The degrees of the polynomials a,c and d are:

a 0 4 20 84 340 1356 ... 
c 1 5 21 85 341 1355 ...
d   3 11 43 171  683 ...
**********************************************************/

#include "mps.h"

/******************************************************
*              SUBROUTINE FNEWTON_USR                 *   
******************************************************* 
******************************************************/
void
fnewton_usr(cplx_t x, double *rad, cplx_t corr, mps_boolean * again)
{

{
  cplx_t xm1,xm2,a,as,ac,aq,c,cs,cc,cq,ap,cp,d,dp,ccp,aap;
  cplx_t xm1s,xm1s2,xm1c,na,nc,ctmp1,ctmp2;
  cplx_t xm1s3;  
  double eps,ea,ec,ed,ax,axm1,axm1s,axm1c,axm2;
  double absa,absas,absac,absaq,absc,abscs,abscc,abscq;
  double rtmp1,rtmp2, ean,ecn,absd;
  double s1,s2,s3,s4,s5;
  int i, m;
 
  /* deg a */
  m=(int) (log(3*n+4.0)/log(4.0)-1+0.4);

  /* deg d */
/*
  m=(int) (log((3*n-1.0)/2.0)/log(4.0)+0.4); 
*/

/* initial conditions */
/* polynomials */
  cplx_sub(xm1,x,cplx_one);
  cplx_sub(xm2,xm1,cplx_one);
  cplx_set(c,xm2);
  cplx_set(a,cplx_one);

/* derivatives */
  cplx_set(ap,cplx_zero);
  cplx_set(cp,cplx_one);

/* error bounds */
  eps= DBL_EPSILON;ea=0;ec=0;
  s1 = 16.15 * DBL_EPSILON;
  s2 = 6.66 * DBL_EPSILON;
  s3 = 5.66 * DBL_EPSILON;
  s4 = 7.66 * DBL_EPSILON;
  s5 = 10.49 * DBL_EPSILON;
  ax= cplx_mod(x);

/* auxiliary terms */
  cplx_mul(xm1s,xm1,xm1);        // (x-1)^2
  cplx_add(xm1s2,xm1s,xm1s);   // 2(x-1)^2
  cplx_add(xm1s3,xm1s2,xm1s);     // 3(x-1)^2
  cplx_mul(xm1c,xm1,xm1s);       // (x-1)^3
  axm1=cplx_mod(xm1);            // |x-1|
  axm1s=axm1*axm1;     // |x-1|^2
  axm1c=axm1*axm1s;    // |x-1|^3
  axm2=cplx_mod(xm2);            // |x-2|

/* computation */
  for (i = 1; i <= m; i++) {
   /* polynomial */
   cplx_mul(as,a,a); cplx_mul(aq,as,as);
   cplx_mul(cs,c,c); cplx_mul(cq,cs,cs);
   cplx_mul(cc,cs,c); cplx_mul(ac,as,a);
   cplx_mul(na,xm1c,aq); cplx_add(na,na,cq);
   cplx_mul(d,cs,xm2);cplx_mul(ctmp1,xm1s2,as);
   cplx_add(d,d,ctmp1);
   cplx_mul(nc,d,cs); 
   
   /* derivative  */
   /* a */
   cplx_mul(aap,a,ap);
   cplx_mul(ccp,c,cp);
   cplx_mul(ap,xm1c,as);
   cplx_mul(ap,ap,aap);
   cplx_mul(ctmp1,cs,ccp);
   cplx_add(ap,ap,ctmp1);cplx_mul_d(ap,ap,4);
   cplx_mul(ctmp1,xm1s3,aq);
   cplx_add(ap,ap,ctmp1); // a'

   /* d */
   cplx_mul(ctmp1,ccp,xm2);cplx_mul(ctmp2,xm1,as);
   cplx_add(ctmp2,ctmp2,ctmp2);
   cplx_add(ctmp1,ctmp1,ctmp2);cplx_mul(ctmp2,xm1s2,aap);
   cplx_add(ctmp1,ctmp1,ctmp2);cplx_add(ctmp1,ctmp1,ctmp1);
   cplx_add(dp,cs,ctmp1); // d'

   /* c */
   cplx_mul(cp,cs,dp); cplx_mul(ctmp1,ccp,d);cplx_add(ctmp1,ctmp1,ctmp1);
   cplx_add(cp,cp,ctmp1);


   /* error bound */
   absa=cplx_mod(a);
   absc=cplx_mod(c);
   absd=cplx_mod(d);
   absas=absa*absa; absac=absas*absa;
   absaq=absas*absas;
   abscs=absc*absc;abscc=abscs*absc;
   abscq=abscs*abscs;

   rtmp1=absac*ea;rtmp1=4*rtmp1;
   ean=absaq*s1;ean=ean+rtmp1;
   ean=ean*axm1c;
   rtmp1=abscq*s2;ean=ean+rtmp1;
   rtmp1=abscc*ec;
   rtmp1=4*rtmp1;ean=rtmp1+ean; //ean

   rtmp1=abscs*s4;ed=absc*ec;
   ed=ed*2;
   ed=ed+rtmp1; ed=ed*axm2;
   rtmp1=absas*s5;rtmp2=absa*ea;
   rtmp2=2*rtmp2;rtmp1=rtmp1+rtmp2;
   rtmp1=rtmp1*axm1s;rtmp1=2*rtmp1;
   ed=ed+rtmp1; //ed

   rtmp1=absd*s3;rtmp1=rtmp1+ed;
   rtmp1=rtmp1*abscs;
   ecn=ec*absd;ecn=ecn*absc;ecn=2*ecn;
   ecn=ecn+rtmp1; //ecn


   /* shift */
   cplx_set(a,na);
   cplx_set(c,nc);
   ea=ean;ec=ecn;
  }


/* Polynomial  a(x) */

  absa=cplx_mod(a); 
  cplx_div(corr,a,ap);
  *again = (ean<absa); //rdpe_lt(ean, absa);
  rtmp1=cplx_mod(ap);
  absa=n*(absa+ean)/rtmp1;
  *rad=absa;
}


}

/******************************************************
*              SUBROUTINE DNEWTON_USR                 *   
******************************************************* 
 DPE computation
******************************************************/
void
dnewton_usr(cdpe_t x, rdpe_t rad, cdpe_t corr, mps_boolean * again)
{
  cdpe_t xm1,xm2,a,as,ac,aq,c,cs,cc,cq,ap,cp,d,dp,ccp,aap;
  cdpe_t xm1s,xm1s2,xm1c,na,nc,ctmp1,ctmp2;
  cdpe_t xm1s3;  
  rdpe_t eps,ea,ec,ed,ax,axm1,axm1s,axm1c,axm2;
  rdpe_t absa,absas,absac,absaq,absc,abscs,abscc,abscq;
  rdpe_t rtmp1,rtmp2, ean,ecn,absd;
  double s1,s2,s3,s4,s5;
  int i, m;
 
  /* deg a */
  m=(int) (log(3*n+4.0)/log(4.0)-1+0.4);

  /* deg d */
/*
  m=(int) (log((3*n-1.0)/2.0)/log(4.0)+0.4); 
*/


/* initial conditions */
/* polynomials */
  cdpe_sub(xm1,x,cdpe_one);
  cdpe_sub(xm2,xm1,cdpe_one);
  cdpe_set(c,xm2);
  cdpe_set(a,cdpe_one);

/* derivatives */
  cdpe_set(ap,cdpe_zero);
  cdpe_set(cp,cdpe_one);

/* error bounds */
  rdpe_set_d(eps, DBL_EPSILON);
  rdpe_set(ea,rdpe_zero);
  rdpe_set(ec,rdpe_zero);
  s1 = 16.15 * DBL_EPSILON;
  s2 = 6.66 * DBL_EPSILON;
  s3 = 5.66 * DBL_EPSILON;
  s4 = 7.66 * DBL_EPSILON;
  s5 = 10.49 * DBL_EPSILON;
  cdpe_mod(ax, x);

/* auxiliary terms */
  cdpe_mul(xm1s,xm1,xm1);        // (x-1)^2
  cdpe_mul_2exp(xm1s2,xm1s,1);   // 2(x-1)^2
  cdpe_add(xm1s3,xm1s2,xm1s);     // 3(x-1)^2
  cdpe_mul(xm1c,xm1,xm1s);       // (x-1)^3
  cdpe_mod(axm1,xm1);            // |x-1|
  rdpe_mul(axm1s,axm1,axm1);     // |x-1|^2
  rdpe_mul(axm1c,axm1,axm1s);    // |x-1|^3
  cdpe_mod(axm2,xm2);            // |x-2|

/* computation */
  for (i = 1; i <= m; i++) {
   /* polynomial */
   cdpe_mul(as,a,a); cdpe_mul(aq,as,as);
   cdpe_mul(cs,c,c); cdpe_mul(cq,cs,cs);
   cdpe_mul(cc,cs,c); cdpe_mul(ac,as,a);
   cdpe_mul(na,xm1c,aq); cdpe_add(na,na,cq);
   cdpe_mul(d,cs,xm2);cdpe_mul(ctmp1,xm1s2,as);
   cdpe_add(d,d,ctmp1);
   cdpe_mul(nc,d,cs); 
   
   /* derivative  */
   /* a */
   cdpe_mul(aap,a,ap);
   cdpe_mul(ccp,c,cp);
   cdpe_mul(ap,xm1c,as);
   cdpe_mul(ap,ap,aap);
   cdpe_mul(ctmp1,cs,ccp);
   cdpe_add(ap,ap,ctmp1);cdpe_mul_2exp(ap,ap,2);
   cdpe_mul(ctmp1,xm1s3,aq);
   cdpe_add(ap,ap,ctmp1); // a'

   /* d */
   cdpe_mul(ctmp1,ccp,xm2);cdpe_mul(ctmp2,xm1,as);
   cdpe_mul_2exp(ctmp2,ctmp2,1);
   cdpe_add(ctmp1,ctmp1,ctmp2);cdpe_mul(ctmp2,xm1s2,aap);
   cdpe_add(ctmp1,ctmp1,ctmp2);cdpe_mul_2exp(ctmp1,ctmp1,1);
   cdpe_add(dp,cs,ctmp1); // d'

   /* c */
   cdpe_mul(cp,cs,dp); cdpe_mul(ctmp1,ccp,d);cdpe_mul_2exp(ctmp1,ctmp1,1);
   cdpe_add(cp,cp,ctmp1);


   /* error bound */
   cdpe_mod(absa,a);
   cdpe_mod(absc,c);
   cdpe_mod(absd,d);
   rdpe_mul(absas,absa,absa); rdpe_mul(absac,absas,absa);
   rdpe_mul(absaq,absas,absas);
   rdpe_mul(abscs,absc,absc); rdpe_mul(abscc,abscs,absc);
   rdpe_mul(abscq,abscs,abscs);

   rdpe_mul(rtmp1,absac,ea);rdpe_mul_2exp(rtmp1,rtmp1,2);
   rdpe_mul_d(ean,absaq,s1);rdpe_add(ean,ean,rtmp1);
   rdpe_mul(ean,ean,axm1c);
   rdpe_mul_d(rtmp1,abscq,s2);rdpe_add(ean,ean,rtmp1);
   rdpe_mul(rtmp1,abscc,ec);
   rdpe_mul_2exp(rtmp1,rtmp1,2);rdpe_add(ean,rtmp1,ean); //ean

   rdpe_mul_d(rtmp1,abscs,s4);rdpe_mul(ed,absc,ec);
   rdpe_mul_2exp(ed,ed,1);
   rdpe_add(ed,ed,rtmp1); rdpe_mul(ed,ed,axm2);
   rdpe_mul_d(rtmp1,absas,s5);rdpe_mul(rtmp2,absa,ea);
   rdpe_mul_2exp(rtmp2,rtmp2,1);rdpe_add(rtmp1,rtmp1,rtmp2);
   rdpe_mul(rtmp1,rtmp1,axm1s);rdpe_mul_2exp(rtmp1,rtmp1,1);
   rdpe_add(ed,ed,rtmp1); //ed

   rdpe_mul_d(rtmp1,absd,s3);rdpe_add(rtmp1,rtmp1,ed);
   rdpe_mul(rtmp1,rtmp1,abscs);
   rdpe_mul(ecn,ec,absd);rdpe_mul(ecn,ecn,absc);rdpe_mul_2exp(ecn,ecn,1);
   rdpe_add(ecn,ecn,rtmp1); //ecn


   /* shift */
   cdpe_set(a,na);
   cdpe_set(c,nc);
   rdpe_set(ea,ean);rdpe_set(ec,ecn);
  }

/*
fprintf(outstr,"x=");cdpe_out_str(outstr,x);fprintf(outstr,"\n");
fprintf(outstr,"a=");cdpe_out_str(outstr,a);fprintf(outstr,"\n");
fprintf(outstr,"c=");cdpe_out_str(outstr,c);fprintf(outstr,"\n");
fprintf(outstr,"d=");cdpe_out_str(outstr,d);fprintf(outstr,"\n");
fprintf(outstr,"ap=");cdpe_out_str(outstr,ap);fprintf(outstr,"\n");
fprintf(outstr,"cp=");cdpe_out_str(outstr,cp);fprintf(outstr,"\n");
fprintf(outstr,"dp=");cdpe_out_str(outstr,dp);fprintf(outstr,"\n");


fprintf(outstr,"ea=");rdpe_out_str(outstr,ean);fprintf(outstr,"\n");
fprintf(outstr,"absa=");rdpe_out_str(outstr,absa);fprintf(outstr,"\n");
fprintf(outstr,"\n");
*/

/* Polynomial  a(x) */

  cdpe_mod(absa,a); 
  cdpe_div(corr,a,ap);
  *again = rdpe_lt(ean, absa);
  cdpe_mod(rtmp1,ap);
  rdpe_add(rad,absa,ea);
  rdpe_div_eq(rad,  rtmp1);
  rdpe_mul_eq_d(rad, (double) n);
}


/******************************************************
*              SUBROUTINE MNEWTON_USR                 *   
******************************************************* 
 multiprecision computation
******************************************************/
void
mnewton_usr(mpc_t x, rdpe_t rad, mpc_t corr, mps_boolean * again)
{
  rdpe_t d0, d1, d2, ap0, ap1,ap2, ax, eps, rtmp, rtmp1;
  int i, m;
  cdpe_t ctmp;
  tmpc_t p0,p1,p2,pp0,pp1,pp2,tmp1,tmp2;

  tmpc_init2(p0, mpwp);
  tmpc_init2(p1, mpwp);
  tmpc_init2(p2, mpwp);
  tmpc_init2(pp0, mpwp);
  tmpc_init2(pp1, mpwp);
  tmpc_init2(pp2, mpwp);
  tmpc_init2(tmp1, mpwp);
  tmpc_init2(tmp2, mpwp);



  rdpe_set(eps, mp_epsilon);
  m = (int) (log(n ) / LOG2);
  m=m+1;
  mpc_get_cdpe(ctmp,x);
  cdpe_mod(ax, ctmp);

/* initial conditions */
  mpc_set_d(p0, 1.0,0.0);    /* p0=1  */
  mpc_set_d(p1, 1.0,0.0);    /* p1=1  */
  mpc_set_d(pp0, 0.0,0.0);  /* p0'=0 */
  mpc_set_d(pp1,0.0,0.0);   /* p1'=0 */

  rdpe_set(ap0, rdpe_one);   /* |p0|=1 */
  rdpe_set(ap1, rdpe_one);  /* |p1|=1 */
  rdpe_set(d0,rdpe_zero);
  rdpe_set(d1,d0);

/* computation */



  for (i = 1; i <= m; i++) {
   /* polynomial */
    mpc_mul(tmp2, p0, p0);
    mpc_mul(pp2,tmp2,tmp2);
    mpc_mul(p2,pp2,x);
    mpc_mul(tmp1,p1,p1);
    mpc_add(p2,p2,tmp1);

   /* derivative  */
    mpc_mul(tmp1,p0,tmp2);
    mpc_get_cdpe(ctmp,tmp1);
    cdpe_mod(rtmp,ctmp);
    mpc_mul(tmp1,tmp1,pp0);
    mpc_mul(tmp1,x,tmp1);

    mpc_add(tmp1,tmp1,tmp1);
    mpc_add(tmp1,tmp1,tmp1);

    mpc_add(pp2,tmp1,pp2);
    mpc_mul(tmp1,p1,pp1);
    mpc_add(tmp1,tmp1,tmp1);
 
    mpc_add(pp2,tmp1,pp2);
    mpc_set(p0,p1);
    mpc_set(p1,p2);
    mpc_set(pp0,pp1);
    mpc_set(pp1,pp2);

   /* error bound */
   /* d2=2 d1 |p1| + 4 |x| d0 |p0|^3 + 
         u( |p2| + 2sqrt 2 |p1|^2 + |x| |p0|^4(4sqrt 2+1) */
    mpc_get_cdpe(ctmp,p2);
    cdpe_mod(ap2,ctmp);
    rdpe_mul(rtmp,rtmp,d0);
    rdpe_mul(rtmp,rtmp,ax);

    rdpe_mul_2exp(rtmp,rtmp,2);  /* 4 |x| |p0|^3 d0 */

    rdpe_mul(d2,d1,ap1);

    rdpe_mul_2exp(d2,d2,1);  /*  2 d1 |p1|  */
    rdpe_add(d2,d2,rtmp);

    rdpe_mul(rtmp1,ap0,ap0);
    rdpe_mul(rtmp1,rtmp1,rtmp1); /*  |p0|^4 */
    rdpe_mul(rtmp1,rtmp1,ax);
    rdpe_mul_d(rtmp1,rtmp1, (double) 6.66); /* |x|(4\sqrt 2 +1)|p0|^4 */
    rdpe_add(rtmp1,rtmp1,ap2);
    rdpe_mul(rtmp,ap1,ap1);
    rdpe_mul_d(rtmp,rtmp, (double) 2.83);
    rdpe_mul_2exp(rtmp,rtmp,1);

    rdpe_mul(rtmp,rtmp,eps);
    rdpe_add(d2,d2,rtmp);

    rdpe_set(d0,d1);
    rdpe_set(d1,d2);
    rdpe_set(ap0,ap1);
    rdpe_set(ap1,ap2);

  }

  mpc_div(corr,p1,pp1);
  *again = rdpe_lt(d2, ap2);
  mpc_get_cdpe(ctmp,pp1);
  cdpe_mod(rtmp,ctmp);
  rdpe_add(rad,ap2,d2);
  rdpe_div_eq(rad,  rtmp);
  rdpe_mul_eq_d(rad, (double) n);

  tmpc_clear(tmp2);
  tmpc_clear(tmp1);
  tmpc_clear(pp2);
  tmpc_clear(pp1);
  tmpc_clear(pp0);
  tmpc_clear(p2);
  tmpc_clear(p1);
  tmpc_clear(p0);
}
