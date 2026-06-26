/* 
 ---------------------------------------------------------------------------
 Copyright (c) 1998-2010, Brian Gladman, Worcester, UK. All rights reserved.

 The redistribution and use of this software (with or without changes)
 is allowed without the payment of fees or royalties provided that:

   source code distributions include the above copyright notice, this
   list of conditions and the following disclaimer;

   binary distributions include the above copyright notice, this list
   of conditions and the following disclaimer in their documentation.

 This software is provided 'as is' with no explicit or implied warranties
 in respect of its operation, including, but not limited to, correctness
 and fitness for purpose.
 ---------------------------------------------------------------------------
*/

#include <malloc.h>
#include "bessel.h"
#include "mathieu.h"

#define NORMALISED 0
/*
 This routine, which is a member function of mathieu, computes radial solutions 
 of Mathieu's equation of the 1st and 2nd kinds which correspond with periodic
 angular solutions.  These solutions are mcm(1)(x,q), mcm(2)(x,q), msm(1)(x,q) 
 and msm(2)(x,q) together with their derivatives. The parameter q and the root
 values, r and q, are bound into mathieu itself so that the function only needs
 the x value as a parameter together with a specifier of the outputs which are
 to be returned.

 Input:	    kc	    function code
            kc = 1 	for computing Mcm(1) (x,q) or Msm(1) (x,q)
            kc = 2	for computing Mcm(2) (x,q) or Msm(2) (x,q)
            kc = 4	for computing Mcm(1)'(x,q) or Msm(1)'(x,q)
            kc = 8	for computing Mcm(2)'(x,q) or Msm(2)'(x,q)

            x	    function argument

  Output:	 rf1	Mcm(1)(x,q) or Msm(1)(x,q)
             rd1	Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
             rf2	Mcm(2)(x,q) or Msm(2)(x,q)
             rd2	Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
*/

void mathieu_radial(int kc, const double x, double *rf1, double *rf2, double *rd1, double *rd2, mathieu_qr w[1])
{
    int		k, x_s, i_s, i_d, n_bf, xf;
    double	u1, u2, f, ff, sn_s, sn_k, sgn_2, s_k, *j1, *j2, *y2;

    if(!w->coeff)
        mathieu_expansion_coeffs(w);

    x_s = w->c_max + 1; 
    sgn_2 = (w->type == type_b ? -1.0 : 1.0); 

#if NORMALISED
    ff = sqrt(2.0 / M_PI) * w->coeff[w->c_max]; 
#else
    ff = w->coeff[w->c_max]; 
#endif

    sn_s = ((w->c_max + w->r / 2) & 1 ? -1.0 : 1.0); 
    sn_k = (w->r & 2 ? -1.0 : 1.0);

    if(!(w->r & 1))
    {
        if(w->type == type_b)
        {		
            x_s++; 
            sn_s = -sn_s; 
            sn_k = -sn_k;
        }
        else
        {
            x_s--; 
            ff = (w->c_max ? ff : ff + ff);
        }
    }

    sn_s /= ff; 
    sn_k /= ff; 
    n_bf = x_s + w->c_len - 1;

    u1 = sqrt(w->q) * exp(-x); 
    u2 = sqrt(w->q) * exp(x);
    
    if(!(w->bf1.b_nof) || u1 != w->bf1.b_arg || n_bf > w->bf1.b_nof)
    {
        if(n_bf > w->bf1.b_nof)
        {
            if(w->bf1.b_jnp)
                free(w->bf1.b_jnp); 
            if(w->bf1.b_ynp)
                free(w->bf1.b_ynp); 
            w->bf1.b_nof = w->c_len;
            w->bf1.b_jnp = malloc((n_bf + 1) * sizeof(double)); 
            w->bf1.b_ynp = malloc((n_bf + 1) * sizeof(double)); 
        }

        bess(u1, n_bf, w->bf1.b_jnp, w->bf1.b_ynp); 
    }

    if(!(w->bf2.b_nof) || u2 != w->bf2.b_arg || n_bf > w->bf2.b_nof)
    {
        if(n_bf > w->bf2.b_nof)
        {
            if(w->bf2.b_jnp)
                free(w->bf2.b_jnp); 
            if(w->bf2.b_ynp)			
                free(w->bf2.b_ynp); 
            w->bf2.b_nof = w->c_len;
            w->bf2.b_jnp = malloc((n_bf + 1) * sizeof(double)); 
            w->bf2.b_ynp = malloc((n_bf + 1) * sizeof(double)); 
        }

        bess(u2, n_bf, w->bf2.b_jnp, w->bf2.b_ynp);
    }
    
    j1 = w->bf1.b_jnp; 
    j2 = w->bf2.b_jnp; 
    y2 = w->bf2.b_ynp;

    if(kc & 1)
    {
        k = 0; 
        f = 0.0; 
        s_k = sn_k; 
        i_s = x_s;

        do
        {
            ff = f; 
            i_d = (k < w->c_max ? w->c_max - k : k - w->c_max);

            f += (k < w->c_max ? sn_s : s_k) * w->coeff[k] 
                        * (j1[i_d] * j2[i_s] + sgn_2 * j1[i_s] * j2[i_d]);
            s_k = -s_k; 
            i_s++;
            xf = (k >= x_s && f == ff ? xf + 1 : 0); 
        }
        while
            (++k < w->c_len && xf < 4);

        *rf1 = f; 
    }

    if(kc & 2)
    {
        k = 0; 
        f = 0.0; 
        s_k = sn_k; 
        i_s = x_s;

        do
        {
            ff = f; 
            i_d = (k < w->c_max ? w->c_max - k : k - w->c_max);

            f += (k < w->c_max ? sn_s : s_k) * w->coeff[k] 
                        * (j1[i_d] * y2[i_s] + sgn_2 * j1[i_s] * y2[i_d]);
            s_k = -s_k;
            i_s++;
            xf = (k >= x_s && f == ff ? xf + 1 : 0); 
        }
        while
            (++k < w->c_len && xf < 4);

        *rf2 = f; 
    }

    if(kc & 4)
    {
        k = 0; 
        f = 0.0; 
        s_k = sn_k; 
        i_s = x_s;

        do
        {
            ff = f; 
            i_d = (k < w->c_max ? w->c_max - k : k - w->c_max);

            f += (k < w->c_max ? sn_s : s_k) * w->coeff[k] * 
                 ((i_s - i_d) * (j1[i_d] * j2[i_s] - sgn_2 * j1[i_s] * j2[i_d]) 
                + u1 * (j1[i_d + 1] * j2[i_s] + sgn_2 * j1[i_s + 1] * j2[i_d])
                - u2 * (j1[i_d] * j2[i_s + 1] + sgn_2 * j1[i_s] * j2[i_d + 1]));
            s_k = -s_k; 
            i_s++;
            xf = (k >= x_s && f == ff ? xf + 1 : 0); 
        }
        while
            (++k < w->c_len && xf < 4);

        *rd1 = f; 
    }

    if(kc & 8)
    {
        k = 0; 
        f = 0.0; 
        s_k = sn_k; 
        i_s = x_s;

        do
        {			
            ff = f; 
            i_d = (k < w->c_max ? w->c_max - k : k - w->c_max);

            f += (k < w->c_max ? sn_s : s_k) * w->coeff[k] * 
                 ((i_s - i_d) * (j1[i_d] * y2[i_s] - sgn_2 * j1[i_s] * y2[i_d]) 
                + u1 * (j1[i_d + 1] * y2[i_s] + sgn_2 * j1[i_s + 1] * y2[i_d])
                - u2 * (j1[i_d] * y2[i_s + 1] + sgn_2 * j1[i_s] * y2[i_d + 1]));
            s_k = -s_k; 
            i_s++;
            xf = (k >= x_s && f == ff ? xf + 1 : 0); 
        }
        while
            (++k < w->c_len && xf < 4);

        *rd2 = f; 
    }

    return;
}
