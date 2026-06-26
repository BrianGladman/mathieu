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

 Compute expansion coefficients for Mathieu functions
 and modified Mathieu functions

 Input	mt_type    function type (see below)
              q	   parameter
              r	   function order
             cv    characteristic value (assumed correct)

 Output	coeff[k]   expansion coefficients of Mathieu
                   functions (r= 0, 1,...,kr - 1)
 
 where the function returns kr, the highest coefficient
 needed.  The w->coeff[w->r] correspond to:

  function	   r   period     coefficients ....
  even (1)  even       pi     A0    A2    A4   ....
   odd (2)   odd     2.pi	  B1    B3    B5   ....
  even (1)   odd     2.pi     A1    A3    A5   ....
   odd (2)  even	   pi     B2    B4    B6   ....
*/
 
#include "mathieu.h"
#include <malloc.h>

double mathieu_init_qr(const mt_type ty, const int r, const double q, mathieu_qr w[1]) 
{
    w->type = ty;
    w->r = r; 
    w->q = q;
    
    w->h.dp = w->h.ep = 0;
    w-> h.nm = 0;
    mathieu_blanch_values((int)w->type, w->q, w->r, w->r + 1, &w->cv); 
    w->coeff = 0;
    
    w->bf1.b_nof = 0; 
    w->bf1.b_jnp = w->bf1.b_ynp = 0;
    w->bf2.b_nof = 0; 
    w->bf2.b_jnp = w->bf2.b_ynp = 0;
    return w->cv;
}

void mathieu_free_qr(mathieu_qr w[1])
{
    if(w->coeff)
        free(w->coeff);

    w->coeff = 0;

    if(w->bf1.b_jnp)
        free(w->bf1.b_jnp);

    if(w->bf1.b_ynp)
        free(w->bf1.b_ynp);

    w->bf1.b_nof = 0; 
    w->bf1.b_jnp = 0;
    w->bf1.b_ynp = 0;

    if(w->bf2.b_jnp)
        free(w->bf2.b_jnp);

    if(w->bf2.b_ynp)
        free(w->bf2.b_ynp);

    free_hill_workspace(&w->h);

    w->bf2.b_nof = 0; 
    w->bf2.b_jnp = 0;
    w->bf2.b_ynp = 0;
}

void mathieu_expansion_coeffs(mathieu_qr w[1])
{	int		i, it, j, jo, k;
    double	u, v, f, f_hi, s_lo, s_hi;

    if(w->coeff)
        return;

    w->c_len = hill_determinant_size(w->q, w->r);
    w->coeff = malloc(w->c_len * sizeof(double));
    w->c_max = 0;

    if(fabs(w->q) <= 1.0e-15)	/* effectively zero */
    {
        for (i = 0; i < w->c_len; i++)
            w->coeff[i] = 0.0;

        u = 1.0;
        w->c_max = w->r / 2;
        if(w->r & 1)
            w->c_max = (w->r - 1) / 2;
        else if(w->type == type_b)
            w->c_max = w->r / 2 - 1;
        else
            u = (w->r ? 1.0 : sqrt(0.5));

        w->coeff[w->c_max] = u;
        return;
    }

    jo = (w->r & 1 ? 3 : (w->type == type_b ? 4 : 2));

    w->coeff[w->c_len - 1] = f = 1.0e-100; 
    u = 0.0; 
    s_hi = 0.0;

    /* perform downward recurrence while values increase */
    for(i = w->c_len - 2, j = 2 * i + jo; i >= 0; --i, j -= 2)
    {
        v = u; 
        u = f; 
        f = (w->cv - j * j) * u / w->q - v; 
        s_hi += u * u; w->coeff[i] = f;

        if(fabs(f) < fabs(u))
        {
            it = i + 1; 
            f_hi = w->coeff[i + 1]; 
            break;
        }
    }

    if(i < 0)	/* downward recurrence terminated at zero	*/
    {
        if((w->r & 1) || w->type == type_b)
            s_lo = s_hi + f * f;
        else
        {
            w->coeff[0] = 0.5 * f; 
            s_lo = s_hi + 0.5 * f * f;
        }

        f = 1.0; 
        it = w->c_len;
    }
    else	/* use forward recurrence for remaining terms	*/
    {	
        w->coeff[0] = u = 1.0e-100; v = 1.0;
        
        if(w->r & 1)		
            f = (w->cv - 1.0 - (((w->type == type_b ? ~w->r : w->r) & 1) ? 1.0 : -1.0) * w->q) / w->q * w->coeff[0];
        else if(w->type == type_b)
            f = (w->cv - 4.0) / w->q * w->coeff[0];
        else
        {
            f = w->cv / w->q * w->coeff[0]; 
            v = 2.0; 
            u = 2.0 * u;
        }

        s_lo = v * w->coeff[0] * w->coeff[0];
        
        for(i = 1, j = jo; i < it; ++i, j += 2)
        {
            if(s_lo > 1.0e+100)
            {
                s_lo *= 1.0e-200;
                v *= 1.0e-100;
                u *= 1.0e-100;
                f *= 1.0e-100;
                for(k = 0; k < i; ++k)
                    w->coeff[k] *= 1.0e-100;
            }
            w->coeff[i] = f; 
            s_lo += f * f; 
            v = u;  
            u = f; 
            f = (w->cv - j * j) * u / w->q - v;
        }

        /* adjust scale factors to join the two sequences	*/
        f /= f_hi; 
        s_lo += s_hi * f * f;
    }

    s_lo = (w->coeff[0] < 0.0 ? -1.0 : 1.0) / sqrt(s_lo); 
    s_hi = f * s_lo;

    /* scale the sequences	*/
    i = 0;
    while(i < it)
        w->coeff[i++] *= s_lo;
    while(i < w->c_len)
        w->coeff[i++] *= s_hi;

    u = fabs(w->coeff[0]);

    for(k = 1; k < w->c_len; k++)
        if((v = fabs(w->coeff[k])) > u)
        {
            w->c_max = k; u = v;
        }
}
