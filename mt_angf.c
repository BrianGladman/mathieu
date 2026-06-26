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

#include "mathieu.h"

/*
 The routine is a member of mathieu and computes periodic angular solutions 
 of Mathieu's equation, cem(x,q) and sem(x,q) and their derivatives.   The
 parameter q and the root values, r and a, are bound into mathieu itself so 
 that the function only requires the x value as a parameter.   The routine
 returns both the function and its derivative.
*/

void mathieu_angular(const double x, double *af, double *ad, mathieu_qr w[1])
{
    int 	k, sw, ki, xf;
    double	f, ff, s, u, v, s2x, c2x;

    if(!w->coeff)
        mathieu_expansion_coeffs(w);

    sw = (w->type == type_a ? 0 : 2)  + (w->r & 1);

    c2x = cos(2.0 * x); 
    s2x = sin(2.0 * x);

    switch(sw)
    {
        case 0:	u = 1.0; v = 0.0; s2x = -s2x; break;
        case 1: u = cos(x); v = sin(x); s2x = -s2x; break;
        case 2: u = sin(2.0 * x); v = cos(2.0 * x); break;
        case 3:	u = sin(x); v = cos(x); break;
    }

    f = 0.0; k = 0;
    do
    {
        ff = f; 
        f += w->coeff[k] * u; 
        s = u * c2x + v * s2x, v = v * c2x - u * s2x, u = s;
        xf = (k > w->c_max && f == ff ? xf + 1 : 0);
    }
    while
        (++k < w->c_len && xf < 4);

    *af = f;

    switch(sw)
    {
        case 0:	ki = 0; u = 0.0; v = -1.0; break;
        case 1: ki = 1; u = -sin(x); v = -cos(x); break;
        case 2: ki = 2; u = cos(2.0 * x); v = sin(2.0 * x); break;
        case 3:	ki = 1; u = cos(x); v = sin(x); break;
    }

    s2x = -s2x; f = 0.0; k = 0;
    do
    {
        ff = f; 
        f += ki * w->coeff[k] * u; 
        ki += 2; 
        s = u * c2x + v * s2x, v = v * c2x - u * s2x, u = s;
        xf = (k > w->c_max && f == ff ? xf + 1 : 0);
    }
    while
        (++k < w->c_len && xf < 4);

    *ad = f;

    return;
};
