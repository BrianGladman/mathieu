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
 
 Calculate one accurate characteristic value using the secant root finding 
 approach if possible but reverting to the use of the Hill Determinant if 
 no starting estimate is available for the specified type, q and r values. 

   input   ty	type (type_a = 1 or type_b = 2)
  			q	parameter
  			r	order
  			a	characteristic value estimate

   output	a	refined characteristic value

 The starting tolerance is 10^-15 but this is relaxed if convergence fails 
 to see if a root can be found. However the tolerance is not reduced below 
 10^-12.
*/

#include <malloc.h>
#include "mathieu.h"

static double mathieu_characteristic_function(int ty, double q, int r, double a)
{	int 	m, m0, ms, mx, mc;
    double 	g1, g2, t, a0, b0, a1, b1, v;

    ms = r / 2 + (int)(3.2 * pow(q, 1.0 / 3.0));
    m0 = 1; 

    /* set starting values for forward iteration			*/

    if(r & 1)	/* start at g3, compute g3, g5,... upwards	*/
        g2 = (a - 1.0) / q + (ty == type_b ? 1.0 : -1.0);
    else if(ty == type_b)   /* start at g4, compute g6,g8,.	*/
    {
        g2 = (a - 4.0) / q; 
        m0 = 2;
    }	
    else		/* start at g2, compute g4,g6,...			*/
        g2 = a / (2.0 * q);

    /* compute g's in upward direction until |g| < 1		*/
    for (m = m0; m < ms && fabs(g2) >= 1.0; ++m)
    {
        t = m + m + (r & 1); 
        g2 = (a - t * t) / q - 1.0 / g2; 
    }
    
    mc = m;
    mx = (ms > m + 10 ? ms : m + 10); 

    /* compute (negative of) the continued fraction tail	*/
    a0 = b1 = 1.0; a1 = b0 = 0.0; m = mx;
    do
    {	t = m + m + (r & 1); 
        t = (a - t * t) / q;
        
        v = t * a1 - a0; 
        a0 = a1; 
        a1 = v;
        
        v = t * b1 - b0; 
        b0 = b1; 
        b1 = v;

    }
    while
        (fabs(a0 * b1 - a1 * b0) > 1.0e-14 * fabs(a1 * b0) && ++m < ms + 20);
    
    /* compute g's in downward direction using tail value   */

    g1 = -a1 / b1;

    for (m = mx - 1; m >= mc; --m)
    {
        t = m + m + (r & 1); 
        g1 = 1.0 / ((a - t * t) / q - g1);
    }

    return g2 - g1;
}

int mathieu_characteristic_root(int ty, double q, int r, double *a)
{	double	x0, x1, f0, f1, inc, etol, elim;
    int		i, it;

    if(q == 0.0)
        return EXIT_OK;

    etol = 1.0e-15;	/* set the initial tolerance	*/

    /* iterate the secant method with a reduced		*/
    /* accuracy if no result in 25 cycles			*/
    for(i = 0; i < 8; ++i)
    {
        it = 0;		/* root finder iteration count	*/
        elim = etol * (fabs(*a) + 0.01);

        /* find two initial function values for the	*/
        /* secant root finding method				*/
        x1 = *a;
        x0 = 1.002 * x1;
        f0 = mathieu_characteristic_function(ty, q, r, x0);
        f1 = mathieu_characteristic_function(ty, q, r, x1);

        /* try to find root in less than 25 cycles	*/
        do
        {	/* estimate root from existing values	*/
            inc = f1 * (x1 - x0) / (f1 - f0);
            x0 = x1; f0 = f1; x1 -= inc; 
            f1 = mathieu_characteristic_function(ty, q, r, x1);
        }
        while
            (++it < 25 && fabs(inc) > elim && f1 != f0);

        *a = x1;  

        if(it < 25) 
            return EXIT_OK;

        etol *= 3.0;
    }

    return EXIT_ERROR;
}

/* This routine calculates the characteristic values for Mathieu
   functions using the approach developed by Blanch 
*/

int mathieu_blanch_values(int ty, double q, int r_min, int r_max, double cv[])
{
    int i, r, r0, nm, m_type;
    double cav, cbv, cva[4], cvb[4];

    if(r_max <= r_min || ty == 0 || ty & type_both == type_b && r_min == 0) 
        return -1;

    i = r_max - r_min;
    if(ty & type_both == type_both)
      i = i + i - (r_min ? 0 : 1);

    r = 0;
    r0 = 10;
    if(r_min > 10) 
    {
        if(r_min * r_min < (int)(0.2e0 * q) - 6 || r_min * r_min > (int)(4.0e0 * q) + 4)
            r = r_min;
        else
        {
            r = sqrt(0.2e0 * q) - 6;
            r = r < 6 ? 6 : r;
        }
        r0 = r + 4;
    }

    for( ; r < r_max ; ++r )
    {
        if(ty & type_a)
        {
            if(r < r0)
                cav = mathieu_characteristic_approx(type_a, q, r);
            else
            {
                cav = 4.0 * cva[0] - 6.0 * cva[1] + 4.0 * cva[2] - cva[3];
                cav = cav > cva[0] ? cav : cva[0];
            }
            cva[3] = cva[2]; cva[2] = cva[1]; cva[1] = cva[0];
            if(mathieu_characteristic_root(type_a, q, r, &cav) == EXIT_ERROR)
                return -1;
            cva[0] = cav;
        }

        if((ty & type_b) && r != 0)
        {
            if(r < r0)
                cbv = mathieu_characteristic_approx(type_b, q, r);
            else
            {
                cbv = 4.0 * cvb[0] - 6.0 * cvb[1] + 4.0 * cvb[2] - cvb[3];
                cbv = cbv > cvb[0] ? cbv : cvb[0];
            }
            cvb[3] = cvb[2]; cvb[2] = cvb[1]; cvb[1] = cvb[0];
            if(mathieu_characteristic_root(type_b, q, r, &cbv) == EXIT_ERROR) 
                return - 1;
            cvb[0] = cbv;
        }

        if(r >= r_min)
        {
            if((ty & type_both) == type_both)
            {
                if(r_min == 0)
                {
                    cv[2 * r] = cav;
                    if(r != 0) 
                        cv[2 * r - 1] = cbv;
                }
                else
                {
                    cv[2 * (r - r_min)] = cbv;
                    cv[2 * (r - r_min) + 1] = cav;
                }
            }
            else if((ty & type_both) == type_a)
                cv[r - r_min] = cav;
            else
                cv[r - r_min] = cbv;
        }
    }
    return 0;
}
