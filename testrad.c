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

#include <stdio.h>
#include <time.h>
#include "mathieu.h"

void main(void)
{
    double 	q, qlo, qhi, qi, x, xlo, xhi, xi, a, b, f1r, f2r, d1r, d2r;
    int 	i, j, rb, rlo, rhi, nq, nx;

    printf("\nLowest and highest roots required? ");
    scanf("%d %d", &rlo, &rhi);

    printf("\nInitial and final q values and the increment? ");
    scanf("%lf %lf %lf", &qlo, &qhi, &qi);

    printf("\nInitial and final x values and the increment? ");
    scanf("%lf %lf %lf", &xlo, &xhi, &xi);

    clock_t start = clock();
    nq = (int)((qhi - qlo + 1.0e-10) / qi);
    nx = (int)((xhi - xlo + 1.0e-10) / xi);

    for(rb = rlo; rb <= rhi ; ++rb)
    {
        for(i = 0, q = qlo; i <= nq; i++, q += qi)
        {
            mathieu_qr s0[1], s1[1];

            a = mathieu_init_qr(type_a, rb, q, s0);
            if(rb)
                b = mathieu_init_qr(type_b, rb, q, s1);

            printf("\norder (r) = %4d, q = %12.6f (a = %26.16e", rb, q, a);
            if(rb)
                printf(", b = %26.16e)", b);
            else
                printf(")");

            printf("\n         x                 Mc1(r, q)                Mc1'(r, q)"
                   "                  Mc2(r,q)                 Mc2'(r,q)");

            for(j = 0, x = xlo; j <= nx; j++, x += xi)
            {
                mathieu_radial(15, x, &f1r, &f2r, &d1r, &d2r, s0);
                printf("\n%10.2f %25.16e %25.16e %25.16e %25.16e", x, 
                    f1r, d1r, f2r, d2r);
            }
            
            if(rb)
            {
                printf("\n         x                 Ms1(r, q)                Ms1'(r, q)"
                       "                  Ms2(r,q)                 Ms2'(r,q)");

                for(j = 0, x = xlo; j <= nx; j++, x += xi)
                {
                    mathieu_radial(15, x, &f1r, &f2r, &d1r, &d2r, s1);
                    printf("\n%10.2f %25.16e %25.16e %25.16e %25.16e", x, 
                        f1r, d1r, f2r, d2r);
                }
            }
            printf("\n");
            mathieu_free_qr(s0);
            if(rb)
                mathieu_free_qr(s1);
        }
    }
    printf("\n");
    start = clock() - start;
    printf("\nTime: %2.1f seconds\n\n\n", ((double)start) / (1.0 * CLOCKS_PER_SEC));

}
