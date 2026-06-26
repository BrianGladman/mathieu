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
 Test driver for characteristic value functions. This driver tests the 
 solution obtained by solving the tri-diagional equations.
*/

#include <stdio.h>
#include "mathieu.h"

void main(void)
{
    double 	q, qlo, qhi, qi, *cv, *cvv;
    int 	i, r, rlo, rhi, nq, ty;
    hill_workspace w[1] = { 0, 0, 0};
    char   *st = { "A(q) (blanch) even         B(q) (blanch) odd        "
                   "A(q)   (hill) even         B(q)   (hill) odd" };

    printf("\nInitial and final q values and the increment? "); 
    scanf_s("%lf %lf %lf", &qlo, &qhi, &qi);

    printf("\nLowest and highest roots required? "); 
    scanf_s("%d %d", &rlo, &rhi);

    nq = (int)((qhi - qlo + 1.0e-10) / qi);
    ty = type_both;

#if 1   /* tables for r values  */ 
    cv  = (double*)malloc((2 * (rhi - rlo) + 3) * sizeof(double));
    cvv = (double*)malloc((2 * (rhi - rlo) + 3) * sizeof(double));

    for(i = 0, q = qlo; i <= nq; i++, q += qi)
    {
        mathieu_blanch_values(ty, q, rlo, rhi + 1, cv);
        mathieu_hill_values(ty, q, rlo, rhi + 1, cvv, w);
        printf("\n\nq = %8.2f     %s", q, st);

        for(r = rlo; r <= rhi ; ++r)
        {
            if(rlo)
            {
                printf("\n%8d %26.16e", r, cv[2 * (r - rlo) + 1]);
                if(r)
                    printf("%26.16e", cv[2 * (r - rlo)]);
                printf("%26.16e", cvv[2 * (r - rlo) + 1]);
                if(r)
                    printf("%26.16e", cvv[2 * (r - rlo)]);
            }
            else
            {
                printf("\n%8d %26.16e", r, cv[2 * (r - rlo)]);
                if(r)
                    printf("%26.16e", cv[2 * (r - rlo) - 1]);
                else
                    printf("%26s", " ");

                printf("%26.16e", cvv[2 * (r - rlo)]);
                if(r)
                {
                    printf("%26.16e", cvv[2 * (r - rlo) - 1]);
                }
                else
                    printf("%26s", " ");
            }
        }
    }
    free(cvv);
    free(cv);

#else   /* tables for q values  */

    for(r = rlo; r <= rhi ; ++r)
    {
        if(r)
        {
            printf("\nr = %8d  (odd)", r);

            for(i = 0, q = qlo; i <= nq; i++, q += qi)
            {
                mathieu_blanch_values(type_a, q, r, r + 1, &scv);
                printf("\n%12.6f %26.16e", q, scv);
            }
        }

        printf("\nr = %8d  (even)", r);

        for(i = 0, q = qlo; i <= nq; i++, q += qi)
        {
            mathieu_blanch_values(type_b, q, r, r + 1, &scv);
            printf("\n%12.6f %26.16e", q, scv);
        }
    }
#endif
    printf("\n");
}
