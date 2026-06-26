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

 This routine calculates a sequence of characteristic values for Mathieu's
 equation for a given function type (ty) and parameter (q) by obtaining the
 eigenvalues of the tri-diagonal matrix for the expansion coefficients.  It
 produces a sequence of characteristic values including solutions that are
 of period pi and 2 * pi by solving the four distinct tri-diagonal equations
 and interleaving the roots obtained in the output array. The parameters are:
  
   Input	ty		function type (type_a = 1, type_b = 2, type_both = 3)
  			q		the parameter value
  			r_min	the lowest root required
  			r_max	one above the highest root required
  
   Output	cv[]	array of characteristic values r_min <= r <= r_max - 1
  					in increasing order of value
  
 The output order for roots is a0, b1, a1, b2, a2,...  Note that r_min cannot
 be 1 if ty is 2 since the lowest root in in this case is for r = 1.  If an 
 error occurs the funtion returns a non zero value, otherwise it returns zero.  

 This routine calculates the characteristic values for Mathieu functions by
 computing the eignevalues of the Hill determinant The hill determinant are 
 of four possible forms as follows:
  
 Type 0: Even solution with period pi (a0, a2, a4, ... cos):
 
 |     0^2 - a   sqrt(2).q          0       ... |
 |   sqrt(2).q     2^2 - a          q       ... |
 |           0           q    4^2 - a     q ... |
 | .....                                        |
 
 Type 1:  Odd solution with period 2.pi (b1, b3, b5, ... sin):
 
 | 1^2 - q - a           q          0       ... |
 |           q     3^2 - a          q       ... |
 |           0           q    5^2 - a     q ... |
 | .....                                        |
 
 Type 2: Even solution with period 2.pi (a1, a3, a5, ... cos):
 
 | 1^2 + q - a           q          0       ... |
 |           q     3^2 - a          q       ... |
 |           0           q    5^2 - a     q ... |
 | .....                                        |
 
 Type 3: Odd solution with period pi (b2, b4, b6, ... sin):
 
 | 2^2 - q - a           q          0       ... |
 |           q     4^2 - a          q       ... |
 |           0           q    6^2 - a     q ... |
 | .....                                        |
 
 The numerical order of the characetristic values for
 q > 0 is a0 < b1 < a1 < b2 < a2 < b3 .....  For q < 0
 the order is: a0 < a1 < b1 < b2 ..... 
*/

#include <malloc.h>
#include "mathieu.h"

int alloc_hill_workspace(size_t nm, hill_workspace *h)
{
    if(h->dp = (double*)malloc(nm * sizeof(double)))
        if(h->ep = (double*)malloc(nm * sizeof(double)))
        {
            h->nm = nm;
            return 0;
        }
        else
            free(h->dp);
    h->nm = 0;
    return EXIT_ERROR;
}

void free_hill_workspace(hill_workspace *h)
{
    if(h->ep)
        free(h->ep);
    if(h->ep)
        free(h->dp);
    h->dp = h->ep = 0;
    h->nm = 0;
}

int mathieu_hill_values(int ty, double q, int r_min, int r_max, 
                                            double cv[], hill_workspace *h)
{	int		i, r, r0, nm, m_type;

    if(r_max <= r_min || ty == type_b && r_min == 0)
        return 0;

    /* we need to expand the matrix for more terms than are asked for	*/
    /* to obtain high accuracy characteristic values. This value is		*/
    /* intended to achieve 14 decimal digit accuracy					*/
    nm = hill_determinant_size(q, r_max);
    alloc_hill_workspace(nm, h);

    /* scan the four distinct forms of solution	*/
    for(m_type = 0; m_type < 4; ++m_type)			
    {	
        /* solution type is m_type & 1 - skip if not needed		*/
        if(!(m_type & 1) && !(ty & type_a) || (m_type & 1) && !(ty & type_b))
            continue;

        /* r0 is the lowest root for this form of solution	*/
        r0 = (m_type < 2 ? m_type : m_type - 1);
                                            
        /* skip even or odd root solutions if not needed	*/
        if(r_max == r_min + 1 && (r0 & 1) != (r_min & 1))
            continue;
        
        /* set up the tridaigonal matrix					*/
        for(i = 0, r = r0; i < nm; ++i, r += 2)
        {
            h->ep[i] = q; 
            h->dp[i] = (double)r * r;
        }

        if(r0 & 1)
            h->dp[0] += ((m_type & 1) ? -q : q);
        else if(!r0)
            h->ep[0] = sqrt(2.0) * q;

        if(q)
            if(tridiag_matrix_eigenvalues(h->dp, h->ep, nm, 0, 0))
                return -1;

        /* output required characteristic values in order	*/
        /* of increasing value								*/
        
        if(ty == (type_a|type_b))	/* both types needed	*/
        {
            i = ((r_min ? (r_min & 1 ? 3 : 1) : 0) + m_type) % 4;
            r = r_min  + (r_min ? i / 2 : r0);
        }
        else		/* either but not both a and b required	*/
        {
            i = ((r_min & 1) == (r0 & 1) ? 0 : 1); 
            r = r_min + i;
        }

        while(r < r_max)	/* move values to output array	*/
        {
            cv[i] = h->dp[(r - r0) / 2]; 
            r += 2; 
            i += (ty == (type_a|type_b) ? 4 : 2);
        }
    }

    free_hill_workspace(h);
    return EXIT_OK;
}
