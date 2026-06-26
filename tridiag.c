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
 This subroutine finds the eigenvalues and eigenvectors of a symmetric 
 tridiagonal matrix by the implicit ql method.
 
 This subroutine is a translation of the algol procedure in Numerical 
 Mathematics 12, p377-383 (1968) by Martin and Wilkinson, as modified 
 in numerical mathematics 15, p 450(1970) by Dubrulle.   Handbook for 
 Automatic Computation, vol.ii-linear algebra, p241-248(1971).
 
 Inputs:

  		d		array containing the diagonal elements of the input 
    			matrix in elements 0..n-1
  		e		array containing the subdiagonal elements of the input 
 	    		matrix in elements 0..n-2; e[n-1] is not used.
  		n		the order of the matrix.
  		v_req	is true if the eigenvectors are required, false otherwise
  		z		is the array for eignevector output; the eigenvectors 0, 
    			1, 2,.. are placed in elements 0..n-1, n..2n-1, 2n..3n-1
 	    		as a single dimension array; if the eignevectors are to 
 		    	be computed z must contain the identity matrix on input
 Outputs:
 
 		d		contains the eigenvalues in ascending order.  if an
 				error exit is made, the eigenvalues are correct but
 				unordered for indices 0,1,...,err - 1.
 
 		e		destroyed.
 
 		z		contains orthonormal eigenvectors of the symmetric
 				tridiagonal matrix.  if an error exit is made, z 
 				contains the eigenvectors associated with the stored
 				eigenvalues.
 
 Return (err)
 
 		0		for normal return,
 
        i	    if the i-th eigenvalue has not been determined after 
 				30 iterations (i = 1,2...)
*/

#include <stdlib.h>
#include <math.h>

#define MAXIT	30

static int double_cmp(const void *a, const void*b)
{
    return (*(double*)a < *(double*)b ? -1 :
            *(double*)a > *(double*)b ?  1 : 0);
}

static double diag(const double x, const double y)
{	double	r;

    if(fabs(x) > fabs(y))
        return (r = y / x, fabs(x) * sqrt(1.0 + r * r));
    else if(y != 0.0)
        return (r = x / y, fabs(y) * sqrt(1.0 + r * r));
    else
        return 0.0;
}

int tridiag_matrix_eigenvalues(double d[], double e[], int n, int ev_req, double z[])
{	int		i, j, k, it, ev_no;
    double	b, c, f, g, p, r, s;

    for(ev_no = 0; ev_no < n; ++ev_no)		/* loop for each eigenvalue		*/
    {
        for(it = 0; it < 30; ++it)			/* iterate for next eigenvalue	*/
        {
            for(i = ev_no; i < n - 1; ++i)	/* find an insignificant		*/
            {								/* off-diagonal element			*/ 
                p = fabs(d[i]) + fabs(d[i + 1]);
                if(p + fabs(e[i]) == p) 
                    break;
            }

            if(i == ev_no) 
                break;			            /* exit if an eigenvalue found	*/
            
            /* form an implicit shift										*/
            g = (d[ev_no + 1] - d[ev_no]) / (2.0 * e[ev_no]);	
            g += (g < 0.0 ? -1.0 : 1.0) * diag(g, 1.0);	
            g = d[i] - d[ev_no] + e[ev_no] / g;	
            s = c = 1.0; 
            p = 0.0; 

            /* perform a plane rotation followed by Givens rotations to		*/
            /* restore tridiagonal form										*/
            for(j = i; j > ev_no; --j)		
            {
                f = s * e[j - 1]; b = c * e[j - 1];
                e[j] = r = diag(f, g);
        
                if(r == 0.0) 
                    break;
        
                s = f / r; c = g / r; 
                g = d[j] - p;
                r = (d[j - 1] - g) * s + 2.0 * b * c;
                p = s * r; 
                d[j] = g + p; g = c * r - b;

                if(ev_req)					/* if eignevalues are required	*/
                    for(k = (j - 1) * n; k < j * n; ++k) 
                    {
                        b = z[k]; f = z[k + n];
                        z[k] = c * b - s * f; 
                        z[k + n] = s * b + c * f;
                    }
            }

            if(j == ev_no)
                e[j] = g;

            d[j] -= p; e[i] = 0.0;
        }

        if(it == MAXIT)					/* if the process fails to converge	*/
            return ev_no + 1;
    }

    /* sort eigenvalues/vectors in numerical order	*/
    qsort(d, n, sizeof(double), double_cmp);
    return 0;
}
