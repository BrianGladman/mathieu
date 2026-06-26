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

 The set of programs described here compute the characteristic values for 
 periodic solutions of Mathieu's equation together with the related angular 
 and radial functions.  

 The parameters used here are for Mathieu's equation in the form:

	      d^2y / dx^2 + (a - 2q cos(2x)) y = 0

as described here: https://mathworld.wolfram.com/MathieuFunction.html

 This suite of programs use the following common parameters:
  
    ty   function type (see below)
     q	 function parameter
     r	 function order

    function behaviour:		    at: 0  pi/2  period coefficients
    ty: 1  r: 0,2,4,.. cem(x,q)  even  even      pi    a[0,2,..]
    ty: 2  r: 1,3,5,.. sem(x,q)   odd  even    2.pi    b[1,3,..]
    ty: 1  r: 1,3,5,.. cem(x,q)  even   odd    2.pi    a[1,3,..]  
    ty: 2  r: 2,4,6,.. sem(x,q)   odd   odd      pi    b[2,4,..]

 This header file contains the definitions for the structures and classes that
 maintains common data for each particular solution.

 For each value of the pair (q, r), a structure is needed to hold the expansion 
 coefficients used for the series for the angular and radial functions. This 
 structure (or class) is then used in the computation of the angular and radial 
 functions. The latter functions also need temporary storage for the bessel 
 functions used to calculate the radial functions. Most of the computation is 
 done in routines lying outside the C++ class definition so they can be used 
 in C.
*/

#ifndef MATHIEU_H
#define MATHIEU_H

#define _USE_MATH_DEFINES

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#define NORMALISED

#if defined(__cplusplus)
extern "C"
{
#endif

#define NaN(x)	((double)(0xffff0000+(x)))
#define mathieu_no_estimate	NaN(1)

#define EXIT_OK     0
#define EXIT_ERROR -1

typedef struct
{
    size_t  nm;
    double *dp;
    double *ep;
} hill_workspace;

typedef enum { type_a = 1, type_b = 2, type_both = 3 }  mt_type;

typedef struct	/* structure to store bessel functions used in computing    */
{				/* radial solutions to Mathieu's equation                   */
    double	b_arg;
    int		b_nof;
    double *b_jnp;
    double *b_ynp;
} bf_str;

typedef struct
{
    mt_type	type;	/* solution type (see above)            */
    int		r;	    /* root (r) number                      */
    double	q;	    /* parameter (q) value                  */
    double	cv;	    /* root (a or b) value                  */
    double *coeff;	/* pointer to coefficients              */
    int		c_len;	/* number of significant coefficients   */
    int		c_max;	/* index of the maximum coefficient     */
    bf_str	bf1;	/* storage for bessel functions for u1  */
    bf_str	bf2;	/* storage for bessel functions for u2  */
    hill_workspace h;
} mathieu_qr;

#define hill_determinant_size(q, r_max) ((r_max) / 2 \
                        + (int)(3.2 * pow((q), 1.0 / 3.0)) + 20)

int alloc_hill_workspace(size_t n, hill_workspace *h);
void free_hill_workspace(hill_workspace *h);

int tridiag_matrix_eigenvalues(double d[], double e[], 
                                    int n, int ev_req, double z[]);

int mathieu_characteristic_root(int ty, double q, int r, double *a);

double mathieu_characteristic_approx(int ty, double q, int r);

int mathieu_blanch_values(int ty, double q, int r_min, int r_max, double cv[]);

int mathieu_hill_values(int ty, double q, int r_min, int r_max, double cv[], 
                                    hill_workspace *h);

double mathieu_init_qr(const mt_type ty, const int r, 
                                    const double q, mathieu_qr w[1]); 

void mathieu_free_qr(mathieu_qr w[1]);

void mathieu_expansion_coeffs(mathieu_qr w[1]);

void mathieu_angular(const double x, double *af, double *ad, mathieu_qr w[1]);

void mathieu_radial(int kc, const double x, double *rf1, double *rf2, 
                                    double *rd1, double *rd2, mathieu_qr w[1]);

#if defined(__cplusplus)
}

class	mathieu
{
    private:
        mathieu_qr  w;

    public:

        mathieu(const mt_type ty, const int root, const double q)
        {
            mathieu_init_qr(ty, root, q, &w); 
        }
       
        ~mathieu(void)
        {
            mathieu_free_qr(&w); 
        }
        
        double characteristic_value(void) { return w.cv; };

        void angular(const double x, double *af, double *afd)
        {
            if(!w.coeff)
                mathieu_expansion_coeffs(&w);

            mathieu_angular(x, af, afd, &w);
        }

        void radial(const int kc, const double x, double *f1r, 
                                double *f2r, double *d1r, double *d2r)
        {
            if(!w.coeff)
                mathieu_expansion_coeffs(&w);
            
            mathieu_radial(kc, x, f1r, f2r, d1r, d2r, &w);
        }
};

#endif
#endif
