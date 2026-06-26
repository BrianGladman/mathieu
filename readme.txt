
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

Implementation Details

This header file contains the definitions for the structures and classes that
maintains common data for each particular solution.

The code contains a self standing bessel function implementation that is needed
for the computation of the mathieu radial functions.

For each value of the pair (q, r), a structure is needed to hold the expansion 
coefficients used for the series for the angular and radial functions. This 
structure (or class) is then used in the computation of the angular and radial 
functions. The latter functions also need temporary storage for the bessel 
functions used to calculate the radial functions. 

All computation is coded in C but a C++ class definition is included for use
if necessaey.

Test Applications

1. tdtest

This outputs the Mathieu Characteristic values for a range of q and r values. The 
output is in the form of a table with the first column being the q value, the 
second and third colums being A(q) and B(q) computed using Blanch's method and
the the fourth and fifth columns being A(q) and B(q) using Hill's method.

2. testang

For each specified range of orders (r) and range of q values this application
outputs tables for ce(a, q, x), ce'(a,q,x), se(a,q,x) and se'(a,q,x) where 
the x value is in degrees.

3. testrad

For each specified range of orders (r) and range of q values this application outputs
two tables: one for Mc1(r, q, z), Mc1'(r, q, z), Mc2(r, q, z) and Mc2'(r, q, z) and a
table for Ms1(r, q, z), Ms1'(r, q, z), Ms2(r, q, z) and Ms2'(r, q, z).
