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
 Power series expansions for characteristic values for low q - the 
 value of a[r] and b[r] are given by:

       r * r + even * q ^ 2 + (ty == type_b ? -1 : 1) * odd * q ^ r
 
 where:
 
       even (ty == 1) = sum{ce[n] * (q * q) ^ n}
       odd  (ty == 2) = sum(co[n] * (q * q) ^ n}
 
 The range of the power series is 0 <= q <= r + 1. Note that the
 so called 'odd' component is only odd when r is odd. These tables
 were produced by using the program to generate characteristic 
 values starting with estimates from power series expansions 
*/

#include "mathieu.h"

static float ce[5][4] = 
{
  {-4.999953e-01f, 5.453104e-02f,-1.174431e-02f, 2.074040e-03f }, //  0
  {-1.250000e-01f,-6.511377e-04f, 8.331945e-05f,-2.549378e-06f }, //  1
  { 1.656757e-01f,-2.290380e-02f, 2.354393e-03f,-1.060710e-04f }, //  2
  { 6.252267e-02f, 5.907108e-04f,-7.064678e-05f, 1.592599e-06f }, //  3
  { 3.333187e-02f, 6.900294e-05f, 4.230433e-07f,-2.731119e-08f }, //  4
};

static float co[4][4] = 
{
  { 1.000001e+00f,-1.563419e-02f, 3.103979e-04f, 6.153545e-07f }, //  1
  { 2.490090e-01f,-2.326540e-02f, 2.357969e-03f,-1.061075e-04f }, //  2
  { 1.563225e-02f,-3.208815e-04f,-2.889445e-07f, 1.448925e-07f }, //  3
  { 4.340139e-04f,-2.874380e-06f,-8.896640e-09f, 1.334570e-10f }, //  4
};

/*
 Rational function approximations for lower characteristic values.  The
 approximation is given by T/B where:
  
        T = c[0] + q * (c[1] + q * (c[2] + q * c[3]))
        B =  1.0 + q * (c[4] + q * (c[5] + q * (c[6] + q * c[7])))
  
 for ql <= q <= qh and the coefficients 'c' from the following set as
 determined by r and ty for 0 <= r <= 12 and ty = 1, 2.   Again these 
 expansions were obtained by using parts of the program opeating with
 a rational approximation routine. 
*/

static float coef[39][7] = 
{
/*  r =  0   ty =  1   ql =     1.00   qh =    19.00   error = 4.56e-02		*/
  { -9.572504e-03f,  2.181183e-01f, -8.957264e-01f,  5.443359e-01f,
    -1.995345e-03f,  2.746719e-05f, -3.059519e-08f },

/*  r =  1   ty =  1   ql =     2.00   qh =    20.00   error = 3.05e-03		*/

  { -4.781774e-01f,  1.558679e+00f, -2.355274e-01f,  1.018049e-01f,
     8.093620e-03f, -4.187466e-04f,  7.357119e-06f },

/*  r =  1   ty =  2   ql =     2.00   qh =     8.00    error = 1.02e-06	*/

  { -5.330372e-03f, -9.885717e-01f, -3.001315e-01f,  1.644902e-01f,
     7.344259e-04f, -6.590539e-05f,  1.920919e-06f },

/*  r =  2   ty =  1   ql =     3.00   qh =    23.00   error = 3.99e-03		*/
  { -8.416122e-01f,  1.026750e+00f, -5.912349e-02f, -8.835663e-02f,
     1.272722e-02f, -5.232156e-04f,  7.827410e-06f },

/*  r =  2   ty =  2   ql =     3.00   qh =    11.00    error = 3.98e-06	*/

  { -8.082837e-03f,  9.241842e-03f, -8.734052e-02f,  8.683684e-03f,
     4.348346e-03f, -2.063813e-04f,  4.000632e-06f },

/*  r =  3   ty =  1   ql =     4.00   qh =    28.00   error = 6.42e-03		*/

  { -1.097120e+00f,  6.199764e-01f, -1.866312e-02f, -1.233746e-01f,
     1.037156e-02f, -3.517304e-04f,  4.579231e-06f },

/*  r =  3   ty =  2   ql =     4.00   qh =    16.00    error = 9.42e-05	*/

  { -1.904479e-01f,  1.972414e-01f, -2.483147e-02f, -8.544548e-02f,
     8.200276e-03f, -3.258739e-04f,  5.181634e-06f },

/*  r =  4   ty =  1   ql =     5.00   qh =    35.00   error = 2.08e-02		*/

  { -6.559897e-01f,  2.439298e-01f,  1.547649e-03f, -1.241765e-01f,
     8.287584e-03f, -2.349010e-04f,  2.912091e-06f },

/*  r =  4   ty =  2   ql =     5.00   qh =    23.00    error = 7.46e-04	*/

  { -2.994338e-01f,  1.857419e-01f, -9.077080e-03f, -1.089814e-01f,
     7.867404e-03f, -2.597814e-04f,  3.356968e-06f },

/*  r =  5   ty =  1   ql =     3.50   qh =    44.00   error = 1.39e-02		*/

  {  4.124948e-02f,  2.371988e-04f,  1.245009e-02f, -8.872630e-02f,
     4.395964e-03f, -8.720352e-05f,  9.057343e-07f },

/*  r =  5   ty =  2   ql =     4.50   qh =    32.00    error = 5.38e-03	*/

  { -2.375640e-01f,  1.239652e-01f, -3.341063e-03f, -1.119682e-01f,
     6.626589e-03f, -1.836208e-04f,  1.994051e-06f },

/*  r =  6   ty =  1   ql =     9.00   qh =    55.00   error = 1.04e-02		*/

  {  3.033573e-01f, -2.665733e-02f,  7.695309e-03f, -6.716398e-02f,
     2.281432e-03f, -3.137019e-05f,  2.201043e-07f },

/*  r =  6   ty =  2   ql =    10.00   qh =    43.00    error = 2.29e-03	*/

  { -4.453535e-01f,  1.342611e-01f, -2.092198e-03f, -7.755135e-02f,
     3.171118e-03f, -6.104305e-05f,  4.786959e-07f },

/*  r =  7   ty =  1   ql =    15.50   qh =    68.00   error = 9.30e-03		*/

  {  6.022800e-01f, -2.964852e-02f,  4.227882e-03f, -5.532242e-02f,
     1.411226e-03f, -1.501976e-05f,  7.388252e-08f },

/*  r =  7   ty =  2   ql =    16.50   qh =    56.00    error = 1.31e-03	*/

  { -6.016477e-01f,  1.273952e-01f, -1.138846e-03f, -5.828675e-02f,
     1.760176e-03f, -2.485785e-05f,  1.466472e-07f },

/*  r =  8   ty =  1   ql =    23.00   qh =    83.00   error = 9.01e-03		*/

  {  9.787534e-01f, -3.099186e-02f,  2.400737e-03f, -4.660754e-02f,
     9.439266e-04f, -8.169672e-06f,  3.036888e-08f },

/*  r =  8   ty =  2   ql =    24.00   qh =    71.00    error = 8.63e-04	*/

  { -6.732349e-01f,  1.131661e-01f, -5.199122e-04f, -4.575569e-02f,
     1.062658e-03f, -1.143680e-05f,  5.224866e-08f },

/*  r =  9   ty =  1   ql =    31.50   qh =   100.00   error = 9.98e-03		*/

  {  1.464729e+00f, -3.499184e-02f,  1.486529e-03f, -3.953768e-02f,
     6.545760e-04f, -4.697334e-06f,  1.388026e-08f },

/*  r =  9   ty =  2   ql =    32.50   qh =    88.00    error = 5.75e-04	*/

  { -6.395463e-01f,  9.638234e-02f, -1.526826e-04f, -3.703403e-02f,
     6.808179e-04f, -5.749646e-06f,  2.077988e-08f },

/*  r = 10   ty =  1   ql =    41.00   qh =   119.00   error = 1.21e-02		*/

  {  2.114565e+00f, -4.302238e-02f,  1.035079e-03f, -3.381067e-02f,
     4.675469e-04f, -2.826382e-06f,  6.877568e-09f },

/*  r = 10   ty =  2   ql =    42.00   qh =   107.00    error = 3.53e-04	*/

  { -4.874417e-01f,  7.933522e-02f,  5.001251e-05f, -3.072697e-02f,
     4.578147e-04f, -3.110670e-06f,  9.057229e-09f },

/*  r = 11   ty =  1   ql =    51.50   qh =   140.00   error = 1.60e-02		*/

  {  2.978135e+00f, -5.470297e-02f,  8.157653e-04f, -2.911458e-02f,
     3.413433e-04f, -1.757393e-06f,  3.601725e-09f },

/*  r = 11   ty =  2   ql =    52.50   qh =   128.00    error = 1.79e-04	*/

  { -1.639266e-01f,  6.100362e-02f,  1.783854e-04f, -2.605267e-02f,
     3.218543e-04f, -1.802452e-06f,  4.321228e-09f },

/*  r = 12   ty =  1   ql =    63.00   qh =   163.00   error = 2.08e-02		*/

  {  3.997476e+00f, -6.579732e-02f,  6.739983e-04f, -2.533854e-02f,
     2.555058e-04f, -1.135118e-06f,  1.992598e-09f },

/*  r = 12   ty =  2   ql =    64.00   qh =   151.00    error = 4.64e-04	*/

  {  2.411235e-01f,  4.640926e-02f,  2.157775e-04f, -2.243795e-02f,
     2.332651e-04f, -1.094795e-06f,  2.181775e-09f },

/*  r = 13   ty =  1   ql =    75.50   qh =   188.00   error = 2.88e-02		*/

  {  5.217540e+00f, -7.741978e-02f,  5.880204e-04f, -2.220874e-02f,
     1.945417e-04f, -7.521728e-07f,  1.144537e-09f },

/*  r = 13   ty =  2   ql =    76.50   qh =   176.00    error = 9.34e-04	*/

  {  7.743380e-01f,  3.301789e-02f,  2.266893e-04f, -1.959277e-02f,
     1.741259e-04f, -6.966677e-07f,  1.172470e-09f },

/*  r = 14   ty =  1   ql =    89.00   qh =   215.00   error = 3.89e-02		*/

  {  6.593237e+00f, -8.792469e-02f,  5.234550e-04f, -1.960669e-02f,
     1.505550e-04f, -5.109017e-07f,  6.805282e-10f },

/*  r = 14   ty =  2   ql =    90.00   qh =   203.00    error = 1.55e-03	*/

  {  1.497652e+00f,  1.904010e-02f,  2.395547e-04f, -1.728428e-02f,
     1.331573e-04f, -4.607459e-07f,  6.663622e-10f },

/*  r = 15   ty =  1   ql =   103.50   qh =   244.00   error = 5.17e-02		*/

  {  8.117087e+00f, -9.729289e-02f,  4.718434e-04f, -1.742290e-02f,
     1.182045e-04f, -3.547045e-07f,  4.170509e-10f },

/*  r = 15   ty =  2   ql =   104.50   qh =   232.00    error = 2.36e-03	*/

  {  2.383969e+00f,  5.786177e-03f,  2.446338e-04f, -1.537950e-02f,
     1.038274e-04f, -3.142558e-07f,  3.953774e-10f },

/*  r = 16   ty =  1   ql =   119.00   qh =   275.00   error = 6.72e-02		*/

  {  9.784390e+00f, -1.056141e-01f,  4.288956e-04f, -1.557457e-02f,
     9.400965e-05f, -2.511344e-07f,  2.625349e-10f },

/*  r = 16   ty =  2   ql =   120.00   qh =   263.00    error = 3.43e-03	*/

  {  3.436906e+00f, -6.741282e-03f,  2.450750e-04f, -1.378236e-02f,
     8.228247e-05f, -2.199199e-07f,  2.432106e-10f },

/*  r = 17   ty =  1   ql =   135.50   qh =   308.00   error = 8.57e-02		*/

  {  1.158976e+01f, -1.129747e-01f,  3.920885e-04f, -1.399819e-02f,
     7.564344e-05f, -1.809810e-07f,  1.692977e-10f },

/*  r = 17   ty =  2   ql =   136.50   qh =   296.00    error = 4.79e-03	*/

  {  4.658811e+00f, -1.854911e-02f,  2.425727e-04f, -1.242590e-02f,
     6.611329e-05f, -1.572874e-07f,  1.542488e-10f },

/*  r = 18   ty =  1   ql =   153.00   qh =   343.00   error = 1.07e-01		*/

  {  1.352735e+01f, -1.194564e-01f,  3.598803e-04f, -1.264425e-02f,
     6.151300e-05f, -1.325381e-07f,  1.115776e-10f },

/*  r = 18   ty =  2   ql =   154.00   qh =   331.00    error = 6.64e-03	*/

  {  6.049977e+00f, -2.963509e-02f,  2.381131e-04f, -1.126105e-02f,
     5.375258e-05f, -1.146045e-07f,  1.004153e-10f },

/*  r = 19   ty =  1   ql =   171.50   qh =   380.00   error = 1.39e-01		*/

  {  1.557423e+01f, -1.249093e-01f,  3.307836e-04f, -1.146912e-02f,
     5.047173e-05f, -9.839666e-08f,  7.491144e-11f },

/*  r = 19   ty =  2   ql =   172.50   qh =   368.00    error = 9.06e-03	*/

  {  7.727912e+00f, -4.167075e-02f,  2.392508e-04f, -1.024494e-02f,
     4.412932e-05f, -8.478660e-08f,  6.690970e-11f }
};

/*
 Compute the value of a rational function with coeficients for the 
 numerator powers 0..L_TOP and the denominator powers 0..L_BOT. The
 coefficients for the numerator are in coef[0..L_TOP] and those for 
 the denominator ae in coef[L_TOP+1..L_TOP+L_BOT].  for the latter
 the first (constant) coefficient is assumed to be 1 (not stored). 
*/

#define L_TOP	2
#define L_BOT	4

float ratval(const float x, const float *coef)
{	float	top = coef[L_TOP], bottom = 0.0f;
    int		i;

    for(i = L_TOP - 1; i >= 0; --i)
        top = x * top + coef[i];

    for(i = L_TOP + L_BOT; i >= L_TOP + 1; --i)

        bottom = x * (bottom + coef[i]);

    return top / (1.0f + bottom);
}

/*
 Use power series above to obtain an estimate when r is low and q is
 in the range 0 <= q <= r + 1.0, otherwise use the low q approximation 
 for moderate q and large r. This provides a good initial estimate for 
 the blanch root finding procedure r < floor(sqrt(q / 5)).
*/

static float loq_est(const int ty, const int r, const double q)
{	float	r2 = (float)(r * r), t2 = (float)(q * q), ve, vo, *c;
    int		i;

    ve = vo = 0.0;

    if(r < 5)
    {
        c = ce[r]; ve = t2 * (c[0] + t2 * (c[1] + t2 * (c[2] + t2 * c[3])));
    
        if(r)
        {
            c = co[r - 1]; vo = c[0] + t2 * (c[1] + t2 * (c[2] + t2 * c[3])); 
            t2 = (float)q;
            
            for(i = r; i; i >>= 1)
            {
                if(i & 1)
                    vo *= t2;
                t2 *= t2;
            }
        }

        ve += (ty == type_b ? -vo : vo);
    }
    else if(q != 0.0)
    {
        ve = 0.5f * (float)q / (r2 - 1.0f); t2 = ve * ve;
        vo = 1.0f + t2 * (1.25f * r2 + 1.75f + t2 * 
            ((4.5f * r2 + 29.0f) * r2 + 14.5f) / (r2 - 9.0f)) / (r2 - 4.0f);
        ve *= ((float)q * vo);
    }

    return r2 + ve;
}

/*
 Asymptotic expansion for large q values.This provides a good initial 
 estimate for the blanch root finding procedure r > floor(sqrt(4 * q)).
*/

static float asy_est(const int ty, const int r, const double q)
{	float	iw2, v, ss;

    v = (2.0f * r + (ty == type_b ? -1.0f : 1.0f)); 
    iw2 = 1.0f / (v * v); 
    v *= 0.125f / (float)sqrt(q);

    ss = (175045.0f + iw2 * (9702612.0f + iw2 * (107798166.0f 
                    + iw2 * (288161796.0f + iw2 * 130610637.0f)))) / 2048.0f;

    ss = v * ss + (9387.0f + iw2 * (388780.0f + iw2 * (2845898.0f 
                    + iw2 * (4021884.0f + iw2 * 506979.0f)))) / 256.0f;

    ss = v * ss + (527.0f + iw2 * (15617.0f 
                    + iw2 * (69001.0f + iw2 * 41607.0f))) / 32.0f;

    ss = v * ss + (63.0f + iw2 * (1260.0f 
                    + iw2 * (2943.0f + iw2 * 486.0f))) / 8.0f;

    ss = v * ss + (33.0f + iw2 * (410.0f + iw2 * 405.0f)) / 8.0f;

    ss = v * ss + (5.0f + iw2 * (34.0f + iw2 * 9.0f)) / 2.0f;

    ss = v * ss + (1.0f + iw2 * 3.0f) * 2.0f;

    ss = v * ss + (1.0f + iw2) * 4.0f;

    return -2.0f * (float)q * (1.0f - v * (8.0f - v * ss));
}

/*
 Produce characteristic value estmates for low order curves
 or low and high q values: q <= 0.5 * r * r and q >= r * r
*/

double mathieu_characteristic_approx(int ty, double q, int r)
{	double	ql, qh;
    qh = r * r; ql = 0.5 * qh;
    if(r < 20)
    {
        ql = (r < 5 ? r + 1 : ql - 9); 
        qh += 31 - 12 * (ty == type_b ? 2 : 1);
    }

    if(q <= ql)						/* use low q expansion		*/
        return loq_est(ty, r, q);
    else if(q >= qh)				/* use asymptotic expansion	*/
        return asy_est(ty, r, q);
    else if(r < 20)					/* use rational expansions	*/
        return r * r + 
            ratval((float)q, coef[r ? 2 * r + (ty == type_b ? 0 : -1) : 0]);
    else							/* return 'no estimate'		*/
        return mathieu_no_estimate;
}
