#include <R.h>
void rpbounds (double *h1l,     double *h1u,
               double *h2l,     double *h2u,
               double *q,       int    *num,
               double *fl,      double *fu,
               double *rplower, double *rpupper)
{
    double hl[*num];
    double hu[*num];

    double p;
    double facl;
    double qfacl;
    double pfacl;

    /* Convolution of h1l and h2l and of h1u and h2u, respectively */
    for (int k = 0; k <= *num; k++)
    {
        /* Initialize */
        hl[k] = 0.0;
        hu[k] = 0.0;

        /* Convolve */
        for (int n = 0; n <= k; n++)
        {
            hl[k] += h1l[n] * h2l[k-n];
            hu[k] += h1u[n] * h2u[k-n];
        }
    }

    /* Define some variables for quantities that will be re-used. */
    facl  = 1.0 / (1.0 - *q * hl[0]);
    p     = 1.0 - *q;
    qfacl = *q * facl;
    pfacl = p * facl;

    /* Fill in the initial values. */
    fl[0] = h1l[0] * p * facl;
    fu[0] = 0.0;

    rplower[0] = 1.0;
    rpupper[0] = 1.0;

    /* Successively fill each component of fl, fu, rplower and rpupper. */
    for (int i = 1; i <= *num; i++)
    {
        fl[i] = 0.0;
        fu[i] = 0.0;

        /* The actual recursion happens here. */
        for (int j = 1; j <= i; j++)
        {
            fl[i] += hl[j] * fl[i-j];
            fu[i] += hu[j] * fu[i-j];
        }

        /* Multiply resulting values with constant factors ... */
        fl[i] *= qfacl;
        fu[i] *= *q;

        /* ... and insert the additive terms. */
        fl[i] += h1l[i] * pfacl;
        fu[i] += h1u[i] * p;

        /* Construct the vectors of bounds, which are actually 1 minus the
         * vector containing the cumulative sums of fl and fu, respectively. */
        rplower[i] = rplower[i-1] - fl[i-1];
        rpupper[i] = rpupper[i-1] - fu[i];
    }
}
