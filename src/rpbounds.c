#include <R.h>
void rpbounds (
        double *h1l, double *h1u,
        double *h2l, double *h2u,
        double *q, int *num,
        double *fl, double *fu,
        double *rplower, double *rpupper
        )
{
    double hl[*num];
    double hu[*num];

    // Convolution of h1l and h2l and of h1u and h2u, respectively
    for (int k = 0; k <= *num; k++)
    {
        hl[k] = 0;
        hu[k] = 0;
        for (int n = 0; n <= k; n++)
        {
            hl[k] += h1l[n] * h2l[k-n];
            hu[k] += h1u[n] * h2u[k-n];
        }
    }

    double facl = 1 / (1 - *q * hl[0]);
    double p = 1 - *q;

    fl[0] = h1l[0] * p * facl;
    fu[0] = 0;

    rplower[0] = 1;
    rpupper[0] = 1;

    for (int i = 1; i <= *num; i++)
    {
        fl[i] = 0; fu[i] = 0;

        for (int j = 1; j <= i; j++)
        {
            fl[i] += hl[j] * fl[i-j];
            fu[i] += hu[j] * fu[i-j];
        }

        fl[i] *= *q * facl;
        fu[i] *= *q;

        fl[i] += h1l[i] * p * facl;
        fu[i] += h1u[i] * p;

        rplower[i] = rplower[i-1] - fl[i-1];
        rpupper[i] = rpupper[i-1] - fu[i];
    }
}
