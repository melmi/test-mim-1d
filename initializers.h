#include "domain.h"
#include <algorithm>

#ifndef INITIALIZERS_H
#define INITIALIZERS_H

class initializers
{
public:
    static void init_zero(domain *d, double depth)
    {
        d->lbc_val = 0;
        d->rbc_val = 0;

        std::fill_n(d->u, d->n, 0);
        std::fill_n(d->ustar, d->n, 0);
        std::fill_n(d->p, d->n, 0);
        std::fill_n(d->d, d->n, depth);
        std::fill_n(d->h, d->n, 0);
        std::fill_n(d->xi, d->n, 0);

        d->calc_dependent_vars();
    }

    static void init_dambreak(domain *d, double depth, double xil, double xir, int n0)
    {
        init_zero(d, depth);
        d->lbc = d->rbc = domain::constp;
        std::fill_n(d->xi, n0, xil);
        std::fill_n(d->xi + n0, d->n - n0, xir);
        d->calc_dependent_vars();
    }
};

#endif
