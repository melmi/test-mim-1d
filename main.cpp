#include <sstream>
#include <iostream>
#include <algorithm>

#include "domain.h"
#include "initializers.h"

int main(int argc, char const *argv[])
{
    domain *d = new domain(1000);
    d->dx = 5;
    d->dt = .6;
    initializers::init_dambreak(d, 4, 2, 0, d->n / 2);

    double tend = 150;
    int step_write = int(1.0 / d->dt);
    int nsteps = int(tend / d->dt);

    double sigma_xi_start, sigma_energy_start;
    d->get_state(sigma_xi_start, sigma_energy_start);

    d->write("./out/flow0.csv");
    for (int i = 0; i < nsteps; ++i)
    {
        d->step_calc();

        if ((i + 1) % step_write == 0)
        {
            std::stringstream ss;
            ss << "./out/flow" << (i + 1) / step_write << ".csv";
            std::cout << "writing time step: " << (i + 1) / step_write << std::endl;
            d->write(ss.str());
        }
    }

    double sigma_xi_end, sigma_energy_end;
    d->get_state(sigma_xi_end, sigma_energy_end);

    std::cout << "start: ( int(xi,dx):" << sigma_xi_start << "    int(energy,dx): " << sigma_energy_start << " )" << std::endl;
    std::cout << "end: ( int(xi,dx):" << sigma_xi_end << "    int(energy,dx): " << sigma_energy_end << " )" << std::endl;

    delete d;
    return 0;
}