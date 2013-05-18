#include <fstream>

#ifndef DOMAIN_H
#define DOMAIN_H

class domain
{
public:
    const double gravity = 9.81;

    enum bctype {constxi, constp};

    double *xi, *d, *p, *h, *u;
    double dx, dt;
    int n;
    bctype lbc, rbc;
    double lbc_val, rbc_val;

    domain(int n)
    {
        this->n = n;

        xi = new double[n];
        d = new double[n];
        p = new double[n];
        h = new double[n];
        u = new double[n];
    }

    double *get_grad(double *phi)
    {
        double *grad = new double[n];
        for (int i = 1; i < n - 1; ++i)
        {
            grad[i] = (phi[i + 1] - phi[i - 1]) / (2.*dx);
        }
        grad[0] = (phi[1] - phi[0]) / dx;
        grad[n - 1] = (phi[n - 1] - phi[n - 2]) / dx;

        return grad;
    }

    void advect(double *phi)
    {
        double *mass = new double[n];
        for (int i = 0; i < n; ++i)
        {
            mass[i] = phi[i] * dx;
        }
        double *grad = get_grad(phi);
        for (int face = 0; face < n - 1; ++face)
        {
            int i = face, ii = face + 1;
            double uf = 0.5 * (u[i] + u[ii]);
            double flux;
            if (uf > 0)
                flux = (phi[i] + (dx / 2. - uf * dt / 2.) * grad[i]) * uf * dt;
            else
                flux = (phi[ii] + (-dx / 2. - uf * dt / 2.) * grad[ii]) * uf * dt;
            mass[i] -= flux;
            mass[ii] += flux;
        }
        for (int i = 0; i < n; ++i)
        {
            phi[i] = mass[i] / dx;
        }
        delete[] grad;
        delete[] mass;
    }

    void advect_p()
    {
        double *mass = new double[n];
        for (int i = 0; i < n; ++i)
            mass[i] = p[i] * dx;

        double *grad = get_grad(p);
        double *grad_hdxidx = get_pressure_grad_narrow(xi);
        double gdt = gravity * dt;
        // for (int i = 0; i < n; ++i)
        //     grad[i] -= gdt * grad_hdxidx[i];

        for (int face = 0; face < n - 1; ++face)
        {
            int i = face, ii = face + 1;
            double uf = 0.5 * (u[i] + u[ii]) - gdt * (u[ii] - u[i]) / dx; //* grad_hdxidx[i] / h[i];
            double flux;
            if (uf > 0)
                flux = (p[i] + (dx / 2. - uf * dt / 2.) * grad[i]) * uf * dt;
            else
                flux = (p[ii] + (-dx / 2. - uf * dt / 2.) * grad[ii]) * uf * dt;
            mass[i] -= flux;
            mass[ii] += flux;
        }
        for (int i = 0; i < n; ++i)
        {
            p[i] = mass[i] / dx;
        }
        delete[] grad;
        delete[] mass;
        delete[] grad_hdxidx;
    }

    void calc_dependent_vars()
    {
        for (int i = 0; i < n; ++i)
        {
            h[i] = d[i] + xi[i];
            u[i] = p[i] / h[i];
        }
    }

    void write(std::string fname)
    {
        std::ofstream file(fname.c_str());
        file << "x, xi, pi" << std::endl;
        for (int i = 0; i < n; ++i)
            file << i *dx << ", " << xi[i] << ", " << p[i] << std::endl;
        file.close();
    }

    void apply_p_bc()
    {
        if (lbc == constp) p[0] = lbc_val;
        if (rbc == constp) p[n - 1] = rbc_val;
    }

    double *get_pressure_grad(double *xi_prev)
    {
        double *grad_xi = get_grad(xi_prev);
        for (int i = 0; i < n; ++i)
            grad_xi[i] *= h[i];
        double *result = get_grad(grad_xi);
        delete[] grad_xi;
        return result;
    }

    double *get_pressure_grad_narrow(double *xi_prev)
    {
        double *pressure_term = new double[n];
        for (int i = 1; i < n - 1; ++i)
        {
            int ie = i, iw = i - 1;
            double e = 0.5 * (h[ie] + h[ie + 1]) * (xi_prev[ie + 1] - xi_prev[ie]) / dx;
            double w = 0.5 * (h[iw] + h[iw + 1]) * (xi_prev[iw + 1] - xi_prev[iw]) / dx;
            // double ue = 0.5 * (u[ie] + u[ie + 1]);
            // double uw = 0.5 * (u[iw] + u[iw + 1]);
            // double e = (ue > 0 ? h[ie] : h[ie + 1]) * (xi_prev[ie + 1] - xi_prev[ie]) / dx;
            // double w = (uw > 0 ? h[iw] : h[iw + 1]) * (xi_prev[iw + 1] - xi_prev[iw]) / dx;
            pressure_term[i] = (e - w) / dx;
        }
        pressure_term[0] = (0.5 * (h[0] + h[1]) * (xi_prev[1] - xi_prev[0]) / dx)  / dx;
        pressure_term[n - 1] = (-0.5 * (h[n - 2] + h[n - 1]) * (xi_prev[n - 1] - xi_prev[n - 2]) / dx) / dx;
        return pressure_term;
    }

    double *calc_xi(double *xi_prev, double *div_p)
    {
        double gdt2 = gravity * dt * dt;
        double *new_xi = new double[n];
        double *pressure_term = get_pressure_grad(xi_prev);
        // double *narrow_pressure_term = get_narrow_pressure_term(xi_prev);//mim

        for (int i = 0; i < n; ++i)
            // new_xi[i] = xi[i] - dt * div_p[i] + pressure_term[i] - narrow_pressure_term[i]; //mim
            new_xi[i] = xi[i] - dt * div_p[i] + gdt2 * pressure_term[i];

        delete[] pressure_term;
        // delete[] narrow_pressure_term;

        if (lbc == constxi) new_xi[0] = lbc_val;
        if (rbc == constxi) new_xi[n - 1] = rbc_val;

        return new_xi;
    }

    void update_p()
    {
        double gdt = gravity * dt;
        double *grad_xi = get_grad(xi);
        for (int i = 0; i < n; ++i)
            p[i] -= gdt * h[i] * grad_xi[i];

        delete[] grad_xi;
    }

    void calc_xi_exp()
    {
        double *div_p = get_grad(p);
        double *new_xi = calc_xi(xi, div_p);
        delete[] div_p;
        delete[] xi;
        xi = new_xi;
    }

    void calc_xi_imp(double eps)
    {
        double *div_p = get_grad(p);
        double *old_xi = new double[n];
        std::copy(xi, xi + n, old_xi);

        for (double err = 1; err > eps;)
        {
            double *new_xi = calc_xi(xi, div_p);

            for (int i = 0; i < n; ++i)
                err += (new_xi[i] - old_xi[i]) * (new_xi[i] - old_xi[i]);
            err /= n;

            delete[] old_xi;
            old_xi = new_xi;
        }

        delete[] div_p;
        delete[] xi;
        xi = old_xi;
    }

    void step_calc()
    {
        // advect(p);
        advect_p();
        apply_p_bc();
        calc_xi_exp();
        // calc_xi_imp(1.e-6);
        // eliminate_xi_nullspace();
        update_p();
        apply_p_bc();
        calc_dependent_vars();
    }

    static double det_2x2(double a11, double a12, double a21, double a22)
    {
        return a11 * a22 - a21 * a12;
    }

    static void solve_2x2(double a11, double a12, double a21, double a22, double b1, double b2, double &x1, double &x2)
    {
        double d = det_2x2(a11, a12, a21, a22);
        x1 = det_2x2(b1, a12, b2, a22) / d;
        x2 = det_2x2(a11, b1, a21, b2) / d;
    }

    void get_elim_coeffs(double &a0, double &a1)
    {
        double a11, a12, a21, a22, b1, b2;
        a11 = a22 = n;
        a12 = a21 = (n + 1) % 2;
        double p1i = -1;
        b1 = b2 = 0;
        for (int i = 0; i < n; ++i)
        {
            p1i *= -1;
            b1 += xi[i];
            b2 += p1i * xi[i];
            // std::cout<<p1i<<" "<<xi[i]<<" "<<b2<<std::endl;
            // std::getchar();
        }
        solve_2x2(a11, a12, a21, a22, b1, b2, a0, a1);

        a0 = -a0;
        a1 = -a1;
    }

    void eliminate_xi_nullspace()
    {
        double a0, a1;
        get_elim_coeffs(a0, a1);

        double p1i = -1;
        for (int i = 0; i < n; ++i)
        {
            p1i *= -1;
            xi[i] += a0 + p1i * a1;
        }
    }

    void get_state(double &sigma_mass, double &sigma_energy)
    {
        //double *grad_u = get_grad(u);
        sigma_mass = sigma_energy = 0;
        for (int i = 0; i < n; ++i)
        {
            sigma_mass += h[i] * dx;
            sigma_energy += (xi[i] + u[i] * u[i] / 2. / gravity) * h[i] * dx;
            // sigma_energy += (xi[i] + (u[i] * u[i] + grad_u[i] * grad_u[i] * dx * dx / 12. ) / 2. / gravity) * dx;
        }
        // delete[] *grad_u;
    }
};

#endif
