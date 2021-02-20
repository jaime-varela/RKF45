// Use the geodesic solver code to create a gsl implementation so you can
// compare the speeds.  The template needs to get close to ,ideally better than,
// exactly to the gsl run times

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

const double alpha = 1.32754125e20;
int odefunc (double x, const double y[], double f[], void *params)
{
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    double r = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);

    f[3] = -(alpha*y[0]) / (r*r*r);
    f[4] = -(alpha*y[1]) / (r*r*r);
    f[5] = -(alpha*y[2]) / (r*r*r);
    return GSL_SUCCESS;
}


int * jac;

int main ()
{
    double AU = 1.496e11; // astronomical unit
    double V = 30000.0; // initial velocity
    // initial conditions
    double y0[] = {0.0,AU,0.0,-V,0.0,0.0}; 
    double T0 = 0.0;
    double Tf = 7.2e6;
    double hi = 7.2e3/3;
    int dim = 6;
    gsl_odeiv2_system sys = {odefunc, NULL, dim, NULL};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, hi, 1e-5, 0.0);
    auto Nruns = 5000;
    auto totalDuration = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    double T;
    for (int i = 1; i <= Nruns; i++)
    {
        T0 = 0.0;
        T = Tf;
        y0[0]=0.0;
        y0[1] = AU;
        y0[2] = 0.0;
        y0[3] = -V;
        y0[4] = 0.0;
        y0[5] = 0.0;
        gsl_odeiv2_driver_apply (d, &T0, T, y0);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    totalDuration += duration;
    std::cout << "GSL runs" << std::endl;
    std::cout << totalDuration/Nruns << std::endl;

    gsl_odeiv2_driver_free (d);
    return 0;
}