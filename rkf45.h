#ifndef RUNGE_KUTTA_45_H
#define RUNGE_KUTTA_45_H

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <functional>
#include <valarray>  
//CONSTANTS; TODO: make tables out of these humbers
double k21 = 1.0/4.0, k2t = 1.0/4.0;
double k31 = 3.0/32.0,k32 = 9.0/32.0, k3t = 3.0/8.0;
double k41 = 1932.0/2197.0,k42 = -7200.0/2197.0,k43 = 7296.0/2197.0, k4t = 12.0/13.0;
double k51 = 439.0/216.0,k52 = -8.0, k53 = 3680.0/513.0, k54 = -845.0/4104.0, k5t = 1.0;
double k61 = -8.0/27.0,k62 = 2.0,k63 = -3544.0/2565.0,k64 = 1859.0/4104.0,k65 = -11.0/40.0,k6t = 0.5;


double y41 = 25.0/216.0;
double y43 = 1408.0/2565.0;
double y44 = 2197.0/4104.0;
double y45 = -1.0/5.0;


double y51 = 16.0/135.0;
double y53 = 6656.0/12825.0;
double y54 = 28561.0/56430.0;
double y55 = -9.0/50.0;
double y56 = 2.0/55.0;
  
//error calculation function
double err_norm(const std::valarray<double> & a)		
{
  double result =0;
  for(int i = 0; i < a.size();i++)
    {
      result += (a[i])*(a[i]);
    }
  return sqrt(result);
}

/*
----------- Class RK45 ---------------------

Adaptive-step runge-kutta method.

Public Methods: 

RK45::RK45(double(*f)(int,double,std::valarray<double> &))

constructor with a function object f

RK45::driver(double t0, double tf, std::valarray<double> & y0,double err, double h_i)

The driver method returns a vector with values  y_i(t_f)

err is the absolute error.
h_i is the initial step.
 */

class RK45
{
public:
  RK45(double(*f)(int,double,std::valarray<double> &)) {deriv = f;}
  
  std::valarray<double> driver(double t0,
			     double tf,
			     const std::valarray<double> & y0,
			       double err,double hi);
private:
  std::function<double(int,double,std::valarray<double> &)> deriv;
  void Fvec(std::valarray<double> & yvec, double T,
			     double & hv, std::valarray<double> & k);
  
};

void RK45::Fvec(std::valarray<double> & yvec, double T,
				 double & hv, std::valarray<double> & k)
{
  for(int i = 0; i< yvec.size();i++)
    k[i] = hv*(deriv(i,T,yvec) );
}


std::valarray<double> RK45::driver(double t0, double tf,
				 const std::valarray<double> & y0,
				   double err, double hi)
{
  int DIM = y0.size();
  double errestimate,errortol;
  double tstep = t0;
  double Tstep;
  std::valarray<double> yf(DIM),erres(DIM),YVALS(DIM);
  std::valarray<double> k1(DIM),k2(DIM),k3(DIM),k4(DIM),k5(DIM),k6(DIM);
  yf = y0;
  double h = hi;
  do
    {
      tstep += h;
      Tstep = tstep-h;//function eval time is h-step away from tstep.
      Fvec(yf,Tstep,h,k1);
      YVALS = yf + k21*k1;
      Fvec(YVALS, Tstep + k2t*h,h,k2);
      YVALS = yf + k31*k1 + k32*k2;
      Fvec(YVALS,Tstep + k3t*h,h,k3);
      YVALS = yf + k41*k1 + k42*k2 + k43*k3;
      Fvec(YVALS,Tstep + k4t*h,h,k4);
      YVALS = yf + k51*k1 + k52*k2 + k53*k3 + k54*k4;
      Fvec(YVALS,Tstep + k5t*h,h,k5);
      YVALS = yf + k61*k1 + k62*k2 + k63*k3 + k64*k4 + k65*k5;
      Fvec(YVALS,Tstep + k6t*h,h,k6);
      
      yf = yf + (y41)*k1 + (y43)*k3 + (y44)*k4 + (y45)*k5;
      erres=  (y51-y41)*k1 +
	(y53-y43)*k3 +
	(y54-y44)*k4 +
	(y55-y45)*k5 +
	(y56)*k6;
      
      errestimate = err_norm(erres);
      //      errortol = err*err_norm(yf); // relative error
      errortol = err;
      if(errestimate/errortol > 1)
	{
	  tstep -= h;
	  h *=0.9*pow(errortol/errestimate,0.2);
	  yf = yf - (y41)*k1 - (y43)*k3 - (y44)*k4 - (y45)*k5;
	}
      else
	  h *=0.9*pow(errortol/errestimate,0.25);
      //final step
      if(tstep + h > tf)
	  h = tf - tstep;
    }while(tstep < tf);  

  return yf;  
}

#endif
