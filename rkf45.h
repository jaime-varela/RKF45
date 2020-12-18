#ifndef RUNGE_KUTTA_45_H
#define RUNGE_KUTTA_45_H

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <concepts>
#include <iterator>
/*
  TODO: 

    1. Fused multiply adds
    2. constant handling
    3. think about implementation
    4. maybe use algorithms such as std::transform

    https://en.cppreference.com/w/cpp/iterator#C.2B.2B20_iterator_concepts

*/

  

namespace RungeKutta
{
  //CONSTANTS; TODO: make tables out of these humbers
  //TODO: think about if you want type info
  typedef double constant_T;
  static constant_T k21 = 1.0/4.0, k2t = 1.0/4.0;
  static constant_T k31 = 3.0/32.0,k32 = 9.0/32.0, k3t = 3.0/8.0;
  static constant_T k41 = 1932.0/2197.0,k42 = -7200.0/2197.0,k43 = 7296.0/2197.0, k4t = 12.0/13.0;
  static constant_T k51 = 439.0/216.0,k52 = -8.0, k53 = 3680.0/513.0, k54 = -845.0/4104.0, k5t = 1.0;
  static constant_T k61 = -8.0/27.0,k62 = 2.0,k63 = -3544.0/2565.0,k64 = 1859.0/4104.0,k65 = -11.0/40.0,k6t = 0.5;
  static constant_T y41 = 25.0/216.0;
  static constant_T y43 = 1408.0/2565.0;
  static constant_T y44 = 2197.0/4104.0;
  static constant_T y45 = -1.0/5.0;
  static constant_T y51 = 16.0/135.0;
  static constant_T y53 = 6656.0/12825.0;
  static constant_T y54 = 28561.0/56430.0;
  static constant_T y55 = -9.0/50.0;
  static constant_T y56 = 2.0/55.0;


  // ---------------- Begin Concepts ---------------------------
  template<class T>
  concept NumT = std::is_arithmetic<T>::value; 

  template<class T>
  concept Index = std::is_integral<T>::value;

  template<class container,class number>
  concept CoordinateContainer = std::is_arithmetic<number>::value &&
  requires(container A,number b) {
    {A[0] < A[0]} -> std::same_as<bool>;
    {A[0] + b} -> std::same_as<number>;
  };

  // ------------------- end concepts ------------------
  
  //error calculation function
  template<NumT number = double,Index index = int,CoordinateContainer<number> coords>
  number err_norm(const coords& a)		
  {
    number result =0;
    for(index i = 0; i < a.size();i++)
    {
      result += (a[i])*(a[i]);
    }
    return sqrt(result);
  }

  template<NumT number = double,Index index= int,CoordinateContainer<number> coords = std::vector<double>>
  class RK45
  {
    public:
      RK45(index numEqns,number(*F_i)(index,number,const coords &))
      {
        fDeriv = F_i;
        fDimSize = numEqns;
      }

      index getSize() const
      {
        return fDimSize;
      }

      coords driver(number t0,number tf,const coords & y0,number err,number hi)
      {
        index DIM = fDimSize;
        number errestimate,errortol;
        number tstep = t0;
        number Tstep;
        coords yf(DIM),erres(DIM),YVALS(DIM);
        coords k1(DIM),k2(DIM),k3(DIM),k4(DIM),k5(DIM),k6(DIM);
        yf = y0;
        number h = hi;
        do
        {
          tstep += h;
          Tstep = tstep-h;//function eval time is h-step away from tstep.
          Fvec(yf,Tstep,h,k1);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + k21*k1[i];
          }
          Fvec(YVALS, Tstep + k2t*h,h,k2);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + k31*k1[i] + k32*k2[i];
          }
          Fvec(YVALS,Tstep + k3t*h,h,k3);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + k41*k1[i] + k42*k2[i] + k43*k3[i];
          }
          Fvec(YVALS,Tstep + k4t*h,h,k4);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + k51*k1[i] + k52*k2[i] + k53*k3[i] + k54*k4[i];
          }
          Fvec(YVALS,Tstep + k5t*h,h,k5);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + k61*k1[i] + k62*k2[i] + k63*k3[i] + k64*k4[i] + k65*k5[i];
          }
          Fvec(YVALS,Tstep + k6t*h,h,k6);
          for(index i = 0;i < fDimSize;++i)
          {
            yf[i] = yf[i] + (y41)*k1[i] + (y43)*k3[i] + (y44)*k4[i] + (y45)*k5[i];
            erres[i]=  (y51-y41)*k1[i] +
  	        (y53-y43)*k3[i] +
  	        (y54-y44)*k4[i] +
  	        (y55-y45)*k5[i] +
  	        (y56)*k6[i];
          }

          errestimate = err_norm(erres);
          errortol = err;
          if(errestimate/errortol > 1)
      	  {
	          tstep -= h;
	          h *=0.9*pow(errortol/errestimate,0.2);
            for(index i = 0;i < fDimSize;++i)
            {
  	          yf[i] = yf[i] - (y41)*k1[i] - (y43)*k3[i] - (y44)*k4[i] - (y45)*k5[i];
            }
	        }
          else
          {
	          h *=0.9*pow(errortol/errestimate,0.25);
          }
          //final step
          if(tstep + h > tf)
          {
      	    h = tf - tstep;
          }
        } while(tstep < tf);  //TODO add a max_iter to avoid infinite loop

      return yf;  
    }

    private:
      std::function<number(index,number,const coords&)> fDeriv;
      index fDimSize;
      /*
        computes k[i] = hv*F_i(T,x);
      */
      void Fvec(const coords& yvec, number T,number hv, coords& k) const
      {
        for(index i = 0; i< fDimSize;i++)
        {
          k[i] = hv*(fDeriv(i,T,yvec) );
        }
      }
  };

} // namespace RK
#endif
