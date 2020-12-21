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
    4. Maybe use ranges

*/

  

namespace RungeKutta
{
  //TODO: think about if you want type info
  typedef double constant_T;
  // intermediate step variables  
  static constant_T kmult [15] = {1.0/4.0,
                                  3.0/32.0,9.0/32.0,
                                  1932.0/2197.0,-7200.0/2197.0,7296.0/2197.0,
                                  439.0/216.0,-8.0,3680.0/513.0,-845.0/4104.0,
                                  -8.0/27.0,2.0,-3544.0/2565.0,1859.0/4104.0,-11.0/40.0};
  // y updates
  static constant_T yupdate [9] = {25.0/216.0,1408.0/2565.0,2197.0/4104.0,-1.0/5.0
                                  ,16.0/135.0,6656.0/12825.0,28561.0/56430.0,-9.0/50.0,2.0/55.0};

  // variables used in time step
  static constant_T ktstep [5] = {1.0/4.0,3.0/8.0,12.0/13.0,1.0,0.5};
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
            YVALS[i] = yf[i] + kmult[0]*k1[i];
          }
          Fvec(YVALS, Tstep + ktstep[0]*h,h,k2);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + kmult[1]*k1[i] + kmult[2]*k2[i];
          }
          Fvec(YVALS,Tstep + ktstep[1]*h,h,k3);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + kmult[3]*k1[i] + kmult[4]*k2[i] + kmult[5]*k3[i];
          }
          Fvec(YVALS,Tstep + ktstep[2]*h,h,k4);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + kmult[6]*k1[i] + kmult[7]*k2[i] + kmult[8]*k3[i] + kmult[9]*k4[i];
          }
          Fvec(YVALS,Tstep + ktstep[3]*h,h,k5);
          for(index i = 0;i < fDimSize;++i)
          {
            YVALS[i] = yf[i] + kmult[10]*k1[i] + kmult[11]*k2[i] + kmult[12]*k3[i] + kmult[13]*k4[i] + kmult[14]*k5[i];
          }
          Fvec(YVALS,Tstep + ktstep[4]*h,h,k6);
          for(index i = 0;i < fDimSize;++i)
          {
            yf[i] = yf[i] + (yupdate[0])*k1[i] + (yupdate[1])*k3[i] + (yupdate[2])*k4[i] + (yupdate[3])*k5[i];
            erres[i]=  (yupdate[4]-yupdate[0])*k1[i] +
  	        (yupdate[5]-yupdate[1])*k3[i] +
  	        (yupdate[6]-yupdate[2])*k4[i] +
  	        (yupdate[7]-yupdate[3])*k5[i] +
  	        (yupdate[8])*k6[i];
          }

          errestimate = err_norm(erres);
          errortol = err;
          if(errestimate/errortol > 1)
      	  {
	          tstep -= h;
	          h *=0.9*pow(errortol/errestimate,0.2);
            for(index i = 0;i < fDimSize;++i)
            {
  	          yf[i] = yf[i] - (yupdate[0])*k1[i] - (yupdate[1])*k3[i] - (yupdate[2])*k4[i] - (yupdate[3])*k5[i];
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
