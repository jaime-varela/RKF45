#ifndef RUNGE_KUTTA_45_H
#define RUNGE_KUTTA_45_H

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <ranges>

#include "rk45Concepts.h"
#include "rk45Constants.h"
#include "rk45utils.h"
/*
  TODO: 

    1. Revamp everything to use a Butcher Tableau inherited class. This will allow one to use more methods than just the 45 method
    2. Put reference formulas in the documentation and use "Solving Ordinary Differential Equations I: Nonstiff Problems" notation

*/
  

namespace RungeKutta
{

  using defaultFunctionPointer = double(*)(int,double,const std::vector<double>&);

  template<NumT number = double,
          Index index= int,
          CoordinateContainer<number> coords = std::vector<double>,
          DerivativeFunction<number,index,coords> derivativeFunc = defaultFunctionPointer>
  class RK45
  {
    public:
      RK45(index numEqns,derivativeFunc F_i) : fDeriv(F_i)
      {
        fDimSize = numEqns;
      }

      index getSize() const
      {
        return fDimSize;
      }
      /**
       * @brief Solves the ODE to obtain y(t_f) given the initial condition.
       *
       * @param t0 initial time.
       * @param tf desired time of solution.
       * @param y0 the initial condition
       * @param err error tolerance
       * @param hi initial step size
       * @return y(tf) for the initial condition.
       */
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
          Tstep = tstep;//function eval time is h-step away from tstep.
          tstep += h;
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
  	        (yupdateDiff1)*k3[i] +
  	        (yupdateDiff2)*k4[i] +
  	        (yupdateDiff3)*k5[i] +
  	        (yupdate[8])*k6[i];
          }

          errestimate = err_norm(erres);
          if(errestimate > err)
      	  {
	          tstep -= h;
	          h *=0.9* integerRootApprox(err/errestimate,5);
            for(index i = 0;i < fDimSize;++i)
            {
  	          yf[i] = yf[i] - (yupdate[0])*k1[i] - (yupdate[1])*k3[i] - (yupdate[2])*k4[i] - (yupdate[3])*k5[i];
            }
	        }
          else
          {
	          h *=0.95*sqrt(sqrt(err/errestimate));
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
      derivativeFunc fDeriv;

      index fDimSize;
      /*
        computes k[i] = hv*F_i(T,x);
      */
      void Fvec(const coords& yvec, number T,number hv, coords& k) const
      {
        for (auto i : std::views::iota(0, fDimSize)) 
        {
          k[i] = hv * fDeriv(i, T, yvec);
        }
      }
  };

} // namespace RK
#endif
