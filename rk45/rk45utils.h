#pragma once

#include "rk45Concepts.h"


namespace RungeKutta
{
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

  //  method returns Nth power of A 
  template<NumT number = double,Index index = int>
  number integerPower(number A, index N) 
  {
    number result = A;
    for(index v = 1;v < N;v++)
    {
      result *= A;
    }
    return result;
  }

//  method returns 1/Nth power of A 
  template<NumT number = double,Index index = int>
  number integerRootApprox(number A, index N) 
  { 
    // intially based on taylor series approximation
    number xPre = (A > 1)? A/2.0 : 1- ((1-A)/(1.0*N)) - (((1-A)*(1-A)*(1.0*N -1))/(2.0*N*N));
  
    //  smaller eps, denotes more accuracy 
    // tunable but don't need full precision as this only does approximate output
    number eps = 1e-4;
  
    // initializing difference between two 
    number delX;
  
    //  xK denotes current value of x 
    number xK; 
  
    //  loop untill we reach desired accuracy 
    do 
    { 
        //  calculating current value from previous 
        // value by newton's method 
        xK = ((N - 1.0) * xPre + A/integerPower(xPre, N-1)) / (1.0*N); 
        delX = abs(xK - xPre); 
        xPre = xK; 
    } while (delX > eps);   
  
    return xK; 
  }

}