#include <stdio.h>
#include <iostream>
#include <chrono>

#include "rkf45.h"

// Newtonion force.
// all quantities in SI units.
const double alpha = 1.32754125e20;
double Force3D(int i, double t, const std::vector<double> & xv){
  if(i < 3)
  {
    // xv[3] = dr/dt,...
    return xv[i+3];
  }
  double r;
  r = sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
  return -(alpha*xv[i-3]) / (r*r*r);
}

int main(){

  double AU = 1.496e11; // astronomical unit
  double V = 30000.0; // initial velocity
  // initial conditions
  std::vector<double> y0 = {0.0,AU,0.0,-V,0.0,0.0}; 

  double T = 7.2e6;
  double hi = 7.2e3/3;
  // solve the RK45 system with
  // t_0 = 0, t_f = T, y0, err = 1e-05, h_i = 7.2e/3
  std::vector<double> result(6);  

  RungeKutta::RK45 Newton(3*2,Force3D); 

  result = Newton.driver(0.0,T,y0,1e-5,hi);

  // Print solution y_0(tf), y_1(tf) ...
  for(auto&& value : result)
    std::cout << value << ",";
  std::cout << std::endl;



  auto Nruns = 5000;
  auto totalDuration = 0;
  auto t1 = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < Nruns; i++)
  {
    Newton.driver(0.0,T,y0,1e-5,hi);
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  totalDuration += duration;    
  std::cout << totalDuration/Nruns << std::endl;

  return 0;
}
