#include <stdio.h>
#include <iostream>
#include "rkf45.h"

// Newtonion force.
// all quantities in SI units.
const double alpha = 1.32754125e20;
double Force3D(int i, double t, std::valarray<double> & xv){
  double r;
  r = sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]); 
  switch(i) {
  case 0:
    return xv[3];
  case 1:
    return xv[4];
  case 2:
    return xv[5];
  case 3:
    return - (alpha*xv[0]) / (r*r*r);
  case 4:
    return -(alpha*xv[1]) / (r*r*r) ;
  case 5:
    return -(alpha*xv[2]) / (r*r*r) ;
  }
}

int main(){

  double AU = 1.496e11; // astronomical unit
  double V = 30000.0; // initial velocity
  // initial conditions
  std::valarray<double> y0 = {0.0,AU,0.0,-V,0.0,0.0}; 

  double T = 7.2e6;
  double hi = 7.2e3/3;
  // initialize the RK45 class
  RK45 Newton(Force3D);

  // solve the RK45 system with
  // t_0 = 0, t_f = T, y0, err = 1e05, h_i = 7.2e/3
  std::valarray<double> result;  
  result = Newton.driver(0.0,T,y0,1e-5,hi);

  // Print solution y_0(tf), y_1(tf) ...
  for(int i = 0; i < result.size(); i++)
    std::cout << result[i] << ",";  
  std::cout << std::endl;
  return 0;
}
