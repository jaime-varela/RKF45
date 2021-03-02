#include "derivativeUpdateFunctions.hpp"

double newtonForce(int i, double t, const std::vector<double> & xv){
  if(i < 3)
  {
    // xv[3] = dr/dt,...
    return xv[i+3];
  }
  double r;
  r = sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
  return -(alpha*xv[i-3]) / (r*r*r);
}