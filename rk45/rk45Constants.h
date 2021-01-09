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

}