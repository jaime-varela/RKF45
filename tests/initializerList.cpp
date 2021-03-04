#include "gtest/gtest.h"
#include "testutils/derivativeUpdateFunctions.hpp"

//system under test
#include "rkf45.h"

TEST(BasicNewtonTest, InitalizerList) {

    double AU = 1.496e11; // astronomical unit
    double V = 30000.0; // initial velocity

    double T = 7.2e6;
    double hi = 7.2e3/3;
    // solve the RK45 system with
    // t_0 = 0, t_f = T, y0, err = 1e-05, h_i = 7.2e/3
    std::vector<double> result(6);  

    RungeKutta::RK45 Newton(3*2,newtonForce);
    result = Newton.driver(0.0,T,{0.0,AU,0.0,-V,0.0,0.0},1e-5,hi);

    ASSERT_NEAR(-149956413569.63306/result[0], 1,1e-10);
    ASSERT_NEAR(21046790566.038616/result[1],1,1e-10);
    ASSERT_NEAR(0, result[2],1e-10);
    ASSERT_NEAR(-4531.5121133203666/result[3],1,1e-10);
}