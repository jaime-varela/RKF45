#include <stdio.h>
#include <iostream>
#include "rkf45.h"

#include <chrono>

int main(){
    std::vector<double> x = {0.12,12.23,12.3,1.3};
    auto val = RungeKutta::err_norm(x);
    std::cout << val << std::endl;

    return 0;
}
