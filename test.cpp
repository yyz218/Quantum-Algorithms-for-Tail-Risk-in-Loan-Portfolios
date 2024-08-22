#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <chrono>
#include <bitset>
#include <boost/math/distributions/normal.hpp>

#include "mc.h"
#include "analytic.h"

int main() {
    // set the profits and losses vector
    int n = 10000000;
    float alpha = 0.001;
    std::vector<float> defaults = {0,1,2};
    std::ofstream outFile("loss.csv");
    simulateBanks_MC("input_file.txt", n, outFile, alpha);
    //simulateBanks_AND("input_file.txt", n, outFile, alpha, defaults);
    //simulateBanks_OR("input_file.txt", n, outFile, alpha, defaults);
    std::ofstream outFile1("cdf.csv");
    analyticBanks_R("input_file.txt", outFile1, 1 - alpha);
    //analyticBanks_AND("input_file.txt", outFile1, 1 - alpha, defaults);
    //analyticBanks_OR("input_file.txt", outFile1, 1 - alpha, defaults);
    
    return 0;
}