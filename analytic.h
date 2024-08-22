#ifndef ANALYTIC_H
#define ANALYTIC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>

// compute P(Z)
float P_Z(float z, float mu, float sigma) {
    if (z <= 0) { return 0; }
    float part1 = 1 / (z * sigma * std::sqrt(2 * M_PI));
    float part2 = -std::pow(std::log(z) - mu, 2) / (2 * std::pow(sigma, 2));
    return part1 * std::exp(part2);
}

// Compute P(A|Z)
float P_AZ(float z, float Z_min, const std::vector<std::vector<float>>& defaults) {
    if (z <= 0) { return 0; }
    float P_AZ = 1;
    if (z < Z_min) {
        for (size_t i = 0; i < defaults.size(); i++) {
            float X = defaults[i][8] - defaults[i][5] * z;
            float z_score = (std::log(X / (1 - defaults[i][5])) - defaults[i][6]) / defaults[i][7];
            float P_default = 0.5 * erfc(-(z_score)/sqrt(2));
            P_AZ = P_AZ * P_default;
        }
    } else {
        P_AZ = 0;
    }
    return P_AZ;
}

// Compute P(A|Z) for conditional OR
double P_AZ_or(float z, float Z_max, const std::vector<std::vector<float>>& defaults) {
    if (z <= 0) { return 0; }
    double P_AZ = 1;
    if (z <= Z_max) {
        for (size_t i = 0; i < defaults.size(); i++) {
            float X = defaults[i][8] - defaults[i][5] * z;
            double P_default = 0;
            if (X > 0) {
                float z_score = (std::log(X / (1 - defaults[i][5])) - defaults[i][6]) / defaults[i][7];
                P_default = 0.5 * erfc(-(z_score)/sqrt(2));
            }
            P_AZ = P_AZ * (1 - P_default);
        }
    } else {
        P_AZ = 0;
    }

    return 1 - P_AZ;
}

// Compute E[X] and Var[X] for P(X<=h|Z)
std::vector<float> X_Z(float z, const std::vector<std::vector<float>>& banks) {
    float mean = 0;
    float variance = 0;
    for (size_t i = 0; i < banks.size(); i++) {
        float X = banks[i][8] - banks[i][5] * z;
        float P_default = 0;
        if (X > 0){
            float z_score = (std::log(X / (1 - banks[i][5])) - banks[i][6]) / banks[i][7];
            P_default = 0.5 * erfc(-(z_score)/sqrt(2));
        }
        float expected_profit = P_default * ( banks[i][4] - 1) * banks[i][2] + (1 - P_default) * banks[i][3] * banks[i][2];
        variance += P_default * std::pow((banks[i][4] - 1) * banks[i][2], 2) + (1 - P_default) * std::pow(banks[i][3] * banks[i][2], 2) - std::pow(expected_profit, 2);
        mean += expected_profit;
    }
    variance = sqrt(variance);
    return {mean, variance};
}

// Compute E[X] and Var[X] for P(X<=h|Z) for conditional OR
std::vector<float> X_Z_or(float z, const std::vector<std::vector<float>>& banks, float P_AZ, const std::vector<std::vector<float>>& defaults) {
    float mean = 0;
    float variance = 0;
    for (size_t i = 0; i < banks.size(); i++) {
        float X = banks[i][8] - banks[i][5] * z;
        float P_default = 0;
        if (X > 0){
            float z_score = (std::log(X / (1 - banks[i][5])) - banks[i][6]) / banks[i][7];
            P_default = 0.5 * erfc(-(z_score)/sqrt(2));
        }
        float expected_profit = P_default * ( banks[i][4] - 1) * banks[i][2] + (1 - P_default) * banks[i][3] * banks[i][2];
        variance += P_default * std::pow(( banks[i][4] - 1) * banks[i][2], 2) + (1 - P_default) * std::pow(banks[i][3] * banks[i][2], 2) - std::pow(expected_profit, 2);
        mean += expected_profit;
    }
    for (size_t i = 0; i < defaults.size(); i++) {
        float X = banks[i][8] - banks[i][5] * z;
        float P_default = 0;
        if (X > 0){
            float z_score = (std::log(X / (1 - banks[i][5])) - banks[i][6]) / banks[i][7];
            P_default = 0.5 * erfc(-(z_score)/sqrt(2));
            P_default = P_AZ == 0 ? P_default : P_default / P_AZ;
        }
        P_default = (P_default < 0) ? 0 : ((P_default > 1) ? 1 : P_default);
        float expected_profit = P_default * ( banks[i][4] - 1) * banks[i][2] + (1 - P_default) * banks[i][3] * banks[i][2];
        variance += P_default * std::pow(( banks[i][4] - 1) * banks[i][2], 2) + (1 - P_default) * std::pow(banks[i][3] * banks[i][2], 2) - std::pow(expected_profit, 2);
        mean += expected_profit;
    }
    variance = sqrt(variance);
    return {mean, variance};
}

// preparation for intermediate steps
void preperation(float a, float b, int n, float sigma, float mu, const std::vector<std::vector<float>>& banks, std::vector<std::vector<float>>& prep) {
    float height = (b - a) / n;
    for (int i = 0; i <= n; i++) {
        prep[i] = X_Z(a + i * height, banks);
        prep[i].push_back(P_Z(a + i * height, mu, sigma));
    }
}

// preparation for intermediate steps of conditional AND/OR
void c_preperation(const std::string type, float a, float b, int n, float sigma, float mu, const std::vector<std::vector<float>>& banks, std::vector<std::vector<float>>& prep, const std::vector<std::vector<float>>& defaults) {
    float height = (b - a) / n;
    std::vector<float> P(n+1);
    if (type == "AND") {
        for (int i = 0; i <= n; i++) {
            prep[i] = X_Z( a + i * height, banks);
            float PZ = P_Z(a + i * height, mu, sigma);
            P[i] = PZ * P_AZ(a + i * height, b, defaults);
        }
    } else if (type == "OR") {
        for (int i = 0; i <= n; i++) {
            float P_AZ = P_AZ_or(a + i * height, b, defaults);
            P[i] = P_Z(a + i * height, mu, sigma) * P_AZ;
            prep[i] = X_Z_or(a + i * height, banks, P_AZ, defaults);
        }
    }

    float sum = P[0] + P[n]; 
    for (int j = 1; j < n; j++) {
        sum += 2 * P[j]; 
    }
    float P_A =  (height / 2) * sum;
    for (int i = 0; i <= n; i++) { 
        prep[i].push_back(P[i]/P_A);
    }
}

// compute the integral
float trapezoidal_integration(float a, float b, int n, float h, const std::vector<std::vector<float>>& prep) {
    float height = (b - a) / n;
    std::vector<float> y(n + 1);
    for (int i = 0; i <= n; i++) {
        float C = prep[i][1] == 0 ? (h - prep[i][0] > 0 ? 1 : 0) : (h - prep[i][0]) / prep[i][1];
        y[i] = 0.5 * erfc(-(C)/sqrt(2)) * prep[i][2];
    }

    float sum = y[0] + y[n]; 
    for (int i = 1; i < n; i++) {
        sum += 2 * y[i]; 
    }
    return (height / 2) * sum;
}

// analytical method for regular situation with no defualts
void analyticBanks_R(const std::string& input_file, std::ofstream& outFile, float alpha) {
    // get the input bank portfolio
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::vector<std::vector<float>> defaults;
    std::string line;

    float upper_limit = 0;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        upper_limit = std::max(upper_limit, bank[8]/bank[5]);
        banks.push_back(bank);
    }

    // set the lower and upper bound of loss
    auto start_time1 = std::chrono::high_resolution_clock::now();
    float lower = 0;
    float upper = 0;
    for (size_t i = 0; i < banks.size(); i++) {
        lower += (banks[i][4] - 1) * banks[i][2];
        upper += banks[i][3] * banks[i][2];
    }

    // set the preperation for integration part
    float lower_limit = 0;
    int interval = 1000;
    std::vector<std::vector<float>> prep(interval+1, std::vector<float>(3));
    preperation(lower_limit, upper_limit, interval, sigma, mu, banks, prep);
    
    // compute f(h) 
    outFile << "F_h,h" << std::endl;
    for (float h = -upper; h < -lower+1000; h += (upper-lower)/7000) {
        float F_h = (h >= -lower) ? 0 : trapezoidal_integration(lower_limit, upper_limit, interval, -h, prep);
        F_h = (F_h < 0) ? 0 : ((F_h > 1) ? 1 : F_h);
        outFile << 1 - F_h << "," << h << std::endl;
    }
    
    // compute the VaR using binary search
    int turn = 0;
    double error = 1e-7;
    while (true && turn < 100) {
        double bound = (upper + lower) / 2;
        float Fh = trapezoidal_integration(lower_limit, upper_limit, interval, bound, prep);
        if (Fh < 1 - alpha + error && Fh > 1 - alpha - error) {
            std::cout << "Analytic Losses VAR at " << alpha * 100 << "% confidence level is: " << -bound << std::endl;
            break;
        } else if (Fh >= 1 - alpha + error) {
            upper = bound;
        } else if (Fh <= 1 - alpha - error) {
            lower = bound;
        }
        turn += 1;
    }

    auto end_time1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed1 = end_time1 - start_time1;
    std::cout << "Time for analytical: " << elapsed1.count() << " seconds" << std::endl;
}

// analytical method for conditional AND
void analyticBanks_AND(const std::string& input_file, std::ofstream& outFile, float alpha, const std::vector<float>& def) {
    // get the input bank portfolio
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::vector<std::vector<float>> defaults;
    std::string line;

    float upper_limit = 1000;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    float P_s = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        if (isDefault(bank[0], def)){
            defaults.push_back(bank);
            upper_limit = std::min(upper_limit, bank[8]/bank[5]);
            P_s += bank[2] * (bank[4] - 1);
        } else {
            banks.push_back(bank);
        }
    }

    // set the lower and upper bound of loss
    auto start_time1 = std::chrono::high_resolution_clock::now();
    float lower = 0;
    float upper = 0;
    for (size_t i = 0; i < banks.size(); i++) {
        lower += (banks[i][4] - 1) * banks[i][2];
        upper += banks[i][3] * banks[i][2];
    }
    for (size_t i = 0; i < defaults.size(); i++) {
        lower += (defaults[i][4] - 1) * defaults[i][2];
    }

    // set the preperation for integration part
    float lower_limit = 0;
    int interval = 1000;
    std::vector<std::vector<float>> prep(interval+1, std::vector<float>(3));
    c_preperation("AND", lower_limit, upper_limit, interval, sigma, mu, banks, prep, defaults);
    
    // compute f(h) 
    outFile << "F_h,h" << std::endl;
    for (float h = -upper; h < -lower+1000; h += (upper-lower)/7000) {
        float F_h = (h + P_s >= -lower) ? 0 : trapezoidal_integration(lower_limit, upper_limit, interval, -h - P_s, prep);
        F_h = (F_h < 0) ? 0 : ((F_h > 1) ? 1 : F_h);
        outFile << 1 - F_h << "," << h << std::endl;
    }

    // compute the VaR using binary search
    int turn = 0;
    double error = 1e-7;
    while (true && turn < 100) {
        double bound = (upper + lower) / 2;
        float Fh = trapezoidal_integration(lower_limit, upper_limit, interval, bound - P_s, prep);
        if (Fh < 1 - alpha + error && Fh > 1 - alpha - error) {
            std::cout << "Conditional-AND Analytic Losses VAR at " << alpha * 100 << "% confidence level is: " << -bound << std::endl;
            break;
        } else if (Fh >= 1 - alpha + error) {
            upper = bound;
        } else if (Fh <= 1 - alpha - error) {
            lower = bound;
        }
        turn += 1;
    }

    auto end_time1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed1 = end_time1 - start_time1;
    std::cout << "Time for analytical: " << elapsed1.count() << " seconds" << std::endl;
}

// analytical method for conditional OR
void analyticBanks_OR(const std::string& input_file, std::ofstream& outFile, float alpha, const std::vector<float>& def) {
    // get the input bank portfolio
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::vector<std::vector<float>> defaults;
    std::string line;

    float upper_limit = 0;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        if (isDefault(bank[0], def)){
            defaults.push_back(bank);
            upper_limit = std::max(upper_limit, bank[8]/bank[5]);
        } else {
            banks.push_back(bank);
        }
    }

    auto start_time1 = std::chrono::high_resolution_clock::now();
    // set the lower and upper bound for h, vector of h values, and vector for F(h)
    float lower = 0;
    float upper = 0;
    for (size_t i = 0; i < banks.size(); i++) {
        lower += (banks[i][4] - 1) * banks[i][2];
        upper += banks[i][3] * banks[i][2];
    }
    for (size_t i = 0; i < defaults.size(); i++) {
        lower += (defaults[i][4] - 1) * defaults[i][2];
    }

    // set the preperation for integration part
    float lower_limit = 0;
    int interval = 1000;
    std::vector<std::vector<float>> prep(interval+1, std::vector<float>(3));
    c_preperation("OR", lower_limit, upper_limit, interval, sigma, mu, banks, prep, defaults); 

    auto end_time1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed1 = end_time1 - start_time1;
    std::cout << "Time for analytical: " << elapsed1.count() << " seconds" << std::endl;
    
    // compute f(h) 
    outFile << "F_h,h" << std::endl;
    for (float h = -upper; h < -lower+1000; h += (upper-lower)/7000) {
        float F_h = (h >= -lower) ? 0 : trapezoidal_integration(lower_limit, upper_limit, interval, -h, prep);
        F_h = (F_h < 0) ? 0 : ((F_h > 1) ? 1 : F_h);
        outFile << 1 - F_h << "," << h << std::endl;
    }

    // compute the VaR using binary search
    int turn = 0;
    double error = 1e-7;
    while (true && turn < 100) {
        double bound = (upper + lower) / 2;
        float Fh = trapezoidal_integration(lower_limit, upper_limit, interval, bound, prep);
        if (Fh < 1 - alpha + error && Fh > 1 - alpha - error) {
            std::cout << "Conditional-OR Analytic Losses VAR at " << alpha * 100 << "% confidence level is: " << -bound << std::endl;
            break;
        } else if (Fh >= 1 - alpha + error) {
            upper = bound;
        } else if (Fh <= 1 - alpha - error) {
            lower = bound;
        }
        turn += 1;
    }
}

#endif
