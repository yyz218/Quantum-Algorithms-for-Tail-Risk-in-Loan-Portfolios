#ifndef MC_H
#define MC_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <chrono>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/random.hpp>
#include <ctime>


std::vector<float> find_var(std::vector<float>& l, float alpha) {
    std::vector<float> losses = l;
    std::vector<float> out(2);
    std::sort(losses.begin(), losses.end());
    int index = static_cast<int>(losses.size() * (1-alpha));
    index = std::max(0, index - 1); 
    out[0] = losses[index];
    float sum = 0;
    float count = 0;
    for (int i = index; i < losses.size(); i++) {
        sum += losses[i];
        count += 1;
    }
    out[1] = sum / count;
    return out;
}

bool isDefault(float id, const std::vector<float>& defaults) {
    return std::find(defaults.begin(), defaults.end(), id) != defaults.end();
}

void simulateBank_MC(const std::vector<std::vector<float>>& banks, int y, float delta_Z, std::vector<float>& losses, std::mt19937& generator, std::lognormal_distribution<float>& distribution_i) {
    for (size_t i = 0; i < banks.size(); i++) {
        std::lognormal_distribution<float> distribution_i(banks[i][6], banks[i][7]);
        float delta_Z_i = banks[i][9] * distribution_i(generator);
        float delta_v_i = (1 - banks[i][5]) * delta_Z_i + banks[i][5] * delta_Z;
        if (banks[i][1] * delta_v_i > banks[i][8]) {
            losses[y] += -banks[i][2] * banks[i][3];
        } else {
            losses[y] += banks[i][2] * (1 - banks[i][4]);
        }
    }
}

void simulateBanks_MC(const std::string& input_file, const int n, std::ofstream& outFile, float alpha) {
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    std::vector<float> losses(n, 0.0);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        banks.push_back(bank);
    }
    
    // set the random generator and lognormal distributions
    std::random_device rd;
    std::mt19937 generator(rd());
    std::lognormal_distribution<float> distribution(mu, sigma);
    std::lognormal_distribution<float> distribution_i;
    
    // run the Monte Carol simulation for n times
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int y = 0; y < n; y++) {
        float delta_Z = Z * distribution(generator);
        simulateBank_MC(banks, y, delta_Z, losses, generator, distribution_i);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed = end_time - start_time;
    std::sort(losses.begin(), losses.end());
    int index = static_cast<int>(losses.size() * (1-alpha));
    index = std::max(0, index - 1); 

    std::map<double, int> prob;
    for (int i = 0; i < losses.size(); i++) {
        prob[losses[i]]++;
    } 
    outFile << "loss,count" << std::endl;
    for (std::map<double, int>::const_iterator it = prob.begin(); it != prob.end(); it++) {
        double element = it->first;
        int count = it->second;
        outFile << element << "," << count << "\n";
    }
    
    std::cout << "Time for simulation is " << elapsed.count() << std::endl;
    std::cout << "Losses VAR is " << losses[index] << std::endl;
}

int simulateBank_AND(const std::vector<std::vector<float>>& banks, int y, float delta_Z, std::vector<float>& losses, std::mt19937& generator, std::lognormal_distribution<float>& distribution_i) {
    for (size_t i = 0; i < banks.size(); i++) {
        std::lognormal_distribution<float> distribution_i(banks[i][6], banks[i][7]);
        float delta_Z_i = banks[i][9] * distribution_i(generator);
        float delta_v_i = (1 - banks[i][5]) * delta_Z_i + banks[i][5] * delta_Z;
        if (banks[i][1] * delta_v_i > banks[i][8] && banks[i][10] == 1) {
            losses[y] = 0;
            return y;
        }
        if (banks[i][1] * delta_v_i > banks[i][8]) {
            losses[y] += -banks[i][2] * banks[i][3];
        } else {
            losses[y] += banks[i][2] * (1 - banks[i][4]);
        }
    }
    return y+1;
}

void simulateBanks_AND(const std::string& input_file, const int n, std::ofstream& outFile, float alpha, const std::vector<float>& defaults) {
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    std::vector<float> losses(n, 0.0);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        if (isDefault(bank[0], defaults)){
            bank.push_back(1);
        } else {
            bank.push_back(0);
        }
        banks.push_back(bank);
    }
    // set the random generator and lognormal distributions
    std::random_device rd;
    std::mt19937 generator(rd());
    std::lognormal_distribution<float> distribution(mu , sigma);
    std::lognormal_distribution<float> distribution_i;

    std::sort(banks.begin(), banks.end(), [](const std::vector<float>& a, const std::vector<float>& b) {
        return a.back() > b.back();
    });

    // run the Monte Carol simulation for n times
    int id = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int y = 0; y < n; y++) {
        float delta_Z = Z * distribution(generator);
        id = simulateBank_AND(banks, id, delta_Z, losses, generator, distribution_i);
    }

    losses.erase(std::remove(losses.begin(), losses.end(), 0), losses.end());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed = end_time - start_time;
    std::sort(losses.begin(), losses.end());
    int index = static_cast<int>(losses.size() * (1-alpha));
    index = std::max(0, index - 1); 

    std::map<double, int> prob;
    for (int i = 0; i < losses.size(); i++) {
        prob[losses[i]]++;
    } 
    outFile << "loss,count" << std::endl;
    for (std::map<double, int>::const_iterator it = prob.begin(); it != prob.end(); it++) {
        double element = it->first;
        int count = it->second;
        outFile << element << "," << count << "\n";
    }
    
    std::cout << "Time for simulation is " << elapsed.count() << std::endl;
    std::cout << "Conditional-AND Losses VAR is " << losses[index] << std::endl;
}

int simulateBank_OR(const std::vector<std::vector<float>>& banks, int y, float delta_Z, std::vector<float>& losses, std::mt19937& generator, std::lognormal_distribution<float>& distribution_i) {
    bool test = false;
    for (int i = 0; i < banks.size(); i++) {
        std::lognormal_distribution<float> distribution_i(banks[i][6], banks[i][7]);
        float delta_Z_i = banks[i][9] * distribution_i(generator);
        float delta_v_i = (1 - banks[i][5]) * delta_Z_i + banks[i][5] * delta_Z;
        if (banks[i][1] * delta_v_i > banks[i][8] && banks[i][10] == 1) {
            if (test == false && banks[i+1][10] == 0) { losses[y] = 0; return y; }
            losses[y] += -banks[i][2] * banks[i][3];
        } else if (banks[i][2] * delta_v_i <= banks[i][8] && banks[i][10] == 1) {
            losses[y] += banks[i][2] - banks[i][4] * banks[i][2];
            test = true;
        } else if ((banks[i][1] * delta_v_i > banks[i][8] && banks[i][10] == 0)) {
            losses[y] += -banks[i][2] * banks[i][3];
        } else if ((banks[i][1] * delta_v_i <= banks[i][8] && banks[i][10] == 0)) {
            losses[y] += banks[i][2] - banks[i][4] * banks[i][2];
        }
    }
    return y+1;
}

void simulateBanks_OR(const std::string& input_file, const int n, std::ofstream& outFile, float alpha, const std::vector<float>& defaults) {
    std::ifstream file(input_file);
    std::vector<std::vector<float>> banks;
    std::string line;
    std::getline(file, line);
    std::istringstream is(line);
    std::vector<float> market(std::istream_iterator<float>(is), {});
    float Z = market[0];
    float mu = market[1];
    float sigma = market[2];
    std::vector<float> losses(n, 0.0);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> bank(std::istream_iterator<float>(iss), {});
        if (isDefault(bank[0], defaults)){
            bank.push_back(1);
        } else {
            bank.push_back(0);
        }
        banks.push_back(bank);
    }

    // set the random generator and lognormal distributions
    std::random_device rd;
    std::mt19937 generator(rd());
    std::lognormal_distribution<float> distribution(mu, sigma);
    std::lognormal_distribution<float> distribution_i;

    std::sort(banks.begin(), banks.end(), [](const std::vector<float>& a, const std::vector<float>& b) {
        return a.back() > b.back();
    });

    // run the Monte Carol simulation for n times
    int id = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int y = 0; y < n; y++) {
        float delta_Z = Z * distribution(generator);
        id = simulateBank_OR(banks, id, delta_Z, losses, generator, distribution_i);
    }

    losses.erase(std::remove(losses.begin(), losses.end(), 0), losses.end());
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed = end_time - start_time;
    std::sort(losses.begin(), losses.end());
    int index = static_cast<int>(losses.size() * (1-alpha));
    index = std::max(0, index - 1); 

    std::map<double, int> prob;
    for (int i = 0; i < losses.size(); i++) {
        prob[losses[i]]++;
    } 
    outFile << "loss,count" << std::endl;
    for (std::map<double, int>::const_iterator it = prob.begin(); it != prob.end(); it++) {
        double element = it->first;
        int count = it->second;
        outFile << element << "," << count << "\n";
    }
    
    std::cout << "Time for simulation is " << elapsed.count() << std::endl;
    std::cout << "Conditional-OR Losses VAR is " << losses[index] << std::endl;
}

#endif