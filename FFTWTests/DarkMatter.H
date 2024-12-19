#ifndef DARKMATTER_H
#define DARKMATTER_H

#include <iostream>     // For standard input/output
#include <cmath>    // For mathematical functions like sqrt, pow
#include <iomanip>      // For output formatting
#include <limits>       // For handling numeric limits
#include <vector>
#include <fftw3.h>
#include <random>
#include <algorithm>

void plot_DM_particles_vtk();

inline void sort_power_spectrum(std::vector<std::pair<double, double>>& data)
{

 // Sort the vector by the first element of the pair
    std::sort(data.begin(), data.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    });

    // Remove duplicates based on the first element of the pair
    auto it = std::unique(data.begin(), data.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first == b.first;
    });

    data.erase(it, data.end());
}

inline double
power_spectrum(double k)
{
    if(k != 0.0) {
        return 1.0/(k*k);
    } else {
        return 0.0;
    }
}
#endif // DARKMATTER_H
