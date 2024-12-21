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

class DMParticle
{
	public:
		double x, y, z;
};


void plot_DM_particles_vtk (std::vector<DMParticle>& dm_particles);

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


#endif // DARKMATTER_H
