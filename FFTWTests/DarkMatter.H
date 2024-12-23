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
#include <fstream>

class DMParticle
{
	public:
		double x, y, z;
};

void plot_DM_particles_vtk (std::vector<DMParticle>& dm_particles);

void writeStructuredGridVTK(const std::string& filename, int nx, int ny, int nz, const std::vector<double>& scalarField);

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

inline void get_physical_wavenumber(const double& L, const int N, 
									const double& kx, const int &ky, const int& kz, 
									std::vector<double>& k_phys, std::vector<int>& k_phys_int)
{
	k_phys_int[0] = kx <= N/2 ? kx : kx-N;
    k_phys_int[1] = ky <= N/2 ? ky : ky-N;
    k_phys_int[2] = kz <= N/2 ? kz : kz-N;

	k_phys[0] = 2.0*M_PI/L*k_phys_int[0];
	k_phys[1] = 2.0*M_PI/L*k_phys_int[1];
	k_phys[2] = 2.0*M_PI/L*k_phys_int[2];
}
#endif // DARKMATTER_H
