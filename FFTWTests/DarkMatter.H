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
#include <functional>


class DMParticle
{
	public:
		double x, y, z;
};

class Node
{
	public:
		double mass;
};

void 
generate_fourier_coefficients_for_field_with_spectrum(const int N, const double L, const std::function<double(double)>& func, fftw_complex *Pk);

void
compute_fourier_coefficients_for_psi(const int N, const double L, const fftw_complex *delta_k,
                                     fftw_complex *psi_x_k, fftw_complex *psi_y_k, fftw_complex *psi_z_k);

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


inline void 
checkCondition(bool condition, const char* errorMessage) {
    if (!condition) {
        throw std::runtime_error(errorMessage);
    }
}

inline int single_index_node(const int& i, const int& j, const int& k, const int& N)
{
    return k * (N+1) * (N+1) + j * (N+1) + i;
}

inline void get_global_node_index_and_indices(const int& id, const int& i, const int& j, const int& k, const int& N, int& gid, std::vector<int>& node_indices)
{
	int inode, jnode, knode;
    if(id == 0){
		inode = i; jnode = j; knode = k;
    } else if (id == 1){
		inode = i+1; jnode = j; knode = k;
    } else if (id == 2){
		inode = i+1; jnode = j+1; knode = k;
    } else if (id == 3){
		inode = i; jnode = j+1; knode = k;
    } else if (id == 4){
		inode = i; jnode = j; knode = k+1;
    } else if (id == 5){
		inode = i+1; jnode = j; knode = k+1;
    } else if (id == 6){
		inode = i+1; jnode = j+1; knode = k+1;
    } else if (id == 7){
		inode = i; jnode = j+1; knode = k+1;
    }
	gid = single_index_node(inode,jnode,knode,N);
	node_indices[0] = inode;
	node_indices[1] = jnode;
	node_indices[2] = knode;
	 
}


#endif // DARKMATTER_H
