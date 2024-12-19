#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <random>


double
power_spectrum(double k)
{
	if(k != 0.0) {
		return 1.0/(k*k);
	} else {
		return 0.0;
	}
}

int main() {
    // Define grid size
    const int N = 8; // Number of points in each dimension
    const double L = 1.0; // Domain size in each dimension
    const double dx = L / N; // Grid spacing

    // Allocate memory for the function and its Fourier transform
    fftw_complex *data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *fft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

    // Create the FFTW plan
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, data, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan inverse_plan = fftw_plan_dft_3d(N, N, N, fft_in, data, FFTW_BACKWARD, FFTW_ESTIMATE);

	double total_points = N * N * N;

	// Random number generation for phases (Gaussian distributed)
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normal_dist(0.0, 1.0); // Mean = 0, Stddev = 1

    // Initialize the input data with f(x, y, z) = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
    for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
				if(ky!=0 or kz!=0) {
					//continue;
				}
				int index = kz * N * N + ky * N + kx; // 3D index flattened
	
                if(kx > N/2 or ky > N/2 or kz > N/2) {
					// For real data, the conjugate transpose of the N - kx coefficient should be set
					// The N - kx coefficient will already have been computed
					int kx_tmp = (kx <= N/2) ? kx : N - kx;
					int ky_tmp = (ky <= N/2) ? ky : N - ky;
					int kz_tmp = (kz <= N/2) ? kz : N - kz;
					
					int index_tmp = kz_tmp * N * N + ky_tmp * N + kx_tmp; 
					fft_in[index][0] =      fft_in[index_tmp][0]; // Real part
                	fft_in[index][1] = -1.0*fft_in[index_tmp][1]; // Imaginary part
				} else {
					int kx_phys = (kx <= N/2) ? kx : kx - N;
					int ky_phys = (ky <= N/2) ? ky : ky - N;
					int kz_phys = (kz <= N/2) ? kz : kz - N;

					// Compute the physical wavenumber magnitude
                	double k_mag = std::sqrt(kx_phys * kx_phys + ky_phys * ky_phys + kz_phys * kz_phys);

                	// Get the power spectrum value
	                double P_k = power_spectrum(k_mag);

    	            // Assign random Fourier coefficients with Gaussian-distributed amplitudes
        	        double amplitude = std::sqrt(P_k / 2.0); // Scale by sqrt(P(k)/2) for real+imag contributions
					double rand = normal_dist(generator);
            	    double real_part = amplitude * 0.0;//std::cos(2.0*M_PI*rand);
                	double imag_part = amplitude * -1.0;//std::sin(2.0*M_PI*rand);

                	fft_in[index][0] = real_part; // Real part
                	fft_in[index][1] = imag_part; // Imaginary part
				}
            }
        }
    }

    // Output the results
    std::cout << "Fourier Transform Output:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                double real_part = fft_in[idx][0];
                double imag_part = fft_in[idx][1];
                //std::cout << "k=(" << i << ", " << j << ", " << k << ") : "
                  //        << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }

	// Perform the inverse Fourier transform
    fftw_execute(inverse_plan);

    // Output the reconstructed real-space data
    std::cout << "\nReconstructed Real-Space Field:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                //std::cout << "x=(" << i << ", " << j << ", " << k << ") : "
                  //        << "Value = " << data[idx][0] << "         " << data[idx][1] << std::endl; // Real part only
            }
        }
    }
	
// Perform the Fourier transform
    fftw_execute(forward_plan);

    // Output the results
    std::cout << "Fourier Transform Output:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                double real_part = fft_out[idx][0]/total_points;
                double imag_part = fft_out[idx][1]/total_points;
                //std::cout << "k=(" << i << ", " << j << ", " << k << ") : "
                  //        << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }	

	std::vector<double> Pk_computed;
	std::vector<double> kvec;
	int count = 1;
	for (int kz = 0; kz <= N/2; ++kz) {
        for (int ky = 0; ky <= N/2; ++ky) {
            for (int kx = 0; kx <= N/2; ++kx) {
				int idx = kz * N * N + ky * N + kx; // 3D index flattened
				double real_part = fft_out[idx][0]/total_points;
                double imag_part = fft_out[idx][1]/total_points;
				
				double kmag = std::sqrt(kx*kx + ky*ky + kz*kz);
				kvec.push_back(kmag);
	
				double tmp = 2.0*(real_part*real_part + imag_part*imag_part);

				Pk_computed.push_back(tmp);
				std::cout << count << " " << kmag << " " << power_spectrum(kmag) << " " << tmp << "\n";
				count++;	
			}
		}
	}
					
    // Free memory and clean up FFTW
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(data);
    fftw_free(fft_in);

    return 0;
}

