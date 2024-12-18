#include <iostream>
#include <cmath>
#include <fftw3.h>

int main() {
    // Define grid size
    const int N = 8; // Number of points in each dimension
    const double L = 1.0; // Domain size in each dimension
    const double dx = L / N; // Grid spacing

    // Allocate memory for the function and its Fourier transform
    fftw_complex *data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

    // Create the FFTW plan
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, data, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan inverse_plan = fftw_plan_dft_3d(N, N, N, fft_out, data, FFTW_BACKWARD, FFTW_ESTIMATE);


    // Initialize the input data with f(x, y, z) = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
    for (int k = 0; k < N; ++k) {
        double z = k * dx; // z-coordinate
        for (int j = 0; j < N; ++j) {
            double y = j * dx; // y-coordinate
            for (int i = 0; i < N; ++i) {
                double x = i * dx; // x-coordinate
                int idx = k * N * N + j * N + i; // 3D index flattened
                data[idx][0] = sin(2 * M_PI * (x+y));// + sin(2 * M_PI * y);//* sin(2 * M_PI * z); // Real part
                data[idx][1] = 0.0; // Imaginary part
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
                double real_part = fft_out[idx][0];
                double imag_part = fft_out[idx][1];
                std::cout << "k=(" << i << ", " << j << ", " << k << ") : "
                          << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }


	// Perform the inverse Fourier transform
    fftw_execute(inverse_plan);

    // Normalize the inverse transform output
    int total_points = N * N * N;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                data[idx][0] /= total_points; // Normalize real part
                data[idx][1] /= total_points; // Normalize imaginary part
            }
        }
    }

    // Output the reconstructed real-space data
    std::cout << "\nReconstructed Real-Space Field:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                std::cout << "x=(" << i << ", " << j << ", " << k << ") : "
                          << "Value=" << data[idx][0] << std::endl; // Real part only
            }
        }
    }

    // Free memory and clean up FFTW
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(data);
    fftw_free(fft_out);

    return 0;
}

