#include <iostream>
#include <fftw3.h>
#include <cmath>

int main() {
    const int N = 8; // Number of points in each dimension
    const double L = 1.0; // Domain size
    const double dx = L / N; // Grid spacing
    const double k_factor = 2 * M_PI / L; // Factor to calculate wave numbers

    // Allocate arrays for input (real) and output (complex)
    double *input = new double[N * N * N];
    fftw_complex *output = new fftw_complex[N * N * N];

    // Fill the input array with sin(2 * pi * x)
    for (int x = 0; x < N; ++x) {
        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                double X = x * dx;
                input[x * N * N + y * N + z] = sin(2 * M_PI * X);
            }
        }
    }

    // Create a plan for forward FFT
    fftw_plan plan_forward = fftw_plan_dft_r2c_3d(N, N, N, input, output, FFTW_ESTIMATE);
    
    // Execute forward FFT
    fftw_execute(plan_forward);
    
    // Output the result of the forward FFT (real and imaginary parts)
    std::cout << "Forward FFT results (kx, ky, kz, real, imag):" << std::endl;
    for (int kx = 0; kx < N; ++kx) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kz = 0; kz < N; ++kz) {
                int index = kx * N * N + ky * N + kz;

                // Print the wave numbers and real/imaginary parts
                std::cout << kx << ", " << ky << ", " << kz << ", "
                          << output[index][0] << ", " << output[index][1] << std::endl;
            }
        }
    }

    // Create a plan for inverse FFT
    fftw_plan plan_inverse = fftw_plan_dft_c2r_3d(N, N, N, output, input, FFTW_ESTIMATE);

    // Execute inverse FFT
    fftw_execute(plan_inverse);
    
    // Normalize the output from the inverse FFT
    for (int i = 0; i < N * N * N; ++i) {
        input[i] /= (N * N * N);
    }

    // Output the result of the inverse FFT
    std::cout << "Inverse FFT results:" << std::endl;
    for (int x = 0; x < N; ++x) {
        for (int y = 0; y < N; ++y) {
            for (int z = 0; z < N; ++z) {
                //std::cout << input[x * N * N + y * N + z] << " ";
            }
            //std::cout << std::endl;
        }
        //std::cout << std::endl;
    }

    // Clean up
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
    delete[] input;
    delete[] output;

    return 0;
}

