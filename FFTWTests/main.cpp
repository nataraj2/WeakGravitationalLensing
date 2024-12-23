#include <DarkMatter.H>

inline double
power_spectrum(double k)
{
    if(k != 0.0) {
        return 1e-5/(k*k);
    } else {
        return 0.0;
    }
}

// Example main function
int main() {

    
    // Define grid size
    const int N = 128; // Number of points in each dimension
    const double L = 1.0; // Domain size in each dimension
    const double dx = L / N; // Grid spacing

    // Allocate memory for the function and its Fourier transform
    fftw_complex *Pk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

    // Create the FFTW plan

	double total_points = N * N * N;


	// Random number generation for phases (Gaussian distributed)
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normal_dist(0.0, 1.0); // Mean = 0, Stddev = 1

	
    // Initialize the input data with f(x, y, z) = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
    for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {

				int index = kz * N * N + ky * N + kx; // 3D index flattened
		
				// Compute the physical wavenumber magnitude
				std::vector<double> k_phys(3,0.0);
				std::vector<int> k_phys_int(3,0);

				get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);
				
				// The -k is 
				double kx_minus = -k_phys_int[0];
				double ky_minus = -k_phys_int[1];
				double kz_minus = -k_phys_int[2];
	
                double k_mag = std::sqrt(k_phys[0] * k_phys[0] + k_phys[1] * k_phys[1] + k_phys[2] * k_phys[2]);

                // Get the power spectrum value
	            double Pk_mag = power_spectrum(k_mag);

    	        // Assign random Fourier coefficients with Gaussian-distributed amplitudes
        	    double amplitude = std::sqrt(Pk_mag / 2.0); // Scale by sqrt(P(k)/2) for real+imag contributions
				double rand = normal_dist(generator);
            	double real_part = amplitude * std::cos(2.0*M_PI*rand);
               	double imag_part = amplitude * std::sin(2.0*M_PI*rand);

                Pk[index][0] = real_part; // Real part
                Pk[index][1] = imag_part; // Imaginary part

				if(kx == N/2 or ky == N/2 or kz == N/2){
					Pk[index][0] = amplitude;
					Pk[index][1] = 0.0;
				}

				
				double kx_mir = kx_minus >= 0 ? kx_minus : N+kx_minus;	
				double ky_mir = ky_minus >= 0 ? ky_minus : N+ky_minus;	
				double kz_mir = kz_minus >= 0 ? kz_minus : N+kz_minus;

				if(kx_mir == -N/2)kx_mir = N/2;
				if(ky_mir == -N/2)ky_mir = N/2;
				if(kz_mir == -N/2)kz_mir = N/2;

				int index_mir = kz_mir * N * N + ky_mir * N + kx_mir;

				Pk[index_mir][0] =  Pk[index][0];
				Pk[index_mir][1] = -Pk[index][1];
            }
        }
    }

    // Output the results
    /*for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                int idx = kz * N * N + ky * N + kx; // 3D index flattened
                double real_part = Pk[idx][0];
                double imag_part = Pk[idx][1];
				std::vector<double> k_phys(3,0.0);
				std::vector<int> k_phys_int(3,0);

				get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);

                std::cout << "k=(" << kx_phys << ", " << ky_phys << ", " << kz_phys << ") : "
                         << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }*/


	// Perform the inverse Fourier transform
    fftw_complex *data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	fftw_plan inverse_plan = fftw_plan_dft_3d(N, N, N, Pk, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(inverse_plan);
	fftw_destroy_plan(inverse_plan);

    // Output the reconstructed real-space data
    std::cout << "\nReconstructed Real-Space Field:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                //std::cout << "x=(" << i << ", " << j << ", " << k << ") : "
                  //      << "Value = " << data[idx][0] << "         " << data[idx][1] << std::endl; // Real part only
            }
        }
    }

// Define scalar field values (example: a simple gradient field)
    std::vector<double> delta(N * N  * N, 0.0);
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
                delta[index] = data[index][0]; // Example: radial distance
            }
        }
    }

    // Write to VTK file
    writeStructuredGridVTK("delta_field.vtk", N, N, N, delta);


	
// Perform the Fourier transform
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, data, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);

    // Output the results
    std::cout << "Fourier Transform Output:" << std::endl;
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = k * N * N + j * N + i; // 3D index flattened
                double real_part = fft_out[idx][0]/total_points;
                double imag_part = fft_out[idx][1]/total_points;
                //std::cout << "k=(" << i << ", " << j << ", " << k << ") : "
                  //       << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }	

	std::vector<std::pair<double, double>> k_and_Pk_unsrt;
	int count = 1;
	for (int kz = 0; kz <= N/2; ++kz) {
        for (int ky = 0; ky <= N/2; ++ky) {
            for (int kx = 0; kx <= N/2; ++kx) {
				int idx = kz * N * N + ky * N + kx; // 3D index flattened
				double real_part = fft_out[idx][0]/total_points;
                double imag_part = fft_out[idx][1]/total_points;

				std::vector<double> k_phys(3,0.0);
				std::vector<int> k_phys_int(3,0);

				get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);
                double k_mag = std::sqrt(k_phys[0] * k_phys[0] + k_phys[1] * k_phys[1] + k_phys[2] * k_phys[2]);
				double tmp = 2.0*(real_part*real_part + imag_part*imag_part);
				tmp = std::fabs(tmp) < 1e-14 ? 0.0 : tmp;

				std::pair<double, double> data_tmp = {k_mag, tmp};
				k_and_Pk_unsrt.emplace_back(data_tmp);
				//std::cout << count << " " << kmag << " " << power_spectrum(kmag) << " " << tmp << "\n";
				count++;	
			}
		}
	}

	sort_power_spectrum(k_and_Pk_unsrt);

	std::cout << "Size is " << k_and_Pk_unsrt.size() << "\n";

	FILE* file_k_vs_Pk;
	file_k_vs_Pk = fopen("k_vs_Pk.txt","w");	
	
	for(const auto& p : k_and_Pk_unsrt) {
		fprintf(file_k_vs_Pk,"%0.15g %0.15g %0.15g\n", p.first, p.second, power_spectrum(p.first));
        //std::cout << "(" << p.first << ", " << p.second << ")\n";
    }	
	fclose(file_k_vs_Pk);

	
	// Define the Fourier coefficients for the displacement vector
    fftw_complex *psi_x_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *psi_y_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *psi_z_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    
 // Initialize the input data with f(x, y, z) = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
    for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                int index = kz * N * N + ky * N + kx; // 3D index flattened

				std::vector<double> k_phys(3,0.0);
				std::vector<int> k_phys_int(3,0);

				get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);

				double kx_phys = k_phys[0];
				double ky_phys = k_phys[1];
				double kz_phys = k_phys[2];

                double k_mag = std::sqrt(k_phys[0] * k_phys[0] + k_phys[1] * k_phys[1] + k_phys[2] * k_phys[2]);

				if(k_mag == 0.0){
					psi_x_k[index][0] = 0.0;
					psi_x_k[index][1] = 0.0;

					psi_y_k[index][0] = 0.0;
					psi_y_k[index][1] = 0.0;
					
					psi_z_k[index][0] = 0.0;
					psi_z_k[index][1] = 0.0;
				} else {
                    // For real data, the conjugate transpose of the N - kx coefficient should be set
                    // The N - kx coefficient will already have been computed
                    psi_x_k[index][0] =  -1.0*kx_phys/(k_mag*k_mag)*Pk[index][1]; // Real part
                    psi_x_k[index][1] =       kx_phys/(k_mag*k_mag)*Pk[index][0]; // Imaginary part
					if(k_phys_int[0] == N/2) {
						psi_x_k[index][1] = 0.0;
					}

					psi_y_k[index][0] =  -1.0*ky_phys/(k_mag*k_mag)*Pk[index][1]; // Real part
                    psi_y_k[index][1] =       ky_phys/(k_mag*k_mag)*Pk[index][0]; // Imaginary part
					if(k_phys_int[1] == N/2) {
						psi_y_k[index][1] = 0.0;
					}

					psi_z_k[index][0] =  -1.0*kz_phys/(k_mag*k_mag)*Pk[index][1]; // Real part
                    psi_z_k[index][1] =       kz_phys/(k_mag*k_mag)*Pk[index][0]; // Imaginary part
					if(k_phys_int[2] == N/2) {
						psi_z_k[index][1] = 0.0;
					}
				}
            }
        }
    }
	
	fftw_complex *psi_x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	inverse_plan = fftw_plan_dft_3d(N, N, N, psi_x_k, psi_x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(psi_x_k);

	fftw_complex *psi_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	inverse_plan = fftw_plan_dft_3d(N, N, N, psi_y_k, psi_y, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(psi_y_k);

	fftw_complex *psi_z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	inverse_plan = fftw_plan_dft_3d(N, N, N, psi_z_k, psi_z, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(psi_z_k);


	std::vector<double> psi_mag;
	// Create a grid and dispalce particles 
    for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                int index = kz * N * N + ky * N + kx; // 3D index flattened
				//std::cout << kx << " " << ky << " " << kz << " " << psi_x_k[index][0] << " " << psi_x_k[index][1] << "\n"; 
				//std::cout << kx << " " << ky << " " << kz << " " << psi_z[index][0] << " " << psi_z[index][1] << "\n"; 
				std::cout << "Values are " << psi_x[index][0] << " " << psi_y[index][0] << " " << psi_z[index][0] << "\n";
				psi_mag.emplace_back(std::sqrt(psi_x[index][0]*psi_x[index][0] + psi_y[index][0]*psi_y[index][0] + psi_z[index][0]*psi_z[index][0]));
			}
		}
	}

	auto it = max_element(psi_mag.begin(),psi_mag.end());
	std::cout << "Max value is " << *it << "\n";
				
	std::vector<DMParticle> dm_particles;

	for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
				double xcen = (i+0.5)*L/N;
				double ycen = (j+0.5)*L/N;
				double zcen = (k+0.5)*L/N;
				
				DMParticle dm_particle_tmp;
				dm_particle_tmp.x = xcen + psi_x[index][0]; 
				dm_particle_tmp.y = ycen + psi_y[index][0]; 
				dm_particle_tmp.z = zcen + psi_z[index][0];

				dm_particles.emplace_back(dm_particle_tmp);
			}
		}
	}
				
	// Create mesh of particles and displace
	plot_DM_particles_vtk(dm_particles);
	
    // Free memory and clean up FFTW
    fftw_free(data);
    fftw_free(Pk);

    return 0;
}

