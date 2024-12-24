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
    fftw_complex *delta_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

    // Create the FFTW plan

	double total_points = N * N * N;


	generate_fourier_coefficients_for_field_with_spectrum(N, L, power_spectrum, delta_k);	

	// Perform the inverse Fourier transform
    fftw_complex *data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	fftw_plan inverse_plan = fftw_plan_dft_3d(N, N, N, delta_k, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(inverse_plan);
	fftw_destroy_plan(inverse_plan);

	// Define scalar field values (example: a simple gradient field)
    std::vector<double> delta(N * N  * N, 0.0);
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
                delta[index] = data[index][0]; // Example: radial distance
                //std::cout << "x=(" << i << ", " << j << ", " << k << ") : "
                  //      << "Value = " << data[idx][0] << "         " << data[idx][1] << std::endl; // Real part only
            }
        }
    }

    // Write to VTK file
    writeStructuredGridVTK("delta_field.vtk", N, N, N, delta);
	
	// Perform the Fourier transform
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, data, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);
    fftw_free(data);

    // Output the results
    std::cout << "Fourier Transform Output:" << std::endl;

	std::vector<std::pair<double, double>> k_and_delta_k_unsrt;
	for (int kz = 0; kz <= N/2; ++kz) {
        for (int ky = 0; ky <= N/2; ++ky) {
            for (int kx = 0; kx <= N/2; ++kx) {
				int idx = kz * N * N + ky * N + kx; // 3D index flattened
				double real_part = fft_out[idx][0]/total_points;
                double imag_part = fft_out[idx][1]/total_points;
	             //std::cout << "k=(" << i << ", " << j << ", " << k << ") : "
                  //       << "Real=" << real_part << ", Imag=" << imag_part << std::endl;

				std::vector<double> k_phys(3,0.0);
				std::vector<int> k_phys_int(3,0);

				get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);
                double k_mag = std::sqrt(k_phys[0] * k_phys[0] + k_phys[1] * k_phys[1] + k_phys[2] * k_phys[2]);
				double tmp = 2.0*(real_part*real_part + imag_part*imag_part);
				tmp = std::fabs(tmp) < 1e-14 ? 0.0 : tmp;

				std::pair<double, double> data_tmp = {k_mag, tmp};
				k_and_delta_k_unsrt.emplace_back(data_tmp);
			}
		}
	}

	sort_power_spectrum(k_and_delta_k_unsrt);

	std::cout << "Size is " << k_and_delta_k_unsrt.size() << "\n";

	FILE* file_k_vs_delta_k;
	file_k_vs_delta_k = fopen("k_vs_delta_k.txt","w");	
	
	for(const auto& p : k_and_delta_k_unsrt) {
		fprintf(file_k_vs_delta_k,"%0.15g %0.15g %0.15g\n", p.first, p.second, power_spectrum(p.first));
        //std::cout << "(" << p.first << ", " << p.second << ")\n";
    }	
	fclose(file_k_vs_delta_k);
	
	// Define the Fourier coefficients for the displacement vector
    fftw_complex *psi_x_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *psi_y_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *psi_z_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

	compute_fourier_coefficients_for_psi(N, L, delta_k, psi_x_k, psi_y_k, psi_z_k);

    // Free memory and clean up FFTW
    fftw_free(delta_k);

    
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

	try {
        // Custom assert with throwing an error
        checkCondition(*it < dx, "Assertion failed: *it < dx");

        std::cout << "Condition satisfied." << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE); // Exit with failure
    }
				
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


	// Assign mass to nodes

	std::vector<Node> nodes;
	nodes.resize((N+1) * (N+1) * (N+1));

	for(auto &v : nodes){
		v.mass = 0.0;
	}

	/*for (int n = 0; n < dm_particles.size(); ++n) {
		double xp = dm_particles[n].x;
		double yp = dm_particles[n].y;
		double zp = dm_particles[n].z;
	
		// Cell indices	
		int i = int(xp/dx);
		int j = int(yp/dx);
		int k = int(zp/dx);

		// Distribute the mass to the 8 nodes of the cell based on bilinear interpolation
	
		for(int l=0; l<8; l++) {
					
		}*/	
		

    return 0;
}

