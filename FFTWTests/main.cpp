#include <DarkMatter.H>

inline double
power_spectrum(double k)
{
    if(k != 0.0) {
        return 1e-8*k;///(1.0 + std::pow(k,4));
    } else {
        return 0.0;
    }
}

// Example main function
int main() {

    
    // Define grid size
    const int N = 256; // Number of points in each dimension
    const double L = 7700.0; // Domain size in each dimension
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
                  //      << "Value = " << data[index][0] << "         " << data[index][1] << std::endl; // Real part only
            }
        }
    }

    // Write to VTK file
    //writeStructuredGridVTK("delta_field.vtk", N, N, N, delta);

	auto it_min = min_element(delta.begin(), delta.end());
	auto it_max = max_element(delta.begin(), delta.end());

	std::cout << "The min and max elements are " << *it_min << " " << *it_max << "\n";
	

	std::vector<std::pair<double, double>> k_and_delta_k_unsrt;
	compute_power_spectrum(N, L, delta, k_and_delta_k_unsrt);
	
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
				//std::cout << "Values are " << psi_x[index][0] << " " << psi_y[index][0] << " " << psi_z[index][0] << "\n";
				psi_mag.emplace_back(std::sqrt(psi_x[index][0]*psi_x[index][0] + psi_y[index][0]*psi_y[index][0] + psi_z[index][0]*psi_z[index][0]));
			}
		}
	}

	auto it = max_element(psi_mag.begin(),psi_mag.end());
	std::cout << "Max value and dx is " << *it << " " << static_cast<double>(L/N) <<  "\n";

	try {
        // Custom assert with throwing an error
        checkCondition(*it < 0.5*dx, "Assertion failed: *it < dx");

        std::cout << "Condition satisfied." << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE); // Exit with failure
    }
				
	// Create mesh of particles and displace
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

				/*dm_particle_tmp.x = xcen + 0.1*psi_x[index][0]; 
				dm_particle_tmp.y = ycen + 0.2*psi_y[index][0]; 
				dm_particle_tmp.z = zcen + 0.8*psi_z[index][0];
				dm_particles.emplace_back(dm_particle_tmp);*/

			}
		}
	}
				
	//plot_DM_particles_vtk(dm_particles);

	// Assign mass to nodes

	const int N_test = N;

	std::vector<Node> nodes;
	nodes.resize((N_test+1) * (N_test+1) * (N_test+1));

	for(auto &v : nodes){
		v.mass = 0.0;
		v.ncontrib = 0;
	}

	const double fac = static_cast<double>((N_test+1e-10)/(N+1e-10));

	std::vector<double> rho(N_test * N_test * N_test, 0.0);
	const double rho_mean = N * N * N/(L * L * L);
	create_cell_centered_rho_field(N_test, dx/fac, dm_particles, nodes, rho_mean, rho);

	std::vector<double> rho_filt(N_test * N_test * N_test, 0.0);
	compute_rho_filtered(N_test, L, dx/fac, rho, rho_filt);
	

	double sum = 0.0;
	for(auto &v : nodes){
		sum += v.mass;
		//std::cout << "nconbtrib is " << v.ncontrib << "\n";
		if(v.ncontrib==0) {
			//std::cout << "There are zero contirbution nodes" << "\n";
		}
		if(v.mass > 1.0){
			 //std::cout << "There are cells with more than one particle " << v.mass  << "\n";
		}
	}

	double rho_avg = 0.0;
	for(auto &v : rho){
		rho_avg += v;
	}

	double delta_mean = 0.0;
	for(auto &v : delta){
		delta_mean += v;
	}
	
	printf("The delta mean is %0.15g\n", delta_mean);
	
	printf("The total mass is %0.15g %0.15g %0.15g %0.15g\n", sum, rho_avg, rho_mean, static_cast<double>(N*N*N));


	std::vector<std::pair<double, double>> k_and_rho_k_unsrt;
	compute_power_spectrum(N_test, L, rho_filt, k_and_rho_k_unsrt);

	FILE* file_k_vs_rho_k;
	file_k_vs_rho_k = fopen("k_vs_rho_k.txt","w");	
	
	for(const auto& p : k_and_rho_k_unsrt) {
		fprintf(file_k_vs_rho_k,"%0.15g %0.15g %0.15g\n", p.first, p.second, power_spectrum(p.first));
        //std::cout << "(" << p.first << ", " << p.second << ")\n";
    }	
	fclose(file_k_vs_rho_k);
	
    return 0;
}
