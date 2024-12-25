#include <DarkMatter.H>

inline double
power_spectrum(double k)
{
    if(k != 0.0) {
        return 1e-8*k/(1.0 + std::pow(k,4));
    } else {
        return 0.0;
    }
}

// Example main function
int main() {

    
    // Define grid size
    const int N = 768; // Number of points in each dimension
    const double L = 2.0*M_PI; // Domain size in each dimension
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
	std::cout << "Max value is " << *it << "\n";

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
			}
		}
	}
				
	//plot_DM_particles_vtk(dm_particles);

	// Assign mass to nodes

	std::vector<Node> nodes;
	nodes.resize((N+1) * (N+1) * (N+1));

	for(auto &v : nodes){
		v.mass = 0.0;
		v.ncontrib = 0;
	}

	std::vector<int>node_indices(3,0);
	int node_index;
	for (int n = 0; n < dm_particles.size(); ++n) {
		double xp = dm_particles[n].x;
		double yp = dm_particles[n].y;
		double zp = dm_particles[n].z;
	
		// Cell indices	
		int i = int(xp/dx);
		int j = int(yp/dx);
		int k = int(zp/dx);

		// Distribute the mass to the 8 nodes of the cell based on bilinear interpolation
		double sum_weights = 0.0;
		for(int id=0; id<8; id++) {
			get_global_node_index_and_indices(id, i, j, k, N, node_index, node_indices);
			// Find the trilinear interpolation coefficetns
			double xnode = node_indices[0]*dx;
			double ynode = node_indices[1]*dx;
			double znode = node_indices[2]*dx;

			double wx, wy, wz;
			double weight = get_interp_coefficients_trilinear(xp, yp, zp, 
									xnode, ynode, znode,
									wx, wy, wz,
									dx);
			/*double weight = get_interp_coefficients_cubic_spline(xp, yp, zp, 
									xnode, ynode, znode,
									wx, wy, wz,
									dx); */

				
			sum_weights += weight;
		}
		for(int id=0; id<8; id++) {
            get_global_node_index_and_indices(id, i, j, k, N, node_index, node_indices);
			double xnode = node_indices[0]*dx;
            double ynode = node_indices[1]*dx;
            double znode = node_indices[2]*dx;

            double wx, wy, wz;
            double weight = get_interp_coefficients_trilinear(xp, yp, zp,
                                    xnode, ynode, znode,
                                    wx, wy, wz,
                                    dx); 
			/*double weight = get_interp_coefficients_cubic_spline(xp, yp, zp, 
									xnode, ynode, znode,
									wx, wy, wz,
									dx); */

			nodes[node_index].mass = nodes[node_index].mass + weight/sum_weights;
			nodes[node_index].ncontrib += 1;	
		}
	}

	std::vector<double> rho(N * N * N, 0.0);
	for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
				for(int id=0; id<8; id++) {
					get_global_node_index_and_indices(id, i, j, k, N, node_index, node_indices);
					rho[index] += nodes[node_index].mass*int(8/nodes[node_index].ncontrib);
				}
				rho[index] = (rho[index]/(8.0) - 1.0)/1.0;
			}
		}
	}

    //writeStructuredGridVTK("rho_field.vtk", N, N, N, rho);


	fftw_complex *rho_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
    fftw_complex *rho_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	

     for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
			
                rho_data[index][0] = rho[index];
                rho_data[index][1] = 0.0;
            }
        }
    }

     // Perform the Fourier transform
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, rho_data, rho_k, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);
    fftw_free(rho_data);

	fftw_complex *rho_filt_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	fftw_complex *rho_filt_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);

	 for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                int index = kz * N * N + ky * N + kx;
			
				std::vector<double> k_phys(3,0.0);
                std::vector<int> k_phys_int(3,0);

                get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);
			
				double k_mag = std::sqrt(k_phys[0] * k_phys[0] + k_phys[1] * k_phys[1] + k_phys[2] * k_phys[2]);

				double sigma_k = 1.0/(2.0*M_PI*dx);
				double Gfac = std::exp(-k_mag*k_mag/(2.0*sigma_k*sigma_k));
				//double Gfac = 0.5*(1.0 + std::tanh((sigma_k - k_mag)/(0.1*N/2*2*M_PI/L)));
				
				//double Gfac = 1.0/(1.0 + std::pow(k_mag/sigma_k,3));
				//double Gfac = std::exp(-k_mag/sigma_k);
				double width = 3.0*dx;
				
				//double Gfac = 3.0/std::pow(k_mag*width,3)*(std::sin(k_mag*width) - k_mag*width*std::cos(k_mag*width));
				if(k_mag == 0){
					Gfac = 1.0;
				}

                //rho_filt_k[index][0] = rho_k[index][0]*Gfac/(N * N * N);
                //rho_filt_k[index][1] = rho_k[index][1]*Gfac/(N * N * N);

				double W_CIC_x = std::sin(M_PI*k_phys_int[0]/(2.0*N/2.0))/(M_PI*k_phys_int[0]/(2.0*N/2.0));
				double W_CIC_y = std::sin(M_PI*k_phys_int[1]/(2.0*N/2.0))/(M_PI*k_phys_int[1]/(2.0*N/2.0));
				double W_CIC_z = std::sin(M_PI*k_phys_int[2]/(2.0*N/2.0))/(M_PI*k_phys_int[2]/(2.0*N/2.0));
			
				if(k_phys_int[0] == 0.0){
					W_CIC_x = 1.0;
				}
				if(k_phys_int[1] == 0.0){
					W_CIC_y = 1.0;
				}
				if(k_phys_int[2] == 0.0){
					W_CIC_z = 1.0;
				}
				double W_CIC = W_CIC_x * W_CIC_y * W_CIC_z;	

				rho_filt_k[index][0] = rho_k[index][0]/(N * N * N);
                rho_filt_k[index][1] = rho_k[index][1]/(N * N * N);

				//rho_filt_k[index][0] = rho_filt_k[index][0]/(W_CIC*W_CIC);
                //rho_filt_k[index][1] = rho_filt_k[index][1]/(W_CIC*W_CIC);

				rho_filt_k[index][0] = rho_filt_k[index][0]*Gfac;
                rho_filt_k[index][1] = rho_filt_k[index][1]*Gfac;

				double k_mag_int = std::sqrt(k_phys_int[0] * k_phys_int[0] + k_phys_int[1] * k_phys_int[1] + k_phys_int[2] * k_phys_int[2]);

				if(k_mag_int > 2.0/3.0*N/2.0) {
					//rho_filt_k[index][0] = 0.0;
					//rho_filt_k[index][1] = 0.0;
				}
            }
        }
    }

    inverse_plan = fftw_plan_dft_3d(N, N, N, rho_filt_k, rho_filt_data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
	fftw_free(rho_filt_k);

	std::vector<double> rho_filt(N * N * N, 0.0);
	for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
				rho_filt[index] = rho_filt_data[index][0];
			}
		}
	}


	double sum = 0.0;
	for(auto &v : nodes){
		sum += v.mass;
	}

	double rho_mean = 0.0;
	for(auto &v : rho){
		rho_mean += v;
	}

	double delta_mean = 0.0;
	for(auto &v : delta){
		delta_mean += v;
	}
	
	printf("The delta mean is %0.15g\n", delta_mean);
	
	printf("The total mass is %0.15g %0.15g %0.15g\n", sum, rho_mean, static_cast<double>(N*N*N));


	std::vector<std::pair<double, double>> k_and_rho_k_unsrt;
	compute_power_spectrum(N, L, rho_filt, k_and_rho_k_unsrt);

	FILE* file_k_vs_rho_k;
	file_k_vs_rho_k = fopen("k_vs_rho_k.txt","w");	
	
	for(const auto& p : k_and_rho_k_unsrt) {
		fprintf(file_k_vs_rho_k,"%0.15g %0.15g %0.15g\n", p.first, p.second, power_spectrum(p.first));
        //std::cout << "(" << p.first << ", " << p.second << ")\n";
    }	
	fclose(file_k_vs_rho_k);
	
    return 0;
}
