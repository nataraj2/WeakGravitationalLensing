#include <DarkMatter.H>

void 
generate_fourier_coefficients_for_field_with_spectrum(const int N, const double L, const std::function<double(double)>& power_spectrum_func, fftw_complex *delta_k)
{
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
                double delta_k_mag = power_spectrum_func(k_mag);

                // Assign random Fourier coefficients with Gaussian-distributed amplitudes
                double amplitude = std::sqrt(delta_k_mag / 2.0); // Scale by sqrt(P(k)/2) for real+imag contributions
                double rand = normal_dist(generator);
                double real_part = amplitude * std::cos(2.0*M_PI*rand);
                double imag_part = amplitude * std::sin(2.0*M_PI*rand);

                delta_k[index][0] = real_part; // Real part
                delta_k[index][1] = imag_part; // Imaginary part

                if(kx == N/2 or ky == N/2 or kz == N/2){
                    delta_k[index][0] = amplitude;
                    delta_k[index][1] = 0.0;
                }


                double kx_mir = kx_minus >= 0 ? kx_minus : N+kx_minus;
                double ky_mir = ky_minus >= 0 ? ky_minus : N+ky_minus;
                double kz_mir = kz_minus >= 0 ? kz_minus : N+kz_minus;

                if(kx_mir == -N/2)kx_mir = N/2;
                if(ky_mir == -N/2)ky_mir = N/2;
                if(kz_mir == -N/2)kz_mir = N/2;

                int index_mir = kz_mir * N * N + ky_mir * N + kx_mir;

                delta_k[index_mir][0] =  delta_k[index][0];
                delta_k[index_mir][1] = -delta_k[index][1];
            }
        }
    }

	 // Output the results
    /*for (int kz = 0; kz < N; ++kz) {
        for (int ky = 0; ky < N; ++ky) {
            for (int kx = 0; kx < N; ++kx) {
                int idx = kz * N * N + ky * N + kx; // 3D index flattened
                double real_part = delta_k[idx][0];
                double imag_part = delta_k[idx][1];
                std::vector<double> k_phys(3,0.0);
                std::vector<int> k_phys_int(3,0);

                get_physical_wavenumber(L, N, kx, ky, kz, k_phys, k_phys_int);

                std::cout << "k=(" << kx_phys << ", " << ky_phys << ", " << kz_phys << ") : "
                         << "Real=" << real_part << ", Imag=" << imag_part << std::endl;
            }
        }
    }*/
}

void 
compute_fourier_coefficients_for_psi(const int N, const double L, const fftw_complex *delta_k, 
									 fftw_complex *psi_x_k, fftw_complex *psi_y_k, fftw_complex *psi_z_k)
{

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
                    psi_x_k[index][0] =  -1.0*kx_phys/(k_mag*k_mag)*delta_k[index][1]; // Real part
                    psi_x_k[index][1] =       kx_phys/(k_mag*k_mag)*delta_k[index][0]; // Imaginary part
                    if(k_phys_int[0] == N/2) {
                        psi_x_k[index][1] = 0.0;
                    }

                    psi_y_k[index][0] =  -1.0*ky_phys/(k_mag*k_mag)*delta_k[index][1]; // Real part
                    psi_y_k[index][1] =       ky_phys/(k_mag*k_mag)*delta_k[index][0]; // Imaginary part
                    if(k_phys_int[1] == N/2) {
                        psi_y_k[index][1] = 0.0;
                    }

                    psi_z_k[index][0] =  -1.0*kz_phys/(k_mag*k_mag)*delta_k[index][1]; // Real part
                    psi_z_k[index][1] =       kz_phys/(k_mag*k_mag)*delta_k[index][0]; // Imaginary part
                    if(k_phys_int[2] == N/2) {
                        psi_z_k[index][1] = 0.0;
                    }

					double k_mag_int = std::sqrt(k_phys_int[0] * k_phys_int[0] + k_phys_int[1] * k_phys_int[1] + k_phys_int[2] * k_phys_int[2]);
                }
            }
        }
    }
}


void compute_power_spectrum(const int&  N, const double& L, const std::vector<double>& field, 
							std::vector<std::pair<double, double>>& k_and_field_k_unsrt)
{

	fftw_complex *data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	fftw_complex *fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * N);
	
	 for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
                data[index][0] = field[index];
                data[index][1] = 0.0;
            }
        }
    }

	 // Perform the Fourier transform
    fftw_plan forward_plan = fftw_plan_dft_3d(N, N, N, data, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan);
    fftw_destroy_plan(forward_plan);
    fftw_free(data);

	int total_points = N * N * N;
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
                //tmp = std::fabs(tmp) < 1e-14 ? 0.0 : tmp;

                std::pair<double, double> data_tmp = {k_mag, tmp};
                k_and_field_k_unsrt.emplace_back(data_tmp);
            }
        }
    }

    sort_power_spectrum(k_and_field_k_unsrt);
    fftw_free(fft_out);
}

void 
create_cell_centered_rho_field(const int& N, const double& dx, const std::vector<DMParticle>& dm_particles,
							   std::vector<Node>& nodes, const double& rho_mean, std::vector<double>& rho)
{

	std::vector<int>node_indices(3,0);
	std::vector<int> num_part_in_cell(N * N * N,0);
    int node_index;
    for (int n = 0; n < dm_particles.size(); ++n) {
        double xp = dm_particles[n].x;
        double yp = dm_particles[n].y;
        double zp = dm_particles[n].z;

        // Cell indices
        int i = int(xp/dx);
        int j = int(yp/dx);
        int k = int(zp/dx);

		int index = k * N * N + j * N + i;
		num_part_in_cell[index] += 1;
	
        // Distribute the mass to the 8 nodes of the cell based on bilinear interpolation
        double sum_weights = 0.0;
		std::vector<double> weight_vec(8,0.0);
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
		if(std::fabs(sum_weights-1.0) > 1e-10){
				std::cout << "The sum of weights is not one " << sum_weights <<  " " << i << " " << j << " " << k << "\n";
				//for(auto &v: weight_vec){
				//	std::cout << v << "\n";
				//}
				exit(0);
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

            nodes[node_index].mass = nodes[node_index].mass + weight;///sum_weights;
            nodes[node_index].ncontrib += 1;
        }
    }

    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
                for(int id=0; id<8; id++) {
                    get_global_node_index_and_indices(id, i, j, k, N, node_index, node_indices);
					if(nodes[node_index].mass!=0.0){
                    	rho[index] += nodes[node_index].mass/nodes[node_index].ncontrib*num_part_in_cell[index];
					}
                }
                rho[index] = (rho[index]/(dx*dx*dx) - rho_mean)/rho_mean;
                //rho[index] = rho[index]/(dx*dx*dx);
            }
        }
    }
}


void
compute_rho_filtered(const int& N, const double& L, const double& dx, 
					 const std::vector<double>& rho, std::vector<double>& rho_filt)
{
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

                double sigma_k = 1.0/(2.0*M_PI*0.2*dx);
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

    fftw_plan inverse_plan = fftw_plan_dft_3d(N, N, N, rho_filt_k, rho_filt_data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(rho_filt_k);

    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                int index = k * N * N + j * N + i;
                rho_filt[index] = rho_filt_data[index][0];
            }
        }
    }
}







