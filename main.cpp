#include <WeakLensing.H>     

int main()
{

	// Constants

	double G_val, M_val, c_val, D_LS, D_OS, D_OL;

	G_val = 6.6743e-11;
	M_val = 1e12*1.989e30;
	c_val = 3e8;
	D_OL = 1.0;
	D_LS = 0.5;
	D_OS = 1.5;

	double D_OL_m = convert_Gpc_to_m(D_OL);
	double D_LS_m = convert_Gpc_to_m(D_LS);
	double D_OS_m = convert_Gpc_to_m(D_OS);

	// z is the direction of shoot axis
	// Shoot till z = z_lens

	int nlens = 1;
	std::vector<double> dist_lens;
	dist_lens.resize(nlens+1,0.0);
	dist_lens[0]= convert_Gpc_to_m(D_OL);
	dist_lens[1]= convert_Gpc_to_m(D_LS);

	double theta_e = std::pow(4.0*G_val*M_val/(c_val*c_val)*D_LS_m/(D_OL_m*D_OS_m),0.5);
	std::cout << "Value of theta_e is " << theta_e << " " << theta_e*180/M_PI*3600;
	//exit(0);

	// Create rays from origin, with the rays being on the surface of a cone and 
	// the for loop increments the half-angle of the cone in each iteration thus covering
	// a large angular diameter

	// Observer is at the origin

	int num_angles_theta = 1000;
	int num_angles_phi = 200;

	int num_rays = num_angles_theta*num_angles_phi;
	std::vector<PolyRay> polyrayvec;
	polyrayvec.resize(num_rays);
	
	double xloc=0.0, yloc=0.0, zloc=0.0;

	double dir[3];

	double theta_max = 2e-5;
	
	// Initialize ray bundles - Vary azimuthal angle phi 
	
	int polyray_count = 0;
	
	for (int itheta=0;itheta<num_angles_theta;itheta++){
		double theta = (itheta+1)*theta_max/num_angles_theta;
		for (int iphi=0;iphi<num_angles_phi;iphi++) {
			polyrayvec[polyray_count].theta_start = theta;
			double phi = iphi*2.0*M_PI/(num_angles_phi-1.0);
			polyrayvec[polyray_count].phi_start = phi;
			dir[0] = std::sin(theta)*std::cos(phi);
			dir[1] = std::sin(theta)*std::sin(phi);
			dir[2] = std::cos(theta);

			polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
			polyray_count++;
		}
	}

	
	// Deflection caused by the lens is given by 
	// alpha = 4*G*M/(c*c*xi) where xi is the distance of the ray on
	// the lens plane	


	// All rays start at the observer	
	for (int i=0;i<num_rays;i++) {
		polyrayvec[i].coord.emplace_back(xloc, yloc, zloc);
	}


	for (int iray=0;iray<num_rays;iray++){
		// Find the intersection on the lens plane
		double x = polyrayvec[iray].coord[0][0] + dist_lens[0]*polyrayvec[iray].dir[0][0];
        double y = polyrayvec[iray].coord[0][1] + dist_lens[0]*polyrayvec[iray].dir[0][1];

		double xi = std::pow(x*x+y*y,0.5);
		double alpha = 4.0*G_val*M_val/(c_val*c_val*xi);

        double theta = polyrayvec[iray].theta_start - alpha;
        double phi = polyrayvec[iray].phi_start;
        dir[0] = std::sin(theta)*std::cos(phi);
        dir[1] = std::sin(theta)*std::sin(phi);
        dir[2] = std::cos(theta);

        polyrayvec[iray].dir.emplace_back(dir[0], dir[1], dir[2]);
    }
	
	// Do plotting in Gpc units
	
	dist_lens[0]= D_OL;
	dist_lens[1]= D_LS;

	// Find the end of each ray section of the poly ray
	// The ray direction remains the same between lenses
	// Hence if there are n lenses, there will be n+1 coordinates
	// for eg. there is 1 lens there will be 2 coordinates

	for(int iray=0;iray<num_rays;iray++){
		for(int i=1;i<=nlens+1;i++) {	
			double x = polyrayvec[iray].coord[i-1][0] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][0];
			double y = polyrayvec[iray].coord[i-1][1] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][1];
			double z = polyrayvec[iray].coord[i-1][2] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][2];
			polyrayvec[iray].coord.emplace_back(x,y,z);	
	 	}
	}

	// Does the end of the polyray intersect the source
	
	std::vector<PolyRay> final_polyrayvec;
	for(int iray=0;iray<num_rays;iray++){
		double x = polyrayvec[iray].coord[nlens+1][0];
		double y = polyrayvec[iray].coord[nlens+1][1];
		//double z = polyrayvec[iray].coord[nlens+1][2];
		std::cout << "Value of " << iray << " " << std::pow(x*x + y*y,0.5) << "\n";
		if(std::pow(x*x + y*y,0.5) < 1e-7){
			std::cout << "Inside here " << std::pow(x*x + y*y,0.5) << "\n";
			//exit(0);
			final_polyrayvec.push_back(polyrayvec[iray]);
		}
	}

	write_polyray_bundle(polyrayvec,"polyray_bundle.vtk");
	if(final_polyrayvec.size() > 0) {	
		write_polyray_bundle(final_polyrayvec, "finalpolyrayvec.vtk");
	}
}



