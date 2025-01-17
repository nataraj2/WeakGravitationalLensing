#include <WeakLensing.H>     

int main()
{
	WeakLensing wl;

	int nlens = 1;

	// Parameters to create light ray bundles	
	int num_angles_theta = 1000;
	int num_angles_phi = 1000;
	double theta_max = 5e-5;

	//int num_angles_theta = 300;
	//int num_angles_phi = 50;
	//double theta_max = 1e-5;

	int num_rays = num_angles_phi;

	std::vector<double> dist_lens_m, dist_lens_Gpc;
	wl.create_lenses(dist_lens_m, dist_lens_Gpc, nlens);

	// Final polyray bundle which intersects the source
	std::vector<PolyRay> final_polyrayvec;

	for (int itheta=0;itheta<num_angles_theta;itheta++){
		double theta = (itheta+1)*theta_max/num_angles_theta;

		std::vector<PolyRay> polyrayvec;
		polyrayvec.resize(num_rays);

		// Create rays from origin, with the rays being on the surface of a cone and 
		// the for loop increments the half-angle of the cone in each iteration thus covering
		// a large angular diameter

		wl.create_polyray_bundle(theta, polyrayvec);
	
		// Deflect the polyray bundle through all the lenses
		wl.deflect_polyray_bundle(dist_lens_m, polyrayvec);

		wl.intersect_lens_polyray_bundle(dist_lens_Gpc, polyrayvec);

		// Check if the end point of the rays after the deflections
		// intersects the source. If it does, then add it to the 
		// final_polyrayvec

		wl.check_intersect_polyray_bundle(polyrayvec, final_polyrayvec);
		std::string filename;
		if(itheta<=9) {
			filename = "polyray_bundle_00" + std::to_string(itheta) + ".vtk";
		} else if(itheta<=99) {
			filename = "polyray_bundle_0" + std::to_string(itheta) + ".vtk";
		} else if(itheta<=999) {
			filename = "polyray_bundle_" + std::to_string(itheta) + ".vtk";
		}
			
		wl.write_polyray_bundle(polyrayvec,filename);
		
	}

	if(final_polyrayvec.size() > 0) {
		std::cout << "Writing final polyrayvec" << "\n";	
		wl.write_polyray_bundle(final_polyrayvec, "finalpolyrayvec.vtk");
		wl.create_polray_bundle_for_movie(final_polyrayvec);	
	}

	wl.write_lens_and_source_objects();
}



