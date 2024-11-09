#include <WeakLensing.H>

void
WeakLensing::create_lenses(std::vector<double>& dist_lens_m,
                           std::vector<double>& dist_lens_Gpc,
						   const int nlens)
{
    dist_lens_m.resize(nlens+1,0.0);
    dist_lens_Gpc.resize(nlens+1,0.0);
	
    dist_lens_m[0]= convert_Gpc_to_m(D_OL);
    dist_lens_m[1]= convert_Gpc_to_m(D_LS);

    dist_lens_Gpc[0]= D_OL;
    dist_lens_Gpc[1]= D_LS;

	double D_OS = D_OL + D_LS;

    double D_OL_m = convert_Gpc_to_m(D_OL);
    double D_LS_m = convert_Gpc_to_m(D_LS);
    double D_OS_m = convert_Gpc_to_m(D_OS);

    double theta_e = std::pow(4.0*G_val*M_val/(c_val*c_val)*D_LS_m/(D_OL_m*D_OS_m),0.5);
    std::cout << "Value of theta_e is " << theta_e << " radians = " << theta_e*180.0/M_PI*3600.0 << " arcseconds " << "\n";
}

void 
WeakLensing::create_polyray_bundle(const double& theta,
					      		   std::vector<PolyRay>& polyrayvec)
{

	int num_rays = polyrayvec.size();

	double dir[3];

	int polyray_count = 0;
	for (int iphi=0;iphi<num_rays;iphi++) {
    	polyrayvec[polyray_count].theta_start = theta;
        double phi = iphi*2.0*M_PI/(num_rays-1.0);
        polyrayvec[polyray_count].phi_start = phi;
        dir[0] = std::sin(theta)*std::cos(phi);
        dir[1] = std::sin(theta)*std::sin(phi);
        dir[2] = std::cos(theta);

        polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
        polyray_count++;
     }

	 // All rays start at the observer
    for (int i=0;i<num_rays;i++) {
        polyrayvec[i].coord.emplace_back(0.0, 0.0, 0.0);
    }
}

void 
WeakLensing::deflect_polyray_bundle(const std::vector<double>& dist_lens_m,
									std::vector<PolyRay>& polyrayvec)
{
	int num_rays = polyrayvec.size();

	double dir[3];
// Deflection caused by the lens is given by
// alpha = 4*G*M/(c*c*xi) where xi is the distance of the ray on
// the lens plane

	 for (int iray=0;iray<num_rays;iray++){
        // Find the intersection on the lens plane
        double x = polyrayvec[iray].coord[0][0] + dist_lens_m[0]*polyrayvec[iray].dir[0][0];
        double y = polyrayvec[iray].coord[0][1] + dist_lens_m[0]*polyrayvec[iray].dir[0][1];

        double xi = std::pow(x*x+y*y,0.5);
        double alpha = 4.0*G_val*M_val/(c_val*c_val*xi);

        double theta = polyrayvec[iray].theta_start - alpha;
        double phi = polyrayvec[iray].phi_start;
        dir[0] = std::sin(theta)*std::cos(phi);
        dir[1] = std::sin(theta)*std::sin(phi);
        dir[2] = std::cos(theta);

        polyrayvec[iray].dir.emplace_back(dir[0], dir[1], dir[2]);
    }
}

void 
WeakLensing::intersect_lens_polyray_bundle(const std::vector<double>& dist_lens_Gpc,
										   std::vector<PolyRay>& polyrayvec)
{

	int num_rays = polyrayvec.size();
	int nlens = polyrayvec[0].dir.size()-1;

 	// Find the end of each ray section of the poly ray
    // The ray direction remains the same between lenses
    // Hence if there are n lenses, there will be n+1 coordinates
    // for eg. there is 1 lens there will be 2 coordinates

    for(int iray=0;iray<num_rays;iray++){
        for(int i=1;i<=nlens+1;i++) {
            double x = polyrayvec[iray].coord[i-1][0] + dist_lens_Gpc[i-1]*polyrayvec[iray].dir[i-1][0];
            double y = polyrayvec[iray].coord[i-1][1] + dist_lens_Gpc[i-1]*polyrayvec[iray].dir[i-1][1];
            double z = polyrayvec[iray].coord[i-1][2] + dist_lens_Gpc[i-1]*polyrayvec[iray].dir[i-1][2];
            polyrayvec[iray].coord.emplace_back(x,y,z);
        }
    }
}

void 
WeakLensing::check_intersect_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
											std::vector<PolyRay>& final_polyrayvec)
{
	int num_rays = polyrayvec.size();
	// Do plotting in Gpc units
    // Does the end of the polyray intersect the source
	int idx = polyrayvec[0].coord.size()-1;	

    for(int iray=0;iray<num_rays;iray++){
        double x = polyrayvec[iray].coord[idx][0];
        double y = polyrayvec[iray].coord[idx][1];
        //double z = polyrayvec[iray].coord[nlens+1][2];
        //std::cout << "Value of " << iray << " " << std::pow(x*x + y*y,0.5) << "\n";
        if(std::pow(x*x + y*y,0.5) < 1e-7){
            //std::cout << "Inside here " << std::pow(x*x + y*y,0.5) << "\n";
            //exit(0);
            final_polyrayvec.push_back(polyrayvec[iray]);
        }
    }
}


