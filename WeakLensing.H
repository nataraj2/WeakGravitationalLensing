#include <iostream>     // For standard input/output
#include <cmath>    // For mathematical functions like sqrt, pow
#include <iomanip>      // For output formatting
#include <limits>       // For handling numeric limits
#include <vector>
#include <PolyRay.H>

inline double convert_ly_to_m(const double& dist_ly) 
{
    return dist_ly*3e8*3600.0*24.0*365.0;
}

inline double convert_Gpc_to_ly(const double& dist_pc) 
{
    return dist_pc*3.26e9;
}

inline double convert_Gpc_to_m(const double& dist_pc) 
{
    double tmp = convert_Gpc_to_ly(dist_pc);
    tmp = convert_ly_to_m(tmp);
    return tmp;
}

inline double dot_product(const Triplet& v1, const Triplet& v2)
{
	double dp;
	dp = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	return dp;
}


inline Triplet cross_product(const Triplet& v1, const Triplet& v2) 
{
	Triplet tmp(0.0,0.0,0.0);
		
	tmp[0] =  v1[1]*v2[2] - v1[2]*v2[1];
	tmp[1] = -v1[0]*v2[2] + v1[2]*v2[0];
	tmp[2] =  v1[0]*v2[1] - v1[1]*v2[0];

	return tmp;
}

inline Triplet make_unit_vector(const Triplet &v) 
{
	Triplet tmp(0.0,0.0,0.0);
	double mag = std::pow(v[0]*v[0] + v[1]*v[1] + v[2]*v[2],0.5);	
	tmp[0] = v[0]/mag;
	tmp[1] = v[1]/mag;
	tmp[2] = v[2]/mag;
	
	return tmp;
}


class WeakLensing{

public:
    // Constants

    const double G_val = 6.6743e-11;
    const double M_val = 1e12*1.989e30;
    const double c_val = 3e8;
    const double D_OL = 1.0;
    const double D_LS = 0.5;

	void create_lenses(std::vector<double>& dist_lens_m,
                       std::vector<double>& dist_lens_Gpc,
					   const int nlens);

	void create_polyray_bundle(const double& theta,
        	              std::vector<PolyRay>& polyrayvec);

	void deflect_polyray_bundle(const std::vector<double>& dist_lens_m,
                                std::vector<PolyRay>& polyrayvec);

	void write_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
                          const std::string filename);

	void intersect_lens_polyray_bundle(const std::vector<double>& dist_lens_Gpc,
                                  std::vector<PolyRay>& polyrayvec);

	void check_intersect_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
                                        std::vector<PolyRay>& final_polyrayvec);

	void create_polray_bundle_for_movie(const std::vector<PolyRay>& polyrayvec);

	void write_lens_and_source_objects();
};


