#include <iostream>     // For standard input/output
#include <cmath>	// For mathematical functions like sqrt, pow
#include <iomanip>      // For output formatting
#include <limits>       // For handling numeric limits
#include <vector>


class Triplet
{
	public:
		Triplet(double val0, double val1, double val2) : 
				m_val0(val0), m_val1(val1), m_val2(val2) {} 

	// Overload operator[] to provide array-like access to x, y, z
    double& operator[](size_t index) {
        if (index == 0) return m_val0;
        if (index == 1) return m_val1;
        if (index == 2) return m_val2;
        throw std::out_of_range("Index out of range for Triplet");
    }

    const double& operator[](size_t index) const {
        if (index == 0) return m_val0;
        if (index == 1) return m_val1;
        if (index == 2) return m_val2;
        throw std::out_of_range("Index out of range for Triplet");
    }

	private:
		double m_val0, m_val1, m_val2;
};

class PolyRay 
{
	public:
		std::vector<Triplet> coord;
		std::vector<Triplet> dir;
};

void write_polyray_bundle(const std::vector<PolyRay>& polyrayvec)
{
	FILE *file_polyrays;
	file_polyrays = fopen("polyray_bundle.vtk","w");
	fprintf(file_polyrays, "%s\n","# vtk DataFile Version 3.0");
	fprintf(file_polyrays, "%s\n","Poly ray bundle");
	fprintf(file_polyrays, "%s\n","ASCII");
	fprintf(file_polyrays, "%s\n","DATASET POLYDATA");

		
	int num_polyrays = polyrayvec.size();
	int num_points = polyrayvec[0].coord.size();
	std::cout << "Num points is " << num_points <<"\n";
	int total_num_points = num_polyrays*polyrayvec[0].coord.size();
	int total_num_lines = num_polyrays*polyrayvec[0].dir.size();

	fprintf(file_polyrays, "%s %d %s\n", "POINTS", total_num_points, "float");

	for(int iray=0;iray<num_polyrays;iray++) {
		for(int i=0;i<num_points;i++){
			fprintf(file_polyrays, "%0.15g %0.15g %0.15g\n", polyrayvec[iray].coord[i][0],
															 polyrayvec[iray].coord[i][1],
															 polyrayvec[iray].coord[i][2]);	
		}
	}

	fprintf(file_polyrays, "%s %d %d\n", "LINES", total_num_lines, total_num_lines*3);
	for(int iray=0;iray<num_polyrays;iray++) {
		for(int i=0;i<num_points-1;i++){
			fprintf(file_polyrays, "%d %d %d\n", 2, iray*num_points + i, iray*num_points + i + 1);
		}
	}
	fclose(file_polyrays);

}


int main()
{

	// Create rays from origin, with the rays being on the surface of a cone and 
	// the for loop increments the half-angle of the cone in each iteration thus covering
	// a large angular diameter

	// Observer is at the origin

	int num_angles_theta = 1;
	int num_angles_phi = 20;
	int nlens = 4;

	int num_rays = num_angles_theta*num_angles_phi;
	std::vector<PolyRay> polyrayvec;
	polyrayvec.resize(num_rays);
	
	double xloc=0.0, yloc=0.0, zloc=0.0;

	double dir[3];

	double theta_max = M_PI/20.0;
	
	// Vary azimuthal angle phi 
	
	int polyray_count = 0;
	
	for (int itheta=0;itheta<num_angles_theta;itheta++){
		double theta = (itheta+1)*theta_max/num_angles_theta;
		for (int iphi=0;iphi<num_angles_phi;iphi++) {
			double phi = iphi*2.0*M_PI/(num_angles_phi-1.0);
			dir[0] = std::sin(theta)*std::cos(phi);
			dir[1] = std::sin(theta)*std::sin(phi);
			dir[2] = std::cos(theta);

			polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
			std::cout << "Dirs are " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
			polyray_count++;
		}
	}

	polyray_count = 0;
	for (int itheta=0;itheta<num_angles_theta;itheta++){
        double theta = (itheta+1)*theta_max/num_angles_theta - M_PI/60.0;
        for (int iphi=0;iphi<num_angles_phi;iphi++) {
            double phi = iphi*2.0*M_PI/(num_angles_phi-1.0);
            dir[0] = std::sin(theta)*std::cos(phi);
            dir[1] = std::sin(theta)*std::sin(phi);
            dir[2] = std::cos(theta);

            polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
            std::cout << "Dirs are " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
            polyray_count++;
        }
    }


	polyray_count = 0;
	for (int itheta=0;itheta<num_angles_theta;itheta++){
        double theta = (itheta+1)*theta_max/num_angles_theta - M_PI/30.0;
        for (int iphi=0;iphi<num_angles_phi;iphi++) {
            double phi = iphi*2.0*M_PI/(num_angles_phi-1.0);
            dir[0] = std::sin(theta)*std::cos(phi);
            dir[1] = std::sin(theta)*std::sin(phi);
            dir[2] = std::cos(theta);

            polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
            std::cout << "Dirs are " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
            polyray_count++;
        }
    }

	polyray_count = 0;
	for (int itheta=0;itheta<num_angles_theta;itheta++){
        double theta = (itheta+1)*theta_max/num_angles_theta - M_PI/10.0;
        for (int iphi=0;iphi<num_angles_phi;iphi++) {
            double phi = iphi*2.0*M_PI/(num_angles_phi-1.0);
            dir[0] = std::sin(theta)*std::cos(phi);
            dir[1] = std::sin(theta)*std::sin(phi);
            dir[2] = std::cos(theta);

            polyrayvec[polyray_count].dir.emplace_back(dir[0], dir[1], dir[2]);
            std::cout << "Dirs are " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
            polyray_count++;
        }
    }





	
	// z is the direction of shoot axis
	// Shoot till z = z_lens

	std::vector<double> dist_lens;
	dist_lens.resize(nlens,0.0);
	dist_lens[0]= 100.0;
	dist_lens[1]= 100.0;
	dist_lens[2]= 100.0;
	dist_lens[3]= 100.0;
	
	// All rays start at the observer	
	for (int i=0;i<num_rays;i++) {
		polyrayvec[i].coord.emplace_back(xloc, yloc, zloc);
	}

	// Find the end of each ray section of the poly ray
	// The ray direction remains the same between lenses
	// Hence if there are n lenses, there will be n+1 coordinates
	// for eg. there is 1 lens there will be 2 coordinates

	for(int iray=0;iray<num_rays;iray++){
		for(int i=1;i<=nlens;i++) {	
			double x = polyrayvec[iray].coord[i-1][0] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][0];
			double y = polyrayvec[iray].coord[i-1][1] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][1];
			double z = polyrayvec[iray].coord[i-1][2] + dist_lens[i-1]*polyrayvec[iray].dir[i-1][2];
			std::cout << "Values are x, y, z " << x << " " << y << " " << z << "\n";
			polyrayvec[iray].coord.emplace_back(x,y,z);	
	 	}
	}

	write_polyray_bundle(polyrayvec);
	
}



