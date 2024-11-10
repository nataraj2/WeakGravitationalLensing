#include <WeakLensing.H>

void
WeakLensing::write_lens_and_source_objects()
{
	FILE *file_lensplane;
    file_lensplane = fopen("lens.vtk","w");
    fprintf(file_lensplane, "%s\n","# vtk DataFile Version 3.0");
    fprintf(file_lensplane, "%s\n","Lens plane");
    fprintf(file_lensplane, "%s\n","ASCII");
    fprintf(file_lensplane, "%s\n","DATASET STRUCTURED_POINTS");
	fprintf(file_lensplane, "%s %d %d %d\n", "DIMENSIONS", 2, 2, 1);
	fprintf(file_lensplane, "%s %0.15g %0.15g %0.15g\n","ORIGIN", -1e-5, -1e-5, D_OL*1e-4);
	fprintf(file_lensplane, "%s %0.15g %0.15g %0.15g\n","SPACING", 2e-5, 2e-5, D_OL*1e-4);
	fprintf(file_lensplane, "%s %d\n","POINT_DATA", 4);
	fprintf(file_lensplane, "%s %d\n", "SCALARS Temperature float", 1);
	fprintf(file_lensplane, "%s\n", "LOOKUP_TABLE default");
	fprintf(file_lensplane, "%0.15g %0.15g\n", 1.0,1.0);
	fprintf(file_lensplane, "%0.15g %0.15g\n", 1.0,1.0);
	fclose(file_lensplane); 


	FILE *file_objects;
    file_objects = fopen("source_and_lens.vtk","w");
    fprintf(file_objects, "%s\n","# vtk DataFile Version 3.0");
    fprintf(file_objects, "%s\n","Lens and source with points");
    fprintf(file_objects, "%s\n","ASCII");
    fprintf(file_objects, "%s\n","DATASET POLYDATA");
	
	fprintf(file_objects, "%s %d %s\n", "POINTS", 2, "float");
	fprintf(file_objects, "%0.15g %0.15g %0.15g\n",0.0,0.0,D_OL*1e-4); 
	fprintf(file_objects, "%0.15g %0.15g %0.15g\n",0.0,0.0,(D_OL+D_LS)*1e-4);
	fclose(file_objects);

}

void 
WeakLensing::write_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
                          		  const std::string filename)
{
    FILE *file_polyrays;
    file_polyrays = fopen(filename.c_str(),"w");
    fprintf(file_polyrays, "%s\n","# vtk DataFile Version 3.0");
    fprintf(file_polyrays, "%s\n","Poly ray bundle");
    fprintf(file_polyrays, "%s\n","ASCII");
    fprintf(file_polyrays, "%s\n","DATASET POLYDATA");


    int num_polyrays = polyrayvec.size();
    int num_points = polyrayvec[0].coord.size();
    int total_num_points = num_polyrays*polyrayvec[0].coord.size();
    int total_num_lines = num_polyrays*(polyrayvec[0].coord.size()-1);

    fprintf(file_polyrays, "%s %d %s\n", "POINTS", total_num_points, "float");

    for(int iray=0;iray<num_polyrays;iray++) {
        for(int i=0;i<num_points;i++){
            fprintf(file_polyrays, "%0.15g %0.15g %0.15g\n", polyrayvec[iray].coord[i][0],
                                                             polyrayvec[iray].coord[i][1],
                                                             polyrayvec[iray].coord[i][2]*1e-4);
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

void 
WeakLensing::create_polray_bundle_for_movie(const std::vector<PolyRay>& polyrayvec)
{

	// Make 10 points on each segment of the rays

    int num_rays = polyrayvec.size();
    int num_points = polyrayvec[0].coord.size();

	int np_seg = 20;
	double tmp_coord[3];
	double min_coord[3], max_coord[3];

	std::vector<PolyRay> tmp_polyrayvec;
    tmp_polyrayvec.resize(num_rays);

	for(int iray=0;iray<num_rays;iray++) {
		for(int i=0;i<num_points-1;i++){
			for(int ii=0; ii<3; ii++){
            	min_coord[ii] = polyrayvec[iray].coord[i][ii];
                max_coord[ii] = polyrayvec[iray].coord[i+1][ii];
            }
			for(int p=0;p<np_seg;p++) {	
				for(int ii=0; ii<3; ii++){
            		double del_coord = (max_coord[ii]-min_coord[ii])/(np_seg-1.0);
            		tmp_coord[ii] = min_coord[ii] + del_coord*p;
            	}
			 	tmp_polyrayvec[iray].coord.emplace_back(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
			}
		}
	}	
	
	num_points = tmp_polyrayvec[0].coord.size();
	for(int i=1;i<num_points;i++){
		std::vector<PolyRay> toplot_polyrayvec;
		toplot_polyrayvec.resize(num_rays);

	 	std::string filename;
		int index = i;
       	if(index<=9){
       		filename = "tmp_polyray_bundle_00" + std::to_string(index) + ".vtk";
       	} else if(index<=99) {
        	filename = "tmp_polyray_bundle_0" + std::to_string(index) + ".vtk";
       	} else if (index<=999) {
      		filename = "tmp_polyray_bundle_" + std::to_string(index) + ".vtk";
        }	
		for(int iray=0;iray<num_rays;iray++){
			for(int p=0;p<=i;p++){
				double x = tmp_polyrayvec[iray].coord[p][0];
				double y = tmp_polyrayvec[iray].coord[p][1];
				double z = tmp_polyrayvec[iray].coord[p][2];
				
				toplot_polyrayvec[iray].coord.emplace_back(x, y, z);	
			}
		}
        write_polyray_bundle(toplot_polyrayvec, filename);
	}
}
