#include <WeakLensing.H>

void write_polyray_bundle(const std::vector<PolyRay>& polyrayvec,
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
    std::cout << "Num points is " << num_points <<"\n";
    int total_num_points = num_polyrays*polyrayvec[0].coord.size();
    int total_num_lines = num_polyrays*polyrayvec[0].dir.size();

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
