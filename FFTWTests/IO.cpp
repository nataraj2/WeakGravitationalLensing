#include <DarkMatter.H>

void
plot_DM_particles_vtk (std::vector<DMParticle>& dm_particles)
{
	FILE* file_turbloc_vtk;
	file_turbloc_vtk = fopen("DM_particles.vtk","w");
	fprintf(file_turbloc_vtk, "%s\n","# vtk DataFile Version 3.0");
	fprintf(file_turbloc_vtk, "%s\n","Dark matter particles");
	fprintf(file_turbloc_vtk, "%s\n","ASCII");
	fprintf(file_turbloc_vtk, "%s\n","DATASET POLYDATA");
	fprintf(file_turbloc_vtk, "%s %ld %s\n", "POINTS", dm_particles.size(), "float");
	for(long unsigned int it=0; it < dm_particles.size(); it++){
	    fprintf(file_turbloc_vtk, "%0.15g %0.15g %0.15g\n", dm_particles[it].x, dm_particles[it].y, dm_particles[it].z);
	}
	fclose(file_turbloc_vtk);
}


// Function to write the scalar field to a VTK file
void writeStructuredGridVTK(const std::string& filename, int nx, int ny, int nz, const std::vector<double>& scalarField)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write VTK file header
    file << "# vtk DataFile Version 3.0\n";
    file << "Structured Grid Example\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";

    // Write grid points
    file << "POINTS " << nx * ny * nz << " float\n";
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                file << i << " " << j << " " << k << "\n";
            }
        }
    }

    // Write scalar data
    file << "POINT_DATA " << nx * ny * nz << "\n";
    file << "SCALARS scalar_field float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const double& value : scalarField) {
        file << value << "\n";
    }

    file.close();
    std::cout << "VTK file written to " << filename << "\n";
}


