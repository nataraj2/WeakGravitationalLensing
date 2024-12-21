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

