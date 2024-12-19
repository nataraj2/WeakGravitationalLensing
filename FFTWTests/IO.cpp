#include <DarkMatter.H>

void
plot_DM_particles_vtk ()
{
	FILE* file_turbloc_vtk;
	file_turbloc_vtk = fopen("DM_particles.vtk","w");
	fprintf(file_turbloc_vtk, "%s\n","# vtk DataFile Version 3.0");
	fprintf(file_turbloc_vtk, "%s\n","Dark matter particles");
	fprintf(file_turbloc_vtk, "%s\n","ASCII");
	fprintf(file_turbloc_vtk, "%s\n","DATASET POLYDATA");
	/*fprintf(file_turbloc_vtk, "%s %ld %s\n", "POINTS", xloc.size(), "float");
	for(int it=0; it<xloc.size(); it++){
	    fprintf(file_turbloc_vtk, "%0.15g %0.15g %0.15g\n", xloc[it], yloc[it], hub_height);
	}*/
	fclose(file_turbloc_vtk);
}

