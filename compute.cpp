// compute voronoi cell of non-rectangular box
//
// Author :: wei xuerui
// Email  :: wxrui2013@gmail.com
// date   :: June 10th 2017

#include "voro++.hh"
using namespace voro;

#include <iostream>
  
using namespace std;  
  
// Set up constants for the container geometry
double tmp_bx,tmp_by,tmp_bz,tmp_bxy,tmp_bxz,tmp_byz;
int 	filenumber;
//Set up the number of blocks that the container is divided into
const int n_x=15,n_y=15,n_z=15;

int main(int argc,char **argv) {
	// reading parameter of shape of sample 
	// parallilepipe  has six argument : bx,by,bz,xy,xz,yz.
	tmp_bx=strtod( argv[1] ),tmp_by=strtod(argv[2]),tmp_bz=strtod(argv[3]),
	tmp_bxy=strtod(argv[4]),tmp_bxy=strtod(argv[5]),tmp_byz=strtod(argv[6]),
	filenumber=atoi(argv[7]);

using namespace voro;
	//create a container with the geometry given above
	container_periodic_poly con(tmp_bx,tmp_bxy,tmp_by,tmp_bxz,
				tmp_byz,tmp_bz,n_x,n_y,n_z,32);
	
	//put particals into the container
	con.import("data1.dat");
	
	//Save the voronoi network of all the particals to text files in 
	//gnuplot formats
	con.draw_cells_gnuplot("particals_pack.gnu");
	
	//compute voronoi cell of all particals in container
	con.compute_all_cells();
	
	//print information about voronoi cell(ID , number of faces, 
	//number of variety face (3,4,5,6,...))
	con.print_custom("%i %s %A  @%i %v @%i %s %n ","packing.dat");
	
}
