/*****************************************************************************
*Class to define simulation matrix (cpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <complex>
#include <cmath>
#include <omp.h>
#include <filesystem>

#include "mptmacros.h"

#include "Point3D.hpp"
#include "Vector3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"


double SimMatrix::calculateBfield(const Point3D& p1, const Point3D& p2, const double r) {
	double distance_p1_p2 = distance(p1, p2);
	double cos_theta = (p2.z() - p1.z()) / distance_p1_p2;
	double B_angled_part = cube(r/distance_p1_p2)*(3 * square(cos_theta) - 1);

	return B_fixed * B_angled_part;

}

double SimMatrix::calculateBfield(const Point3D& p1, const Point3D& p2, const double r, const double magM) {
	double distance_p1_p2 = distance(p1, p2);
	double cos_theta = (p2.z() - p1.z()) / distance_p1_p2;
	(void) r;
	double B_angled_part = cube(1./distance_p1_p2)*(3 * square(cos_theta) - 1);
	double B_fixed_new = D_MO/(4.*D_PI) * magM;

	return B_fixed_new * B_angled_part;

}


void SimMatrix::calculateBfield(std::vector< std::vector< std::vector<double> > >& BfieldGrid, std::vector<Ionp>& ionpInstances, int numThreads, double magnetic_moment) {

	std::cout << "Starting to calculate the B-Field Grid...." << std::endl;

	omp_set_num_threads( numThreads );

	double half_xlength = -m_x_length/2.;
	double half_ylength = -m_y_length/2.;
	double half_zlength = -m_z_length/2.;

	double stepsize_x = m_x_length/(m_nx_nodes-1);
	double stepsize_y = m_y_length/(m_ny_nodes-1);
	double stepsize_z = m_z_length/(m_nz_nodes-1);

	#pragma omp parallel for
	for (unsigned int z = 0; z < SC_UI(m_nz_nodes); ++z) {
    	for (unsigned int x = 0; x < SC_UI(m_nx_nodes); ++x) {
      		for (unsigned int y = 0; y < SC_UI(m_ny_nodes); ++y) {

				Point3D gridPoint;
				gridPoint.set((half_xlength + stepsize_x * x) ,(half_ylength + stepsize_y * y) , (half_zlength + stepsize_z * z) );

				double bfield_tot = 0.;
				if (magnetic_moment != 0.){
					for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
						bfield_tot += calculateBfield(gridPoint, ionpInstances[l].xt(), ionpInstances[l].radius(), magnetic_moment);
					}
				}
				else{
					for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
						bfield_tot += calculateBfield(gridPoint, ionpInstances[l].xt(), ionpInstances[l].radius());
					}
				}

				BfieldGrid[x][y][z] = bfield_tot;
      		}
    	}
  	}


	std::cout << "B-Field Grid calculation finished...." << std::endl;

}

//Trilinear interpolation
double SimMatrix::interpolate(const Point3D& p1, std::vector< std::vector< std::vector<double> > >& BfieldGrid){

	unsigned int x_indexUP = SC_UI(ceil((p1.x() / (m_x_length/(m_nx_nodes-1))) + (m_nx_nodes-1)/2.));
	unsigned int x_indexDOWN = SC_UI(floor((p1.x() / (m_x_length/(m_nx_nodes-1))) + (m_nx_nodes-1)/2.));

	unsigned int y_indexUP = SC_UI(ceil((p1.y() / (m_y_length/(m_ny_nodes-1))) + (m_ny_nodes-1)/2.));
	unsigned int y_indexDOWN = SC_UI(floor((p1.y() / (m_y_length/(m_ny_nodes-1))) + (m_ny_nodes-1)/2.));

	unsigned int z_indexUP = SC_UI(ceil((p1.z() / (m_z_length/(m_nz_nodes-1))) + (m_nz_nodes-1)/2.));
	unsigned int z_indexDOWN = SC_UI(floor((p1.z() / (m_z_length/(m_nz_nodes-1))) + (m_nz_nodes-1)/2.));

	double c000, c010, c110, c100, c001, c011, c111, c101;
	double xd, yd, zd;

	c000 = BfieldGrid[x_indexDOWN][y_indexDOWN][z_indexDOWN];
	c010 = BfieldGrid[x_indexDOWN][y_indexUP][z_indexDOWN];
	c110 = BfieldGrid[x_indexUP][y_indexUP][z_indexDOWN];
	c100 = BfieldGrid[x_indexUP][y_indexDOWN][z_indexDOWN];

	c001 = BfieldGrid[x_indexDOWN][y_indexDOWN][z_indexUP];
	c011 = BfieldGrid[x_indexDOWN][y_indexUP][z_indexUP];
	c111 = BfieldGrid[x_indexUP][y_indexUP][z_indexUP];
	c101 = BfieldGrid[x_indexUP][y_indexDOWN][z_indexUP];

	xd = (p1.x() - (-m_x_length/2. + (m_x_length/(m_nx_nodes-1))*x_indexDOWN)) / ((-m_x_length/2. + (m_x_length/(m_nx_nodes-1))*x_indexUP) - (-m_x_length/2. + (m_x_length/(m_nx_nodes-1))*x_indexDOWN));
	yd = (p1.y() - (-m_y_length/2. + (m_y_length/(m_ny_nodes-1))*y_indexDOWN)) / ((-m_y_length/2. + (m_y_length/(m_ny_nodes-1))*y_indexUP) - (-m_y_length/2. + (m_y_length/(m_ny_nodes-1))*y_indexDOWN));
	zd = (p1.z() - (-m_z_length/2. + (m_z_length/(m_nz_nodes-1))*z_indexDOWN)) / ((-m_z_length/2. + (m_z_length/(m_nz_nodes-1))*z_indexUP) - (-m_z_length/2. + (m_z_length/(m_nz_nodes-1))*z_indexDOWN));

	double c00 = c000*(1-xd) + c100*xd;
	double c01 = c001*(1-xd) + c101*xd;
	double c10 = c010*(1-xd) + c110*xd;
	double c11 = c011*(1-xd) + c111*xd;

	double c0 = c00*(1-yd) + c10*yd;
	double c1 = c01*(1-yd) + c11*yd;

	double c = c0*(1-zd) + c1*zd;


	return c;
}
