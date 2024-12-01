/*****************************************************************************
*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <complex>
#include <omp.h>
#include <filesystem>

#include "mptmacros.h"

#include "RandomList.hpp"
#include "Point3D.hpp"
#include "Vector3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"
#include "MRsequence.hpp"
#include "SimSpace.hpp"


int main( ) {

	//==================== Define simulation space (parameters) - SI-Units ==========================//

	double matrix_x0, matrix_y0, matrix_z0; //matrix origin
	double matrix_xlength, matrix_ylength, matrix_zlength; //matrix edge length
	int matrix_nxnodes, matrix_nynodes, matrix_nznodes; //nodes for B-Field calculation
	double radius_ionp, radius_ionp_std; //particle radii + std for lognormal distribution
	double radius_ionp_1, radius_ionp_std_1, percentage_ionp; //second particle radius + std for bimodal lognormal distribution; percentage of particles with first radius
	double radius_agg, radius_agg_std; //radius of spherical aggregate (std for lognormal distribution)
	double thickness_ionp, thickness_ionp_std; //radius of particle with coating (std for lognormal distribution)
	int ionp_number, agg_number, numThreads; //set amount of IONP, Aggregates and Threads for OpenMP
	double magnetic_moment; //Can be defined in ini file or in include/mptmacros.h
	bool DLS_file; //Yes - There is an DLS-Histogram file under /ini/DLS_File/
	Point3D Cube_orig, Cube_size; //Option to define a Cube where the signal is calculated in bigger volume

	//==================== Define diffusion (parameters) ==========================//

	int proton_number, bound_proton_number; //number of simulated protons
	int sequence_num; // Sequence defined in ini
	double echo_time, time_step, echo_spacing;
	double velocityIONP; //options for simulating moving IONP
	Point3D directionIONP;
	double background_T2; //additional background T2 
	double concentration; //in µmol/L (only unit that is not in SI-units)

	//=================== Find configuration files for parameters ======================//

	std::string path_to_ini = "../../../ini";

	for (const auto & entry : std::filesystem::directory_iterator(path_to_ini)){
    std::string path = entry.path();
    if (path.find(".txt") != std::string::npos){

		auto start = std::chrono::system_clock::now();

		std::time_t start_time = std::chrono::system_clock::to_time_t(start);

		Data dat(path);

		//=================== Redirect cout to log file ======================//

		std::string path_to_data = "../../../data/tmp/";
  		std::string path_to_log = "../../../data/log/";
  		int id = 0;
  		int logid = 0;

		for (const auto & entry1 : std::filesystem::directory_iterator(path_to_data)){
    		(void) entry1;
    		id++;
  		}

  		for (const auto & entry2 : std::filesystem::directory_iterator(path_to_log)){
    		(void) entry2;
    		logid++;
  		}

  		std::string log_out_name = path_to_log + std::to_string( logid ) + "_" + std::to_string( id ) + "_log.txt";

  		std::ofstream out(log_out_name);
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  		std::cout.rdbuf(out.rdbuf()); //redirect std::cout to log.txt!

		std::cout << "Started computation at: " << std::ctime(&start_time) << "\n";

		std::cout << path << std::endl;


		//=================== Read ini files and load parameters ======================//

		dat.readData();

		matrix_x0 = dat.get_matrix_x0() ;
		matrix_y0 = dat.get_matrix_y0() ;
		matrix_z0 = dat.get_matrix_z0() ;
		matrix_xlength = dat.get_matrix_xlength() ;
		matrix_ylength = dat.get_matrix_ylength() ;
		matrix_zlength = dat.get_matrix_zlength() ;
		matrix_nxnodes = dat.get_matrix_nxnodes() ;
		matrix_nynodes = dat.get_matrix_nynodes() ;
		matrix_nznodes = dat.get_matrix_nznodes() ;
		radius_ionp = dat.get_radius_ionp() ;
		radius_ionp_1 = dat.get_radius_ionp_1() ;
		radius_agg = dat.get_radius_agg() ;
		radius_ionp_std = dat.get_radius_ionp_std() ;
		radius_ionp_std_1 = dat.get_radius_ionp_std_1() ;
		radius_agg_std = dat.get_radius_agg_std() ;
		ionp_number = dat.get_ionp_number() ;
		agg_number = dat.get_agg_number() ;
		proton_number = dat.get_proton_number();
		echo_time = dat.get_echo_time() ;
		time_step = dat.get_time_step() ;
		echo_spacing = dat.get_echo_spacing() ;
		numThreads = dat.get_numThreads() ;
		sequence_num = dat.get_sequence_num() ;
		velocityIONP = dat.get_velocity_IONP() ;
		directionIONP = dat.get_direction_IONP() ;
		thickness_ionp = dat.get_thickness_IONPcoating() ;
		thickness_ionp_std = dat.get_thickness_IONPcoating_std() ;
		bound_proton_number = dat.get_bound_proton_number() ;
		background_T2 = dat.get_background_T2() ;
		concentration = dat.get_concentration() ;
		magnetic_moment = dat.get_mag_moment() ;
		percentage_ionp = dat.get_ionp_percentage() ;
		Cube_orig = dat.get_Cube_Origin() ;
		Cube_size = dat.get_Cube_size() ;
		DLS_file = dat.get_DLS_file() ;

		//=================== Calculate number of simulation steps either given time steps or from IONP size ======================//

		int steps_number;

		if (time_step == 0){
			if ( thickness_ionp == 0 ){
				steps_number = SC_I(round(echo_time / (square(radius_ionp) / (6 * D_D_CONST)))); // number of steps
			}
			else{
				steps_number = SC_I(round(echo_time / (square(thickness_ionp) / (6 * D_D_CONST)))); // number of steps
			}
			std::cout << "Steps_n: " << steps_number <<"\n" << std::endl;
		}
		else{
			steps_number = SC_I(round(echo_time / time_step)); // number of steps
			std::cout << "Steps_n: " << steps_number <<"\n" << std::endl;
		}

		//=================== Initialization of simulation matrix (SimMatrix.cpp) ======================//

		SimMatrix sMatrix = SimMatrix(matrix_x0, matrix_y0, matrix_z0, matrix_xlength, matrix_ylength, matrix_zlength, matrix_nxnodes, matrix_nynodes, matrix_nznodes);

		//Already define vectors that will contain IONP/Aggregate informations (optional/not fixed)

		std::vector<Ionp> ionps;
		ionps.reserve(SC_UI(ionp_number));

		std::vector<Ionp> ionpMove;
		ionpMove.reserve(SC_UI(ionp_number));

		std::vector<Aggregate> aggregates;
		aggregates.reserve(SC_UI(agg_number));

		//=================== Placement and definition of Properties of IONP/Aggregates (SimSpace.cpp) ======================//

		if (ionp_number == 0){
			//Option to manually define IONP position a position file at /ini/IONP_pos/ is needed
			if (matrix_xlength == 0.){
				std::cerr << "Size of Simulation volume has to be given, if IONPs are placed by hand!" << std::endl;
				exit(1);
			}

			std::string path_to_ionp = "../../../ini/IONP_pos";

			for (const auto & entryIONP : std::filesystem::directory_iterator(path_to_ionp)){
    		std::string pathIONP = entryIONP.path();
    		if (pathIONP.find(".txt") != std::string::npos || pathIONP.find(".dat") != std::string::npos ){
				dat.readIONP( pathIONP, ionps, sMatrix, radius_ionp);
			}
			else{
				std::cerr << "No IONP number given and no IONP position file found!" << std::endl;
				exit(1);
			}
			}

		}
		else{
			if ( agg_number == 0 && radius_agg == 0){
				SimSpace::simulationSpace(ionps, sMatrix, ionp_number, concentration, radius_ionp, radius_ionp_1, radius_ionp_std, radius_ionp_std_1, thickness_ionp, thickness_ionp_std, percentage_ionp, DLS_file);
			}
			else if ( agg_number != 0 && radius_agg == 0 ){
				SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, concentration, radius_ionp, radius_ionp_1, radius_ionp_std, radius_ionp_std_1, thickness_ionp, percentage_ionp);
			}
			else if ( agg_number == 0 && radius_agg != 0 ){
				std::string path_to_agg = "../../../ini/Agg_pos";

				for (const auto & entryAgg : std::filesystem::directory_iterator(path_to_agg)){
    				std::string pathAgg = entryAgg.path();
    				if (pathAgg.find(".txt") != std::string::npos || pathAgg.find(".dat") != std::string::npos ){
						dat.readAgg( pathAgg, aggregates, sMatrix, radius_agg);
					}
					else{
						std::cerr << "No Aggregates number given and no Aggregates position file found!" << std::endl;
						exit(1);
					}
				}

				SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, concentration, radius_agg, radius_ionp, radius_ionp_1, radius_ionp_std, radius_ionp_std_1, radius_agg_std, thickness_ionp, thickness_ionp_std, percentage_ionp);
			}
			else{
				SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, concentration, radius_agg, radius_ionp, radius_ionp_1, radius_ionp_std, radius_ionp_std_1, radius_agg_std, thickness_ionp, thickness_ionp_std, percentage_ionp);
			}
		}

		//Already define vectors that will contain Proton/Time/MRISignal informations (optional/not fixed)

		std::vector<Proton> protons;
		protons.reserve(SC_UI(agg_number * ionp_number));

		std::vector<Proton> protonPos;
		protonPos.reserve(SC_UI(agg_number * ionp_number));

		std::vector<std::complex<double>> signalV;
		signalV.reserve(SC_UI(steps_number + 1));

		std::vector<std::complex<double>> signalV_MSE;
		signalV_MSE.reserve(SC_UI(steps_number + 1));

		std::vector<double> timeV;
		timeV.reserve(SC_UI(steps_number + 1));

		std::vector<double> timeV_MSE;
		timeV_MSE.reserve(SC_UI(steps_number + 1));

		std::vector<double> phaseV;
		phaseV.reserve(SC_UI(steps_number + 1));

		//Define if variable Time steps and a B-Field grid is used for Simulation depending on Parameters given

		bool Grid, var_Steps = false;

		//No number of nodes defined = no B-Field grid
		if (matrix_nxnodes == 0 && matrix_nynodes == 0 && matrix_nznodes == 0){
			//No time steps = variable time steps used (smaller time steps if protons are closer than R*8 to IONP)
			if ( time_step == 0. ) {
				std::cout << "Variable diffusion step length (R or R/8) used" << std::endl;
				if ( thickness_ionp == 0 ){
					time_step = square(radius_ionp) / (6 * D_D_CONST);
				}
				else{
					time_step = square(thickness_ionp) / (6 * D_D_CONST);
				}
				dat.set_time_step(time_step);
				Grid = false;
				var_Steps = true;
			}
			else{
				Grid = false;
				var_Steps = false;
			}
		}
		else {
			if ( time_step == 0. ) {
				std::cout << "Variable diffusion step length (R or R/8) AND B-Field Grid with interpolation used" << std::endl;
				if ( thickness_ionp == 0 ){
					time_step = square(radius_ionp) / (6 * D_D_CONST);
				}
				else{
					time_step = square(thickness_ionp) / (6 * D_D_CONST);
				}
				dat.set_time_step(time_step);
				Grid = true;
				var_Steps = true;
			}
			else{
				std::cout << "B-Field Grid with interpolation used" << std::endl;
				Grid = true;
				var_Steps = false;
			}
		}

		//=================== Simulation of Protons + Calculation of Phase and MR Signal (MRsequence.cpp) ======================//

		MRsequence MRI(sMatrix, sequence_num, proton_number, bound_proton_number, radius_ionp, echo_time, time_step, echo_spacing, numThreads, matrix_nxnodes, matrix_xlength, directionIONP, velocityIONP, thickness_ionp, background_T2, magnetic_moment, Cube_orig, Cube_size, Grid, var_Steps);

		MRI.Sequence(protons, ionps, signalV, timeV, signalV_MSE, timeV_MSE);


		//============ Outfile magnetization =================================//

		dat.writeData(signalV, timeV);

		if (sequence_num == 2){
			dat.writeData(signalV_MSE, timeV_MSE);
		}

		dat.writeData(ionps);

		if ( agg_number != 0 || radius_agg != 0 ){
			dat.writeData(aggregates);
		}

		auto end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		std::cout << "finished computation at " << std::ctime(&end_time)
			<< "elapsed time: " << elapsed_seconds.count() << "s\n";


		std::cout.rdbuf(coutbuf); //reset to standard output again

    }
	}

	return 0;
}
