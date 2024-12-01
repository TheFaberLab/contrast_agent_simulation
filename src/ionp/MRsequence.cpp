/*****************************************************************************
*Class to apply different MR sequences

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


//Function to perform proton diffusion and check new position
void MRsequence::replaceProton(std::vector<Proton>& protonInstances, std::vector<Ionp>& ionpInstances, double dt, uint i, uint j, RandomList& RandomNormalX, RandomList& RandomNormalY, RandomList& RandomNormalZ, RandomList& RandomDiff){

	/*At the moment:
	outside while loop to check if proton diffuses into IONP; if yes repeat diffusion step until proton and IONP don't overlap; this is checked last
	1) Check if proton diffuses outside of Voxel -> periodocity
	2) Check if proton ends up inside IONP -> repeat diffusion step
	*/

	Proton currentProtonInstance = Proton(0., 0.);
	while (true) {
		bool repeat = false;
		double extra_length = 0.;
		Point3D lastPosition;
		if (!j) {

			protonInstances[i].placeRandom(RandomNormalX.getNumber(), RandomNormalY.getNumber(), RandomNormalZ.getNumber(), true);

			lastPosition = protonInstances[i].xt();

			currentProtonInstance = protonInstances[i];
		}
		else {
			lastPosition = protonInstances[i].xt();
			currentProtonInstance = protonInstances[i];
			currentProtonInstance.diffuse(D_D_CONST, dt, RandomDiff.getNumber(), RandomDiff.getNumber(), RandomDiff.getNumber());
		}

		if (abs(currentProtonInstance.xt().x()) > simMatrix.xLength() / 2.) { //check if proton outside simulation volume -> periodocity (placed opposite side)
			extra_length = abs(currentProtonInstance.xt().x()) - simMatrix.xLength() / 2.;
			currentProtonInstance.setXt((currentProtonInstance.xt().x() > 0 ? extra_length - simMatrix.xLength() / 2. : simMatrix.xLength() / 2. - extra_length), currentProtonInstance.xt().y(), currentProtonInstance.xt().z());
		}
		if (abs(currentProtonInstance.xt().y()) > simMatrix.yLength() / 2.) {
			extra_length = abs(currentProtonInstance.xt().y()) - simMatrix.yLength() / 2.;
			currentProtonInstance.setXt(currentProtonInstance.xt().x(), (currentProtonInstance.xt().y() > 0 ? extra_length - simMatrix.yLength() / 2. : simMatrix.yLength() / 2. - extra_length), currentProtonInstance.xt().z());
		}
		if (abs(currentProtonInstance.xt().z()) > simMatrix.zLength() / 2.) {
			extra_length = abs(currentProtonInstance.xt().z()) - simMatrix.zLength() / 2.;
			currentProtonInstance.setXt(currentProtonInstance.xt().x(), currentProtonInstance.xt().y(), (currentProtonInstance.xt().z() > 0 ? extra_length - simMatrix.zLength() / 2. : simMatrix.zLength() / 2. - extra_length));
		}

		for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
			if (distance(currentProtonInstance.xt(), ionpInstances[k].xt()) < (ionpInstances[k].radius() + ionpInstances[k].coating())) {
				repeat = true;
				break;
			}
		}


		if (!repeat) { //correct position found
			protonInstances[i].setDiffTime(j * dt); //j=0 so at t=0, no time
			break;
		}

	}
	protonInstances[i] = currentProtonInstance;
}

//Function to replace IONP to simulate their movement
void MRsequence::replaceIONP( std::vector<Ionp>& ionpInstances ){
	for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
		ionpInstances[k].move(directionIONP, vIONP, time_step);
	}
}


//Main function here - perform proton simulation and phase/signal calculation
void MRsequence::Sequence(std::vector<Proton>& protonInstances, std::vector<Ionp>& ionpInstances_orig, std::vector<std::complex<double>>& signalVector, std::vector<double>& timeV, std::vector<std::complex<double>>& signalVectorMSE, std::vector<double>& timeVMSE) {

	//Use maximum number available threads if no number is given - for parallelization
	if (numThreads == 0){
    	numThreads = omp_get_max_threads( );
		std::cout << "Number of max Threads: " << omp_get_max_threads( ) << std::endl;
  	}

	//Calculate Magnetic field on set num of nodes
	std::vector< std::vector< std::vector<double> > > BfieldGrid(SC_UI(matrix_nxnodes) , std::vector< std::vector<double> > (SC_UI(matrix_nxnodes), std::vector<double> (SC_UI(matrix_nxnodes)) ) );

	if (Grid == true){
		simMatrix.calculateBfield(BfieldGrid, ionpInstances_orig, numThreads, magnetic_moment);
	}

	std::cout << "Real Time Steps: " << time_step << std::endl;
	std::cout << "Real Radius with coating: " << ionpInstances_orig[0].radius() + ionpInstances_orig[0].coating() << std::endl;
	std::cout << "Magnetization: " << D_M_IONP << std::endl;

	//==================== Proton diffusion ==========================//

	omp_set_num_threads( numThreads );

	const std::complex<double> i_compl(0, 1); // i

	unsigned int all_proton_num = SC_UI(protons_n) + SC_UI(ceil(SC_D(bound_protons_n)/SC_D(ionpInstances_orig.size()))) * SC_UI(ionpInstances_orig.size());

	std::cout << "Total number of protons simulated: " << all_proton_num << std::endl;

	//Function to determine ideal std for specific T2 of Water/Agar
	double standarddeviation_back = (0.14218444427715876 * pow((time_step/0.00001), 0.5014107086422922)) * pow(background_T2*1e3, -0.5014107088655116);

	std::vector<std::vector<int>> protonPosition;

	for (unsigned int i = 0; i < all_proton_num; ++i) {
		protonInstances.push_back(Proton(0, 0));
		protonPosition.push_back(std::vector<int>());
	}

	//Vector where every phase of every proton at every time is saved
	std::vector< std::vector<double> > phaseVector;

	if ( cube_s.x() != 0. ){
		if ( abs(cube_o.x() - cube_s.x()/2.) > simMatrix.xLength() / 2. || abs(cube_o.x() + cube_s.x()/2.) > simMatrix.xLength() / 2. || abs(cube_o.y() - cube_s.y()/2.) > simMatrix.yLength() / 2. || abs(cube_o.y() + cube_s.y()/2.) > simMatrix.yLength() / 2. || abs(cube_o.z() - cube_s.z()/2.) > simMatrix.zLength() / 2. || abs(cube_o.z() + cube_s.z()/2.) > simMatrix.zLength() / 2.) {
			std::cout << "Cube is outside of Simulation Volume" << std::endl;
		}
	}


	for (unsigned int i = 0; i < all_proton_num; ++i) {
		phaseVector.push_back(std::vector<double>());
	}


	//Proton loop (every proton is simulated seperately)

	std::cout << "Now protons are propagating..." << std::endl;

	#pragma omp parallel
	{
		//Define random numbers depending on desired distribution 
		RandomList RandomNormalX( -simMatrix.xLength() / 2., simMatrix.xLength() / 2., false );
		RandomList RandomNormalY( -simMatrix.yLength() / 2., simMatrix.yLength() / 2., false );
		RandomList RandomNormalZ( -simMatrix.zLength() / 2., simMatrix.zLength() / 2., false );
		RandomList RandomDiff( -1., 1., false );

		#pragma omp for
		for (unsigned int i = 0; i < SC_UI(protons_n); ++i) {

			//================================= Start Initialization =================================//

			//Initiate Proton specific variables
			double dephase = 0., phase = 0., time = 0.;
			unsigned int j = 0, leapfrog = 0;
			int counter = 1;

			std::vector<double> phaseOneProt;
			std::vector<int> positionOneProt;

			//Needed if IONP move because every proton is viewed seperately
			std::vector<Ionp> ionpInstances;
			ionpInstances = ionpInstances_orig;

			//Initiate Gaussian distribution for background T2
			RandNumberGaussian rand_phase_back(0, standarddeviation_back);

			//Place proton random
			replaceProton(protonInstances, ionpInstances, time_step, i, j, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff);


			if ( cube_s.x() != 0. ){
				positionOneProt.push_back(0);
				if ( protonInstances[i].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[i].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[i].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[i].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[i].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[i].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
					positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
				}
			}


			//Start Signal at t=0

			phaseOneProt.push_back(phase);
			if (!i) {
				timeV.push_back(time);
			}

			//================================= End Initialization =================================//

			//================================= Start Proton (IONP) propagation =================================//
			while(time < echo_time){

				//Reset Time step size for every iteration if a variable Time step is used
				double dt = time_step;
				bool smallerTimeStep = false;
				double bfield_tot = 0.;
				double agar_relax = 0.;


				if (leapfrog % 2 == 0){

					if ( var_Steps == true ) {
						for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
							if (distance(protonInstances[i].xt(), ionpInstances[k].xt()) < 8. * (ionpInstances[k].radius() + ionpInstances[k].coating())) {
								dt /= 64.;
								smallerTimeStep = true;
								j++;
								time += time_step;
								for (unsigned int s = 1; s <= 64; ++s){

									replaceProton(protonInstances, ionpInstances, dt, i, j, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff);

									bfield_tot = 0.;
									agar_relax = 0.;

									if (magnetic_moment != 0.){
										for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
											bfield_tot += simMatrix.calculateBfield(protonInstances[i].xt(), ionpInstances[l].xt(), ionpInstances[l].radius(), magnetic_moment);
										}
									}
									else{
										for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
											bfield_tot += simMatrix.calculateBfield(protonInstances[i].xt(), ionpInstances[l].xt(), ionpInstances[l].radius());
										}
									}


									if (background_T2 != 0.){
										agar_relax = rand_phase_back();
									}


									dephase = D_GI* bfield_tot * dt + agar_relax; //SI

									phase += dephase;
								}
								if ( sequence_num == 1 ){
									if (j == SC_UI(round(echo_spacing / (2*time_step)))){
        								phase *= -1.;
		 							}
								}
								else if ( sequence_num == 2){
									for (unsigned l = 0; l < SC_UI(round(echo_time/echo_spacing)); l++){
										if (j == SC_UI(round((2*l+1)*echo_spacing / (2*time_step)))){
        									phase *= -1.;
		 								}
									}
								}

								if ( time_step > 0.00001 && echo_time <= 0.1){
									phaseOneProt.push_back(phase);
									if (!i) {
										timeV.push_back(time);
									}
									if ( cube_s.x() != 0. ){
										positionOneProt.push_back(0);
										if ( protonInstances[i].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[i].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[i].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[i].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[i].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[i].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
											positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
										}
									}
								}
								else if (j == SC_UI(round(counter*0.00001 / time_step * echo_time * 10. * all_proton_num / 10000. ))){
									phaseOneProt.push_back(phase);
									if (!i) {
										timeV.push_back(time);
									}
									if ( cube_s.x() != 0. ){
										positionOneProt.push_back(0);
										if ( protonInstances[i].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[i].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[i].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[i].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[i].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[i].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
											positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
										}
									}
									counter ++;
								}

								break;
							}
						}
					}


					if (smallerTimeStep == false){
						time += time_step;
						j++;

						replaceProton(protonInstances, ionpInstances, time_step, i, j, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff);

						//============= Calculation of the Bfield at proton position =================================//
						bfield_tot = 0;

						if (Grid == true){
							bfield_tot = simMatrix.interpolate(protonInstances[i].xt(), BfieldGrid);
						}
						else{
							if (magnetic_moment != 0.){
								for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
									bfield_tot += simMatrix.calculateBfield(protonInstances[i].xt(), ionpInstances[l].xt(), ionpInstances[l].radius(), magnetic_moment);
								}
							}
							else{
								for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
									bfield_tot += simMatrix.calculateBfield(protonInstances[i].xt(), ionpInstances[l].xt(), ionpInstances[l].radius());
								}
							}
						}

						agar_relax = 0.;

						if (background_T2 != 0.){
							agar_relax = rand_phase_back();
						}

						dephase = D_GI* bfield_tot * time_step + agar_relax;

						phase += dephase; //phase at t

						if ( sequence_num == 1 ){
							if (j == SC_UI(round(echo_spacing / (2*time_step)))){
        						phase *= -1.;
		 					}
						}
						else if ( sequence_num == 2){
							for (unsigned k = 0; k < SC_UI(round(echo_time/echo_spacing)); k++){
								if (j == SC_UI(round((2*k+1)*echo_spacing / (2*time_step)))){
        							phase *= -1.;
		 						}
							}
						}


						if ( time_step > 0.00001 && echo_time <= 0.1){
							phaseOneProt.push_back(phase);
							if (!i) {
								timeV.push_back(time);
							}
							if ( cube_s.x() != 0. ){
								positionOneProt.push_back(0);
								if ( protonInstances[i].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[i].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[i].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[i].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[i].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[i].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
									positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
								}
							}
						}
						else if (j == SC_UI(round(counter*0.00001 / time_step * echo_time * 10. * all_proton_num / 10000. ))){
							phaseOneProt.push_back(phase);
							if (!i) {
								timeV.push_back(time);
							}
							if ( cube_s.x() != 0. ){
								positionOneProt.push_back(0);
								if ( protonInstances[i].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[i].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[i].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[i].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[i].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[i].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
									positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
								}
							}
							counter ++;
						}

					}

				}
				else{
					//IONP movement can only be simulated if the B-Field is calculated directly
					if (vIONP != 0. && Grid == false){
						replaceIONP(ionpInstances);
					}
				}

				if (vIONP != 0. && Grid == false){
					leapfrog ++;
				}

			}

			phaseVector[i] = phaseOneProt;
			protonPosition[i] = positionOneProt;

		}

	}

	//================================= Now dephasing of bound protons at the shell of the IONPs =================================//

	unsigned int bound_protons_per_ionp = SC_UI(ceil(SC_D(bound_protons_n)/SC_D(ionpInstances_orig.size())));

	#pragma omp parallel for
	for (unsigned int i = 0; i < bound_protons_per_ionp; ++i) { //protons loop
		//Needed if IONP move because every proton is viewed seperately
		std::vector<Ionp> ionpInstances;
		ionpInstances = ionpInstances_orig;

		if ( i == 0){
			std::cout << "There are bound protons" << std::endl;
		}

		for (unsigned int v = 0; v < SC_UI(ionpInstances.size()); ++v) {

			//================================= Start Initialization =================================//

			//Initiate Proton specific variables
			double dephase = 0., phase = 0., time = 0.;
			unsigned int j = 0, leapfrog = 0;//, sigTime = 0;
			int counter = 1;

			std::vector<double> phaseOneProt;
			std::vector<int> positionOneProt;

			//Initiate Gaussian distribution for background T2
			RandNumberGaussian rand_phase_back(0, standarddeviation_back);

			//Place proton
			protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].placeRandomSphere(ionpInstances[v].x0().x(), ionpInstances[v].x0().y(), ionpInstances[v].x0().z(), ionpInstances[v].radius(), (ionpInstances[v].radius() + ionpInstances[v].coating()), true);

			if ( cube_s.x() != 0. ){
				positionOneProt.push_back(0);
				if ( protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
					positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
				}
			}

			//Start Signal at t=0

			phaseOneProt.push_back(phase);
			if (protons_n == 0){
				if (!i) {
					timeV.push_back(time);
				}
			}

			//================================= End Initialization =================================//


			//================================= Start Proton (IONP) propagation =================================//
			while(time < echo_time){

				//Reset Time step size for every iteration if a variable Time step is used
				double dt = time_step;
				bool smallerTimeStep = false;
				double bfield_tot = 0.;
				double agar_relax = 0.;

				if (leapfrog % 2 == 0){

					if ( var_Steps == true ) {
						for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
							if (distance(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), ionpInstances[k].xt()) < 8. * (ionpInstances[k].radius() + ionpInstances[k].coating())) {
								dt /= 64.;
								smallerTimeStep = true;
								j++;
								time += time_step;
								for (unsigned int s = 0; s < 64; ++s){

									bfield_tot = 0;
									agar_relax = 0.;

									if (magnetic_moment != 0.){
										for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
											bfield_tot += simMatrix.calculateBfield(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), ionpInstances[l].xt(), ionpInstances[l].radius(), magnetic_moment);
										}
									}
									else{
										for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
											bfield_tot += simMatrix.calculateBfield(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), ionpInstances[l].xt(), ionpInstances[l].radius());
										}
									}

									if (background_T2 != 0.){
										agar_relax = rand_phase_back();
									}

									dephase = D_GI* bfield_tot * dt + agar_relax; //SI

									phase += dephase;
								}

								if ( time_step > 0.00001 && echo_time <= 0.1){
									phaseOneProt.push_back(phase);
									if (protons_n == 0){
										if (!i) {
											timeV.push_back(time);
										}
									}
									if ( cube_s.x() != 0. ){
										positionOneProt.push_back(0);
										if ( protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
											positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
										}
									}
								}
								else if (j == SC_UI(round(counter*0.00001 / time_step * echo_time * 10. * all_proton_num / 10000. ))){
									phaseOneProt.push_back(phase);
									if (protons_n == 0){
										if (!i) {
											timeV.push_back(time);
										}
									}
									if ( cube_s.x() != 0. ){
										positionOneProt.push_back(0);
										if ( protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
											positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
										}
									}
									counter ++;
								}

								break;
							}
						}
					}

					(void) dt;
					(void) smallerTimeStep;

					time += time_step;
					j++;

					//============= Calculation of the Bfield at proton position =================================//
					bfield_tot = 0;

					if (Grid == true){
						bfield_tot = simMatrix.interpolate(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), BfieldGrid);
					}
					else{
						if (magnetic_moment != 0.){
							for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
								bfield_tot += simMatrix.calculateBfield(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), ionpInstances[l].xt(), ionpInstances[l].radius(), magnetic_moment);
							}
						}
						else{
							for (unsigned l = 0; l < SC_UI(ionpInstances.size()); ++l) {
								bfield_tot += simMatrix.calculateBfield(protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt(), ionpInstances[l].xt(), ionpInstances[l].radius());
							}
						}
					}

					agar_relax = 0.;

					if (background_T2 != 0.){
						agar_relax = rand_phase_back();
					}

					dephase = D_GI* bfield_tot * time_step + agar_relax;

					phase += dephase; //phase at t

					if ( sequence_num == 1 ){
						if (j == SC_UI(round(echo_spacing / (2*time_step)))){
        					phase *= -1.;
		 				}
					}
					else if ( sequence_num == 2){
						for (unsigned k = 0; k < SC_UI(round(echo_time/echo_spacing)); k++){
							if (j == SC_UI(round((2*k+1)*echo_spacing / (2*time_step)))){
        						phase *= -1.;
		 					}
						}
					}

					if ( time_step > 0.00001 && echo_time <= 0.1){
						phaseOneProt.push_back(phase);
						if (protons_n == 0){
							if (!i) {
								timeV.push_back(time);
							}
						}
						if ( cube_s.x() != 0. ){
							positionOneProt.push_back(0);
							if ( protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
								positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
							}
						}
					}
					else if (j == SC_UI(round(counter*0.00001 / time_step * echo_time * 10. * all_proton_num / 10000. ))){
						phaseOneProt.push_back(phase);
						if (protons_n == 0){
							if (!i) {
								timeV.push_back(time);
							}
						}
						if ( cube_s.x() != 0. ){
							positionOneProt.push_back(0);
							if ( protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() > (cube_o.x() - cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().x() < (cube_o.x() + cube_s.x()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() >(cube_o.y() - cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().y() < (cube_o.y() + cube_s.y()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() > (cube_o.z() - cube_s.z()/2.) && protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].xt().z() < (cube_o.z() + cube_s.z()/2.)) {
								positionOneProt[SC_UI(positionOneProt.size()-1)] = 1;
							}
						}
						counter ++;
					}
				}
				else{
					//IONP movement can only be simulated if the B-Field is calculated directly
					if (vIONP != 0. && Grid == false){
						replaceIONP(ionpInstances);
						protonInstances[SC_UI(protons_n) + i + bound_protons_per_ionp*v].move(directionIONP, vIONP, time_step);
					}
				}

				leapfrog ++;

			}
			phaseVector[SC_UI(protons_n) + i + bound_protons_per_ionp*v] = phaseOneProt;
			protonPosition[SC_UI(protons_n) + i + bound_protons_per_ionp*v] = positionOneProt;
		}
	}

	std::cout << "Magnetization calculation start!" << std::endl;

	for (unsigned int i = 0; i < all_proton_num; ++i) {
		for (unsigned int j = 0; j < phaseVector[0].size(); ++j){
			if (!i) {
				signalVector.push_back(std::exp(i_compl * phaseVector[i][j])/SC_D(all_proton_num));
			}
			else {
				signalVector[j] += std::exp(i_compl * phaseVector[i][j])/SC_D(all_proton_num);
			}
		}
	}

	if ( sequence_num == 2){
		for (unsigned int j = 0; j < signalVector.size(); ++j){
			for (unsigned k = 0; k < SC_UI(round(echo_time/echo_spacing)); k++){
				if ( time_step <= 0.00001 ){
					if (j == SC_UI(round(k*echo_spacing / 0.00001))){
						signalVectorMSE.push_back(signalVector[j]);
						timeVMSE.push_back(timeV[j]);
					}
				}
				else{
					if (j == SC_UI(round(k*echo_spacing / time_step))){
						signalVectorMSE.push_back(signalVector[j]);
						timeVMSE.push_back(timeV[j]);
					}
				}
			}
		}
	}

	//get the Magnitude, so calculate abs of complex signal
	std::for_each(signalVector.begin(), signalVector.end(), [](std::complex<double>& c) { c = std::abs(c); });
	std::for_each(signalVectorMSE.begin(), signalVectorMSE.end(), [](std::complex<double>& c) { c = std::abs(c); });
	std::cout << "Magnetization calculated." << std::endl;

	//option if a smaller cube inside the simulation volume is defined
	if ( cube_s.x() != 0. ){
		std::vector<int> protPos1D;

		if (protonPosition.size() != 0){
			for (unsigned int j = 0; j < protonPosition.size(); ++j){
				for (unsigned int i = 0; i < protonPosition[0].size(); ++i){
					protPos1D.push_back(protonPosition[j][i]);
				}
			}
		}

		std::vector<std::complex<double>> signalShellProtons(phaseVector[0].size(), 0.);
		long int numberShellProtons = std::count(protPos1D.begin(), protPos1D.end(), 1);

		std::cout << "Number of Protons*TimeStep in defined Voxel: " << numberShellProtons << std::endl;

		for (unsigned int i = 0; i < all_proton_num; ++i) {
			for (unsigned int j = 0; j < phaseVector[0].size(); ++j){
				if ( protonPosition[i][j] == 1 ){
					signalShellProtons[j] += std::exp(i_compl * phaseVector[i][j])/SC_D(all_proton_num);
				}
			}
		}

		std::for_each(signalShellProtons.begin(), signalShellProtons.end(), [](std::complex<double>& c) { c = std::abs(c); });


		std::string path_to_data = "../../../data/tmp/";
		int id = 0;

		for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
		    (void) entry;
		    id++;
		 }

		 std::string protonData_out_name;

		 if ( sequence_num == 0 ){
		   protonData_out_name = path_to_data +  std::to_string( id ) + "_SignalCube_FID_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
		 }
		 else if ( sequence_num == 1 ){
		   protonData_out_name = path_to_data +  std::to_string( id ) + "_SignalCube_SE_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
		 }
		 else if ( sequence_num == 2 ){
		   protonData_out_name = path_to_data +  std::to_string( id ) + "_SignalCube_MSE_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
		 }

		 std::cout << protonData_out_name << std::endl;

		 std::ofstream outfile(protonData_out_name, std::ios::out | std::ios::binary);
		 outfile.precision(15);
		 if (!outfile.is_open()) {
			std::cerr << "Unable to open output file: SignalShell" << std::endl;
			exit(1);
		}
		for (unsigned int i = 0; i < signalShellProtons.size(); ++i) {
			outfile << real(signalShellProtons[i]) << "\t" << timeV[i] << std::endl;
		}
		outfile.close();
	}

}
