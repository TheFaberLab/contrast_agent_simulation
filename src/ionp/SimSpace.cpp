/*****************************************************************************
*Class to place IONP in SimVolume in different spacial distributions

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

#include "Point3D.hpp"
#include "Vector3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"

#include "SimSpace.hpp"


void SimSpace::simulationSpace(std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int ionp_n, double conc, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double thickness_coating, double thickness_coating_std, double percentage_ionp, bool DLS_file ) {


	std::cout << "Creating simulation space...." << std::endl;

	std::vector<double> rand_DLSrad;

	if ( DLS_file == true ){
		//Read the DLS measurements File of IONPs
		std::string path_to_DLS = "../../../ini/DLS_File";

		Data dat(path_to_DLS);

		for (const auto & entryDLS : std::filesystem::directory_iterator(path_to_DLS)){
			std::string pathDLS = entryDLS.path();
			if (pathDLS.find(".txt") != std::string::npos || pathDLS.find(".dat") != std::string::npos ){
				dat.readHist( pathDLS, rand_DLSrad );
			}
			else{
				std::cerr << "No DLS file found!" << std::endl;
				exit(1);
			}
		}

		for (unsigned int i = 0; i < SC_UI(rand_DLSrad.size()) - 1; ++i){
			rand_DLSrad[i] = rand_DLSrad[i] / 2. * 1.0e-9;
		}
	}

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberLognorm rand_coat(log(thickness_coating*2.), thickness_coating_std);
	RandNumberBetween rand_norm(0., 1.);

	if ( conc != 0. ){
		double concentration;

		concentration = conc * 5.5845e-5 * (simMatrix.xLength()*simMatrix.yLength()*simMatrix.zLength()*1000.);

		double actual_c = 0.;
		unsigned int i = 0, small = 0, big = 0;

		while ( actual_c <= concentration ){
			if ( DLS_file == true ){
				double coat;
				r_ionp = rand_DLSrad[i];
				if ( thickness_coating == 0. ) {
					coat = 0.;
				}
				else{
					if (thickness_coating_std != 0.){
						double coat_rand = rand_coat()/2.;
						if ( r_ionp >= coat_rand )
						{
							coat = 0.;
						}
						else
						{
							coat = coat_rand - r_ionp;
						}
					}
					else{
						if ( r_ionp >= thickness_coating )
						{
							coat = 0.;
						}
						else
						{
							coat = thickness_coating - r_ionp;
						}
					}
				}
				ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
			}
			else{
				if ( r_ionp_std == 0. ){
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						if (thickness_coating_std != 0.){
							double coat_rand = rand_coat()/2.;
							if ( r_ionp >= coat_rand )
							{
								coat = 0.;
							}
							else
							{
								coat = coat_rand - r_ionp;
							}
						}
						else{
							if ( r_ionp >= thickness_coating )
							{
								coat = 0.;
							}
							else
							{
								coat = thickness_coating - r_ionp;
							}
						}
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
				}
				else{
					if (r_ionp_1 != 0.){
						double prop = rand_norm();
						if (prop <= percentage_ionp){
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
							small += 1;
						}
						else{
							double r = rand1()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
							big += 1;
						}
					}
					else{
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							if (thickness_coating_std != 0.){
								double coat_rand = rand_coat()/2.;
								if ( r >= coat_rand )
								{
									coat = 0.;
								}
								else
								{
									coat = coat_rand - r;
								}
							}
							else{
								if ( r >= thickness_coating )
								{
									coat = 0.;
								}
								else
								{
									coat = thickness_coating - r;
								}
							}
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
				}
			}




				while (true) { // random placement and validation of position
					bool hit_ionp = false; //flag if proton enters IONP

					ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.xLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), -(simMatrix.yLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.yLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), -(simMatrix.zLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.zLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), true);

					for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
						if (i == 0) {
							break;//skip this validation for the first ionp
						}
						if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (ionpInstances[i].radius() + ionpInstances[i].coating() + ionpInstances[k].radius() + ionpInstances[k].coating())) {
							hit_ionp = true;
							break;
						}
					}
					if (!hit_ionp) {
						break;
					}
				}

			actual_c += 4./3. * D_PI * 5.24 * (167.4/231.4) * cube(ionpInstances[i].radius()*100.);
			i++;
		}
		std::cout << "Small: " << small << std::endl;
		std::cout << "Big: " << big << std::endl;
	}
	else{
		for (unsigned int i = 0; i < SC_UI(ionp_n); ++i) {

			if ( DLS_file == true ){
				double coat;
				r_ionp = rand_DLSrad[i];
				if ( thickness_coating == 0. ) {
					coat = 0.;
				}
				else{
					if (thickness_coating_std != 0.){
						double coat_rand = rand_coat()/2.;
						if ( r_ionp >= coat_rand )
						{
							coat = 0.;
						}
						else
						{
							coat = coat_rand - r_ionp;
						}
					}
					else{
						if ( r_ionp >= thickness_coating )
						{
							coat = 0.;
						}
						else
						{
							coat = thickness_coating - r_ionp;
						}
					}
				}
				ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
			}
			else{
				if ( r_ionp_std == 0. ){
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						if (thickness_coating_std != 0.){
							double coat_rand = rand_coat()/2.;
							if ( r_ionp >= coat_rand )
							{
								coat = 0.;
							}
							else
							{
								coat = coat_rand - r_ionp;
							}
						}
						else{
							if ( r_ionp >= thickness_coating )
							{
								coat = 0.;
							}
							else
							{
								coat = thickness_coating - r_ionp;
							}
						}
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
				}
				else{
					if (r_ionp_1 != 0.){
						double prop = rand_norm();
						if (prop <= percentage_ionp){
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						}
						else{
							double r = rand1()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						}
					}
					else{
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							if (thickness_coating_std != 0.){
								double coat_rand = rand_coat()/2.;
								if ( r >= coat_rand )
								{
									coat = 0.;
								}
								else
								{
									coat = coat_rand - r;
								}
							}
							else{
								if ( r >= thickness_coating )
								{
									coat = 0.;
								}
								else
								{
									coat = thickness_coating - r;
								}
							}
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
				}
			}


				while (true) { // random placement and validation of position
					bool hit_ionp = false; //flag if proton enters IONP

					ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.xLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), -(simMatrix.yLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.yLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), -(simMatrix.zLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), (simMatrix.zLength() / 2. - (ionpInstances[i].radius() + ionpInstances[i].coating())), true);

					for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
						if (i == 0) {
							break;//skip this validation for the first ionp
						}
						if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (ionpInstances[i].radius() + ionpInstances[i].coating() + ionpInstances[k].radius() + ionpInstances[k].coating())) {
							hit_ionp = true;
							break;
						}
					}
					if (!hit_ionp) {
						break;
					}
				}
		}
	}



	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}

void SimSpace::surroundingSpace(std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int ionp_n, double r_ionp ) {
	std::cout << "Creating additional 26 simulation volumes for periodicity" << std::endl;

	for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));

		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));

		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() + simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() + simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() + simMatrix.zLength(), r_ionp));
		ionpInstances.push_back(Ionp(ionpInstances[j].xt().x() - simMatrix.xLength(),ionpInstances[j].xt().y() - simMatrix.yLength(),ionpInstances[j].xt().z() - simMatrix.zLength(), r_ionp));
	}
}


void SimSpace::simulationSpace(std::vector<Aggregate>& aggregateInstances, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int agg_n, int ionp_n, double conc, double r_agg, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double r_agg_std, double thickness_coating, double thickness_coating_std, double percentage_ionp) {
	std::cout << "Creating simulation space containing spherical aggregates...." << std::endl;

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberLognorm rand_coat(log(thickness_coating*2.), thickness_coating_std);
	RandNumberBetween rand_norm(0., 1.);
	RandNumberLognorm randAgg(log(r_agg*2.), r_agg_std);

	double r_agg_rand;

	if (aggregateInstances.empty()){
		unsigned int ionp_agg_num = 0;
		for (unsigned int i = 0; i < SC_UI(agg_n); ++i) {
			unsigned int ionp_curr_agg = 0;
			if ( r_agg_std == 0. ){
				aggregateInstances.push_back(Aggregate(r_agg, ionp_n));
				aggregateInstances[i].placeRandom(-(simMatrix.xLength() / 2. - r_agg), (simMatrix.xLength() / 2. - r_agg), -(simMatrix.yLength() / 2. - r_agg), (simMatrix.yLength() / 2. - r_agg), -(simMatrix.zLength() / 2. - r_agg), (simMatrix.zLength() / 2. - r_agg), true);
			}
			else{
				r_agg_rand = randAgg()/2.;
				aggregateInstances.push_back(Aggregate(r_agg_rand, ionp_n));
				aggregateInstances[i].placeRandom(-(simMatrix.xLength() / 2. - r_agg_rand), (simMatrix.xLength() / 2. - r_agg_rand), -(simMatrix.yLength() / 2. - r_agg_rand), (simMatrix.yLength() / 2. - r_agg_rand), -(simMatrix.zLength() / 2. - r_agg_rand), (simMatrix.zLength() / 2. - r_agg_rand), true);
			}

			if ( conc != 0. ){
				double concentration;
				concentration = conc * 5.5845e-5 * (simMatrix.xLength()*simMatrix.yLength()*simMatrix.zLength()*1000.); 
				double actual_c = 0.;
				unsigned int j = 0, small = 0, big = 0;

				std::cout << "Aggregate " << i << ":" << std::endl;

				while ( actual_c <= concentration/agg_n ){
					if ( r_ionp_std == 0. ){
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							if (thickness_coating_std != 0.){
								double coat_rand = rand_coat()/2.;
								if ( r_ionp >= coat_rand )
								{
									coat = 0.;
								}
								else
								{
									coat = coat_rand - r_ionp;
								}
							}
							else{
								if ( r_ionp >= thickness_coating )
								{
									coat = 0.;
								}
								else
								{
									coat = thickness_coating - r_ionp;
								}
							}
						}
						ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp, coat));
					}
					else{
						if (r_ionp_1 != 0.){
							double prop = rand_norm();
							if (prop <= percentage_ionp){
								double r = rand()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									if (thickness_coating_std != 0.){
										double coat_rand = rand_coat()/2.;
										if ( r >= coat_rand )
										{
											coat = 0.;
										}
										else
										{
											coat = coat_rand - r;
										}
									}
									else{
										if ( r >= thickness_coating )
										{
											coat = 0.;
										}
										else
										{
											coat = thickness_coating - r;
										}
									}
								}
								ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
								small += 1;
							}
							else{
								double r = rand1()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									if (thickness_coating_std != 0.){
										double coat_rand = rand_coat()/2.;
										if ( r >= coat_rand )
										{
											coat = 0.;
										}
										else
										{
											coat = coat_rand - r;
										}
									}
									else{
										if ( r >= thickness_coating )
										{
											coat = 0.;
										}
										else
										{
											coat = thickness_coating - r;
										}
									}
								}
								ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
								big += 1;
							}
						}
						else{
							double r;
							r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
							small += 1;
						}
					}

					while (true) { // random placement and validation of position
						bool hit_ionp = false; //flag if proton enters IONP ->update each time enters the loop for repositioning
						ionpInstances[ionp_agg_num].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., aggregateInstances[i].radius(), true);

						for (unsigned int k = 0; k < ionp_curr_agg; ++k) { //not inside IONP // don't check final entry(is this ionp)
							if (ionp_curr_agg == 0) {
								break;//skip this validation for the first ionp
							}
							if (distance(ionpInstances[ionp_agg_num].xt(), ionpInstances[ionp_agg_num - (k + 1)].xt()) < (ionpInstances[ionp_agg_num].radius() + ionpInstances[ionp_agg_num].coating() + ionpInstances[ionp_agg_num - (k + 1)].radius() + ionpInstances[ionp_agg_num - (k + 1)].coating())) {
								hit_ionp = true;
								break;
							}
						}
						if (!hit_ionp) {
							break;
						}
					}
					actual_c += 4./3. * D_PI * 5.24 * (167.4/231.4) * cube(ionpInstances[ionp_agg_num].radius()*100.);
					j++;
					ionp_agg_num++;
					ionp_curr_agg++;
				}
				std::cout << "Small: " << small << std::endl;
				std::cout << "Big: " << big << std::endl;
			}
			else{
				for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
					if ( r_ionp_std == 0. ){
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							if (thickness_coating_std != 0.){
								double coat_rand = rand_coat()/2.;
								if ( r_ionp >= coat_rand )
								{
									coat = 0.;
								}
								else
								{
									coat = coat_rand - r_ionp;
								}
							}
							else{
								if ( r_ionp >= thickness_coating )
								{
									coat = 0.;
								}
								else
								{
									coat = thickness_coating - r_ionp;
								}
							}
						}
						ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp, coat));
					}
					else{
						if (r_ionp_1 != 0.){
							double prop = rand_norm();
							if (prop <= percentage_ionp){
								double r = rand()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									if (thickness_coating_std != 0.){
										double coat_rand = rand_coat()/2.;
										if ( r >= coat_rand )
										{
											coat = 0.;
										}
										else
										{
											coat = coat_rand - r;
										}
									}
									else{
										if ( r >= thickness_coating )
										{
											coat = 0.;
										}
										else
										{
											coat = thickness_coating - r;
										}
									}
								}
								ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
							}
							else{
								double r = rand1()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									if (thickness_coating_std != 0.){
										double coat_rand = rand_coat()/2.;
										if ( r >= coat_rand )
										{
											coat = 0.;
										}
										else
										{
											coat = coat_rand - r;
										}
									}
									else{
										if ( r >= thickness_coating )
										{
											coat = 0.;
										}
										else
										{
											coat = thickness_coating - r;
										}
									}
								}
								ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
							}
						}
						else{
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								if (thickness_coating_std != 0.){
									double coat_rand = rand_coat()/2.;
									if ( r >= coat_rand )
									{
										coat = 0.;
									}
									else
									{
										coat = coat_rand - r;
									}
								}
								else{
									if ( r >= thickness_coating )
									{
										coat = 0.;
									}
									else
									{
										coat = thickness_coating - r;
									}
								}
							}
							ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
						}
					}

					while (true) { // random placement and validation of position
						bool hit_ionp = false; //flag if proton enters IONP ->update each time enters the loop for repositioning
						ionpInstances[i * SC_UI(ionp_n) + j].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., aggregateInstances[i].radius(), true);

						for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
							if (j == 0 && i == 0) {
								break;//skip this validation for the first ionp
							}
							if (distance(ionpInstances[i * SC_UI(ionp_n) + j].xt(), ionpInstances[k].xt()) < (ionpInstances[i * SC_UI(ionp_n) + j].radius() + ionpInstances[i * SC_UI(ionp_n) + j].coating() + ionpInstances[k].radius() + ionpInstances[k].coating())) {
								hit_ionp = true;
								break;
							}
						}
						if (!hit_ionp) {
							break;
						}
					}
				}
			}
		}
	}
	else{
		for (unsigned int i = 0; i <  SC_UI(aggregateInstances.size()); ++i) {

		if (r_agg_std != 0.){
			r_agg_rand = randAgg();
		}

		for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
			if ( r_ionp_std == 0. ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else{
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., rand()));
			}

			while (true) { // random placement and validation of position
				bool hit_ionp = false; //flag if proton enters IONP ->update each time enters the loop for repositioning

				if ( r_agg_std == 0.){
					ionpInstances[i * SC_UI(ionp_n) + j].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_agg, true);
				}
				else{
					ionpInstances[i * SC_UI(ionp_n) + j].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_agg_rand, true);
				}


				for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
					if (j == 0 && i == 0) {
						break;//skip this validation for the first ionp
					}
					if (distance(ionpInstances[i * SC_UI(ionp_n) + j].xt(), ionpInstances[k].xt()) < (ionpInstances[i * SC_UI(ionp_n) + j].radius() + ionpInstances[k].radius())) {
						hit_ionp = true;
						break;
					}
				}
				if (!hit_ionp) {
					break;
				}
			}

		}

		}
	}

	std::cout << "Aggregates created: " << aggregateInstances.size() << std::endl;
	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}


//At the moment no manual placement of dense aggregates possible
void SimSpace::simulationSpace(std::vector<Aggregate>& aggregateInstances, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int agg_n, int ionp_n, double conc, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double thickness_coating, double percentage_ionp){
	std::cout << "Creating simulation space spherical dense aggregates...." << std::endl;

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberBetween rand_norm(0., 1.);

	double concentration;
	double number_parts;
	
	concentration = conc * 5.5845e-5 * (simMatrix.xLength()*simMatrix.yLength()*simMatrix.zLength()*1000.);
	if ( conc != 0.){
		if (r_ionp_1 != 0.){
			number_parts = concentration / (4./3. * D_PI * 5.24 * (167.4/231.4) * cube((r_ionp+r_ionp_1)/2. *100.) * agg_n);
		}
		else{
			number_parts = concentration / (4./3. * D_PI * 5.24 * (167.4/231.4) * cube(r_ionp*100.) * agg_n);
		}
	}
	else{
		number_parts = ionp_n;
	}


	int num_layers = SC_I(ceil(std::pow( number_parts, 1./3.)));

	double r_agg;

	if (r_ionp_1 != 0.){
			r_agg = (num_layers/2. + 2.) * r_ionp_1;
		}
		else{
			r_agg = (num_layers/2. + 2.) * r_ionp;
		}


	std::vector< std::vector<Ionp> > ionpVector;

	if (aggregateInstances.empty()){
		for (unsigned int i = 0; i < SC_UI(agg_n); ++i) {

			ionpVector.push_back(std::vector<Ionp>());

			aggregateInstances.push_back(Aggregate(r_agg, ionp_n));


			while (true) {
				bool hit_agg = false; //flag if agg overlap

				aggregateInstances[i].placeRandom(-(simMatrix.xLength() / 2. - aggregateInstances[i].radius()), (simMatrix.xLength() / 2. - aggregateInstances[i].radius()), -(simMatrix.yLength() / 2. - aggregateInstances[i].radius()), (simMatrix.yLength() / 2. - aggregateInstances[i].radius()), -(simMatrix.zLength() / 2. - aggregateInstances[i].radius()), (simMatrix.zLength() / 2. - aggregateInstances[i].radius()), true);

				for (unsigned int k = 0; k < SC_UI(aggregateInstances.size()) - 1; ++k) {
					if (i == 0) {
						break;//skip this validation for the first ionp
					}
					if (distance(aggregateInstances[i].xt(), aggregateInstances[k].xt()) < (aggregateInstances[i].radius() + aggregateInstances[k].radius())) {
						hit_agg = true;
						break;
					}
				}
				if (!hit_agg) {
					break;
				}
			}

		}

		if ( thickness_coating == 0. ) {
			thickness_coating = r_ionp;
		}

		for (unsigned int i = 0; i < SC_UI(agg_n); ++i) {
			std::vector<Ionp> ionpVector_oneAgg;
			if ( conc != 0. ){
				double actual_c = 0.;
				unsigned int j = 0, small = 0, big = 0;

				int layer = -num_layers/2;
				int row = -num_layers/2;
				int col = -num_layers/2;
				int count_layer = 1;
				int count_row = 1;
				int count_col = 1;

				while ( actual_c <= concentration/agg_n ){
					if ( r_ionp_std == 0. ){
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = thickness_coating - r_ionp;
						}
						ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp, coat));
					}
					else{
						if (r_ionp_1 != 0.){
							double prop = rand_norm();
							if (prop <= percentage_ionp){
								double r = rand()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									coat = thickness_coating - r;
								}
								ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
								small += 1;
							}
							else{
								double r = rand1()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									coat = thickness_coating - r;
								}
								ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
								big += 1;
							}
						}
						else{
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = thickness_coating - r;
							}
							ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
						}
					}

					if(j == 0){
						ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].setXt(aggregateInstances[i].x0().x() + col * 2. * thickness_coating, aggregateInstances[i].x0().y() + row * 2. * thickness_coating, aggregateInstances[i].x0().z() + layer * 2. * thickness_coating);
					}
					else if ( count_row >= num_layers ){
						layer++;
						count_layer++;
						row = -num_layers/2;
						count_row = 1;
						count_col = 1;
						ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].setXt(aggregateInstances[i].x0().x() + col * 2. * thickness_coating, aggregateInstances[i].x0().y() + row * 2. * thickness_coating, aggregateInstances[i].x0().z() + layer * 2. * thickness_coating);
					}
					else if ( count_col >= num_layers ){
						row++;
						count_row++;
						count_col = 1;
						ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].setXt(aggregateInstances[i].x0().x() + col * 2. * thickness_coating, aggregateInstances[i].x0().y() + row * 2. * thickness_coating, aggregateInstances[i].x0().z() + layer * 2. * thickness_coating);
					}
					else{
						ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].setXt(ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 2)].xt().x() + 2. * thickness_coating, aggregateInstances[i].x0().y() + row * 2. * thickness_coating, aggregateInstances[i].x0().z() + layer * 2. * thickness_coating);
						count_col++;
					}


					for (unsigned int k = 0; k < SC_UI(ionpVector_oneAgg.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
						if (j == 0 && i == 0) {
							break;//skip this validation for the first ionp
						}
						if (distance(ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].xt(), ionpVector_oneAgg[k].xt()) + 1e-9 < (ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].radius() + ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].coating() + ionpVector_oneAgg[k].radius() + ionpVector_oneAgg[k].coating())) {
							std::cout << "Something overlaps..." << std::endl;
							std::cout << SC_UI(ionpVector_oneAgg.size() - 1) << " " << k << std::endl;
							std::cout << ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].xt().x() << " " << ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].xt().y() << " " << ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].xt().z() << std::endl;
							std::cout << ionpVector_oneAgg[k].xt().x() << " " << ionpVector_oneAgg[k].xt().y() << " " << ionpVector_oneAgg[k].xt().z() << std::endl;
						}
					}


					actual_c += 4./3. * D_PI * 5.24 * (167.4/231.4) * cube(ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].radius()*100.);
					j++;
				}
				std::cout << "Aggregate " << i << ":" << std::endl;
				std::cout << "Small: " << small << std::endl;
				std::cout << "Big: " << big << std::endl;
				std::cout << "Total: " << ionpVector_oneAgg.size() << std::endl;
			}
			else{
				for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
					if ( r_ionp_std == 0. ){
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = thickness_coating - r_ionp;
						}
						ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp, coat));
					}
					else{
						if (r_ionp_1 != 0.){
							double prop = rand_norm();
							if (prop <= percentage_ionp){
								double r = rand()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									coat = thickness_coating - r;
								}
								ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
							}
							else{
								double r = rand1()/2.;
								double coat;
								if ( thickness_coating == 0. ) {
									coat = 0.;
								}
								else{
									coat = thickness_coating - r;
								}
								ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
							}
						}
						else{
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = thickness_coating - r;
							}
							ionpVector_oneAgg.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r, coat));
						}
					}

					while (true) { // random placement and validation of position
						bool hit_ionp = false; //flag if proton enters IONP ->update each time enters the loop for repositioning
						ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., aggregateInstances[i].radius(), true);

						//At the moment: checks if IONP overlaps with any other IONP independent of aggregate (faster is only inside one aggregate is checkt because the aggregates cannot overlap anymore)
						for (unsigned int k = 0; k < SC_UI(ionpVector_oneAgg.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
							if (j == 0 && i == 0) {
								break;//skip this validation for the first ionp
							}
							if (distance(ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].xt(), ionpVector_oneAgg[k].xt()) < (ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].radius() + ionpVector_oneAgg[SC_UI(ionpVector_oneAgg.size() - 1)].coating() + ionpVector_oneAgg[k].radius() + ionpVector_oneAgg[k].coating())) {
								hit_ionp = true;
								break;
							}
						}
						if (!hit_ionp) {
							break;
						}
					}
				}
			}

			ionpVector[i] = ionpVector_oneAgg;
		}

		for (unsigned int i = 0; i <  SC_UI(aggregateInstances.size()); ++i){
			for (unsigned int j = 0; j <  SC_UI(ionpVector[i].size()); ++j){
				ionpInstances.push_back(ionpVector[i][j]);
			}
		}
	}
	else{
		for (unsigned int i = 0; i <  SC_UI(aggregateInstances.size()); ++i) {

		for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
			if ( r_ionp_std == 0. ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else{
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., rand()));
			}

			while (true) { // random placement and validation of position
				bool hit_ionp = false; //flag if proton enters IONP ->update each time enters the loop for repositioning

				ionpInstances[i * SC_UI(ionp_n) + j].placeRandomSphere(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_agg, true);

				for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
					if (j == 0 && i == 0) {
						break;//skip this validation for the first ionp
					}
					if (distance(ionpInstances[i * SC_UI(ionp_n) + j].xt(), ionpInstances[k].xt()) < (ionpInstances[i * SC_UI(ionp_n) + j].radius() + ionpInstances[k].radius())) {
						hit_ionp = true;
						break;
					}
				}
				if (!hit_ionp) {
					break;
				}
			}

		}

		}
	}

	std::cout << "Dense spherical Aggregates created: " << aggregateInstances.size() << std::endl;
	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}
