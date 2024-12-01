/*****************************************************************************
*Class to read configuration files and write output files

*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <complex>
#include <filesystem>

#include "mptmacros.h"
#include "Point3D.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"

//Function to read all defined parameters from ini file in /ini/
void Data::readData (  )
{
  std::cout << "Entering readFunc...." << std::endl;
	std::string line_in, name, number_str;


	std::ifstream file_in(filename, std::ios::in);

	if (!file_in.is_open()) {
		std::cerr << "Unable to open input file: " << filename << std::endl;
		exit(1);
	}

	while (std::getline(file_in, line_in, '\n')) {

    std::size_t pos = line_in.find("=");

    number_str = line_in.substr (pos+1);

    if (!isspace(line_in[pos-1])){
      name = line_in.substr (0,pos);
    }
    else{
      name = line_in.substr (0,pos-1);
    }

    if (name == "IONP_r" || name == "IONP_radius" || name == "RIONP" || name == "R_IONP" || name == "radiusIONP"){
      radius_ionp = std::stod(number_str);
    }
    else if (name == "Agg_r" || name == "Agg_radius" || name == "Aggregate_radius" || name == "Aggregate_r" || name == "R_Agg" || name == "radiusAgg"){
      radius_agg = std::stod(number_str);
    }
    else if (name == "Matrix_x0" || name == "matrix_x0"){
      matrix_x0 = std::stod(number_str);
    }
    else if (name == "Matrix_y0" || name == "matrix_y0"){
      matrix_y0 = std::stod(number_str);
    }
    else if (name == "Matrix_z0" || name == "matrix_z0"){
      matrix_z0 = std::stod(number_str);
    }
    else if (name == "Origin"){
      matrix_x0 = std::stod(number_str);
      matrix_y0 = std::stod(number_str);
      matrix_z0 = std::stod(number_str);
    }
    else if (name == "Matrix_xlength" || name == "matrix_xlength" || name == "xlength" || name == "xLength"){
      matrix_xlength = std::stod(number_str);
    }
    else if (name == "Matrix_ylength" || name == "matrix_ylength" || name == "ylength" || name == "yLength"){
      matrix_ylength = std::stod(number_str);
    }
    else if (name == "Matrix_zlength" || name == "matrix_zlength" || name == "zlength" || name == "zLength"){
      matrix_zlength = std::stod(number_str);
    }
    else if (name == "Length"){
      matrix_xlength = std::stod(number_str);
      matrix_ylength = std::stod(number_str);
      matrix_zlength = std::stod(number_str);
    }
    else if (name == "Matrix_xnodes" || name == "matrix_xnodes" || name == "xnodes" || name == "xNodes"){
      matrix_nxnodes = std::stoi(number_str);
    }
    else if (name == "Matrix_ynodes" || name == "matrix_ynodes" || name == "ynodes" || name == "yNodes"){
      matrix_nynodes = std::stoi(number_str);
    }
    else if (name == "Matrix_znodes" || name == "matrix_znodes" || name == "znodes" || name == "zNodes"){
      matrix_nznodes = std::stoi(number_str);
    }
    else if (name == "Nodes"){
      matrix_nxnodes = std::stoi(number_str);
      matrix_nynodes = std::stoi(number_str);
      matrix_nznodes = std::stoi(number_str);
    }
    else if (name == "IONP_n" || name == "IONP_num" || name == "IONP_number" || name == "ionp_n" || name == "ionp_num" || name == "ionp_number"){
      ionp_number = std::stoi(number_str);
    }
    else if (name == "Agg_n" || name == "Agg_num" || name == "Agg_number" || name == "agg_n" || name == "agg_num" || name == "agg_number"){
      agg_number = std::stoi(number_str);
    }
    else if (name == "Prot_n" || name == "Prot_num" || name == "Prot_number" || name == "Proton_n" || name == "Proton_num" || name == "Proton_number"){
      proton_number = std::stoi(number_str);
    }
    else if (name == "echo_time" || name == "T_E" || name == "T_Echo" || name == "T_echo" || name == "Echo_time" || name == "Echo_Time"){
      echo_time = std::stod(number_str);
    }
    else if (name == "time_step" || name == "dt" || name == "DeltaT" || name == "Time_Step" || name == "Time_step"){
      time_step = std::stod(number_str);
    }
    else if (name == "echo_spacing" || name == "Echo_spacing" || name == "SE_time" || name == "SE_Time"){
      echo_spacing = std::stod(number_str);
    }
    else if (name == "Threads" || name == "numThreads" || name == "num_Threads" || name == "number_Threads" || name == "Number_Threads"){
      numThreads = std::stoi(number_str);
    }
    else if (name == "vIONP" || name == "velocityIONP" || name == "VelocityIONP" || name == "Velocity_IONP" || name == "velocity_IONP"){
      velocityIONP = std::stod(number_str);
    }
    else if (name == "Direction" || name == "direction" || name == "MoveDirection" || name == "MovementDirection" || name == "Move Direction" || name == "Movement direction"){
      std::size_t elementsXY = number_str.find(" ");
      if (elementsXY == 0){
        elementsXY = number_str.find(" ", elementsXY+1);
      }
      std::size_t elementsYZ = number_str.find(" ", elementsXY+1);
      std::string x_str = number_str.substr (0, elementsXY);
      std::string y_str = number_str.substr (elementsXY+1, elementsYZ-elementsXY);
      std::string z_str = number_str.substr (elementsYZ+1);
      double x = std::stod(x_str);
      double y = std::stod(y_str);
      double z = std::stod(z_str);
      directionIONP.set(x, y, z);
    }
    else if (name == "Sequence" || name == "MR sequence" || name == "MRsequence" ){
      if (number_str == " FID" || number_str == "FID"){
        sequence_num = 0;
        std::cout << "Sequence FID is used" << std::endl;
      }
      else if (number_str == " SE" || number_str == " Spin Echo" || number_str == " spin echo" || number_str == "SE" || number_str == "Spin Echo" || number_str == "spin echo"){
        sequence_num = 1;
        std::cout << "Sequence Spin Echo is used" << std::endl;
      }
      else if (number_str == " MSE" || number_str == " Multi Spin Echo" || number_str == " multi spin echo" || number_str == "MSE" || number_str == "Multi Spin Echo" || number_str == "multi spin echo"){
        sequence_num = 2;
        std::cout << "Sequence Multi Spin Echo is used" << std::endl;
      }
      else{
        std::cerr << "Unknown sequence: " << number_str << std::endl;
  		exit(1);
      }
    }
    else if (name == "DLS_file" || name == "DLS_File" || name == "DLS file" || name == "DLS File" || name == "dls file"  || name == "dls_file" ){
      if (number_str == "Yes" || number_str == " Yes" || number_str == "yes" || number_str == " yes" || number_str == "ja" || number_str == " ja" || number_str == "Ja" || number_str == " Ja"){
        DLS_file = true;
      }
      else if (number_str == "No" || number_str == " No" || number_str == "no" || number_str == " no" || number_str == "nein" || number_str == " nein" || number_str == "Nein" || number_str == " Nein"){
        DLS_file = false;
      }
      else{
        std::cerr << "Neither Yes or No: " << number_str << std::endl;
  		exit(1);
      }
    }
    else if (name == "Coating" || name == "Shell" || name == "Coating_Thickness" || name == "Coating Thickness" || name == "Coating_thickness"){
      thickness_IONPcoating = std::stod(number_str);
    }
    else if (name == "Coating_std" || name == "Shell_std" || name == "Coating_Thickness_std" || name == "Coating_Std" || name == "Shell_Std" || name == "Coating_Thickness_Std" || name == "Coating Thickness Std" || name == "Coating_thickness_std"){
      thickness_IONPcoating_std = std::stod(number_str);
    }
    else if (name == "Bound_Prot_n" || name == "Bound_Prot_num" || name == "Bound_Prot_number" || name == "Bound_Proton_n" || name == "Bound_Proton_num" || name == "Bound_Proton_number"){
      bound_proton_number = std::stoi(number_str);
    }
    else if (name == "T2" || name == "Background T2" || name == "Background T_2" || name == "T2_out"){
      background_T2 = std::stod(number_str);
    }
    else if (name == "IONP_r_std" || name == "IONP_radius_std" || name == "RIONP_std" || name == "R_IONP_std" || name == "radiusIONP_std"){
      radius_ionp_std = std::stod(number_str);
    }
    else if (name == "IONP_r_1" || name == "IONP_radius_1" || name == "RIONP_1" || name == "R_IONP_1" || name == "radiusIONP_1"){
      radius_ionp_1 = std::stod(number_str);
    }
    else if (name == "IONP_r_std_1" || name == "IONP_radius_std_1" || name == "RIONP_std_1" || name == "R_IONP_std_1" || name == "radiusIONP_std_1"){
      radius_ionp_std_1 = std::stod(number_str);
    }
    else if (name == "Agg_r_std" || name == "Agg_radius_std" || name == "Aggregate_radius_std" || name == "Aggregate_r_std" || name == "R_Agg_std" || name == "radiusAgg_std"){
      radius_agg_std = std::stod(number_str);
    }
    else if (name == "Concentration" || name == "Conc"){
      concentration = std::stod(number_str);
    }
    else if (name == "Magnetic Moment" || name == "Magnetic_Moment" || name == "Magnetic moment" || name == "Magnetic_moment" || name == "magnetic moment" || name == "magnetic_moment"){
      mag_moment = std::stod(number_str);
    }
    else if (name == "IONP_percentage" || name == "IONP_per" || name == "IONP_Percentage" || name == "IONP_Per" || name == "IONP percentage" || name == "IONP Percentage"){
      ionp_percentage = std::stod(number_str);
    }
    else if (name == "Cube Origin" || name == "Cube_Origin" || name == "Cube_Orig" || name == "CubeOrigin" || name == "CubeOrig" || name == "Cube_origin"){

      std::size_t spaceBefore = 0;
      std::size_t spaceAfter = 0;

      std::vector<double> vec(3);

      for (unsigned int i = 0; i < 3; ++i) {
        spaceAfter = number_str.find(" ", spaceBefore+1);

        std::string x_str = number_str.substr (spaceBefore+1, spaceAfter-spaceBefore);
        vec[i] = std::stod(x_str);

        spaceBefore = spaceAfter;

      }
      Cube_Origin.set(vec[0], vec[1], vec[2]);
    }
    else if (name == "Cube Size" || name == "Cube_Size" || name == "Cube_size" || name == "CubeSize"){

      std::size_t spaceBefore = 0;
      std::size_t spaceAfter = 0;

      std::vector<double> vec(3);

      for (unsigned int i = 0; i < 3; ++i) {
        spaceAfter = number_str.find(" ", spaceBefore+1);

        std::string x_str = number_str.substr (spaceBefore+1, spaceAfter-spaceBefore);
        vec[i] = std::stod(x_str);

        spaceBefore = spaceAfter;

      }
      Cube_size.set(vec[0], vec[1], vec[2]);
    }
    else {
      std::cerr << "A line in the config file has an incorrect expression: " << name << std::endl;
  		exit(1);
    }


	}
	file_in.close();

  if (DLS_file == false){
    std::cout<<"\n"<<"IONP radius: "<< radius_ionp <<std::endl;
    if (radius_ionp_std != 0.){
        std::cout<<"IONP radius std (var. radius): "<< radius_ionp_std <<std::endl;
    }
  }
  else{
    std::cout << "An IONP radius is used drawn randomly from file in ../../../ini/DLS_File" << "\n" << std::endl;
  }

  if (thickness_IONPcoating != 0.){
    std::cout<<"IONP coating thickness: "<< thickness_IONPcoating <<std::endl;
  }
  if (thickness_IONPcoating_std != 0.){
    std::cout<<"IONP coating thickness std (var. coating): "<< thickness_IONPcoating_std <<std::endl;
  }
  if (agg_number != 0. || radius_agg != 0.)
  {
    std::cout<<"Aggregate radius: "<< radius_agg <<std::endl;
    std::cout<<"Aggregate number: "<< agg_number <<std::endl;
    if (radius_agg_std != 0.){
      std::cout<<"Aggregate radius std (var. radius): "<< radius_agg_std <<std::endl;
    }
  }
  std::cout<<"Proton number: "<< proton_number <<std::endl;
  if (bound_proton_number != 0.){
    std::cout<<"Bound Proton number: "<< bound_proton_number << std::endl;
  }
  std::cout << std::endl;

  if ( mag_moment != 0. ){
    std::cout<<"Magnetic Moment: "<< mag_moment << std::endl;
  }

  if ( radius_ionp_1 != 0. && radius_ionp != 0. ){
    if ( ionp_percentage == 0. ){
      ionp_percentage = 70.;
    }
    std::cout<<"Bimodal radii distribution" << std::endl;
    std::cout<<"Original radius " << ionp_percentage << "% probability and second radius " << (100. - ionp_percentage) << "%" << std::endl;
    std::cout<<"IONP radius 1: " << radius_ionp_1 << std::endl;
    std::cout<<"IONP radius 1 std: " << radius_ionp_std_1 << "\n" << std::endl;
    ionp_percentage *= 0.01;
  }

  if (concentration != 0.){
    std::cout<<"Concentration: "<< concentration <<"µmol/l\n" << std::endl;
  }
  else{
    std::cout<<"IONP number: "<< ionp_number <<std::endl;
  }

  if (background_T2 != 0.){
    std::cout<<"Set T2 of background: "<< background_T2  << std::endl;
  }
  std::cout << std::endl;

  std::cout<<"Origin: "<<matrix_x0<<" "<<matrix_y0<<" "<<matrix_z0<<std::endl;
  std::cout<<"Nodes: "<<matrix_nxnodes<<" "<<matrix_nynodes<<" "<<matrix_nznodes<<std::endl;
  std::cout<<"Dimensions given: "<<matrix_xlength<<" "<<matrix_ylength<<" "<<matrix_zlength<<std::endl;


  if( matrix_xlength == 0. || matrix_ylength == 0. || matrix_zlength == 0.){
    std::cout<<"No Simulation Volume size given ... using constant volume fraction "<<std::endl;
    std::cout<<"Volume fraction: "<<D_VOL_FRAC<<std::endl;

    double mtx_dimensions;

    if (agg_number == 0. || radius_agg == 0.){
      mtx_dimensions = pow(ionp_number*(4./3.)*(D_PI)*(1/D_VOL_FRAC), 1./3.)*radius_ionp;
    }
    else {
      mtx_dimensions = pow(ionp_number*agg_number*(4./3.)*(D_PI)*(1/D_VOL_FRAC), 1./3.)*radius_ionp;
    }

    std::cout<<"Dimensions given: "<<mtx_dimensions<<"\n" <<std::endl;

    matrix_xlength = mtx_dimensions;
    matrix_ylength = mtx_dimensions;
    matrix_zlength = mtx_dimensions;
  }

  if ( Cube_size.x() != 0. ){
    std::cout << "Cube Origin: " << Cube_Origin.x() << " " << Cube_Origin.y() << " " << Cube_Origin.z() << std::endl;
    std::cout << "Cube Size: " << Cube_size.x() << " " << Cube_size.y() << " " << Cube_size.z() << "\n\n";
  }

  if (numThreads == 0){
    std::cout<<"Set threads: "<< numThreads <<std::endl;
    std::cout<<"Maximum available threads are used!!" <<std::endl;
  }
  if (velocityIONP != 0){
    std::cout << "The IONPs are moving during the Simulation!" << std::endl;
    std::cout << "Velocity: " << velocityIONP << "m/s" << std::endl;
    std::cout << "Direction: " << directionIONP.x() << " " << directionIONP.y() << " " << directionIONP.z() << std::endl;
  }

  std::cout<<"Echo Time: "<< echo_time <<std::endl;
  std::cout<<"Time Step: "<< time_step <<std::endl;
  if (echo_spacing != 0.){
    std::cout<<"Echo Spacing (SE/MSE): "<< echo_spacing <<std::endl;
  }

}

//Function to read an Histogram File that contains a radii distribution (often DLS Files)
void Data::readHist ( std::string& path_to_ionp_data, std::vector<double>& randomHist)
{
  std::cout << "Beginning to read Histogram/DLS File" << std::endl;


  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		std::cerr << "Unable to open input file: " << path_to_ionp_data << std::endl;
		exit(1);
	}

  std::vector<double> bins;
  std::vector<int> frequencies;

  double bin;
  double frequency;

  while(!file_in.eof())
  {
    file_in >> bin >> frequency;

    bins.push_back(bin);
    frequencies.push_back(SC_I(frequency));

  }

  RandNumberHist rand(bins, frequencies);

  int lenVector = 10000;

  randomHist = rand.randomVector(lenVector);

  writeData(randomHist);


  std::cout << "Random Vector according to Histogram: " << randomHist.size() << std::endl;

}

//Function to read IONP position files
void Data::readIONP ( std::string& path_to_ionp_data, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, double r_ionp )
{
  std::cout << "Beginning to read IONP position file" << std::endl;

  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		std::cerr << "Unable to open input file: " << path_to_ionp_data << std::endl;
		exit(1);
	}

  while(!file_in.eof())
  {
    double x, y, z;
    file_in >> x >> y >> z;

    ionpInstances.push_back(Ionp(x,y,z, 0., r_ionp));

  }

  for (unsigned int i = 0; i < SC_UI(ionpInstances.size()); ++i) {
		for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
      if (i != k){
        if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (2. * r_ionp)) {
					std::cerr << "IONPs are placed too close to each other: " << i << " and " << k << std::endl;
		      exit(1);
				}
			}
    }
    if (abs(ionpInstances[i].xt().x()) > simMatrix.xLength() / 2.) {
      std::cerr << "IONP " << i  <<" is outside of Simulation Volume in x direction: " << ionpInstances[i].xt().x() << " vs. " << simMatrix.xLength() / 2. << std::endl;
		  exit(1);
		}
		if (abs(ionpInstances[i].xt().y()) > simMatrix.yLength() / 2.) {
			std::cerr << "IONP " << i  <<" is outside of Simulation Volume in y direction: " << ionpInstances[i].xt().y() << " vs. " << simMatrix.yLength() / 2. << std::endl;
		  exit(1);
		}
		if (abs(ionpInstances[i].xt().z()) > simMatrix.zLength() / 2.) {
			std::cerr << "IONP " << i  <<" is outside of Simulation Volume in z direction: " << ionpInstances[i].xt().z() << " vs. " << simMatrix.zLength() / 2. << std::endl;
		  exit(1);
		}

    std::cout << ionpInstances[i].xt().x() << " " << ionpInstances[i].xt().y() << " " << ionpInstances[i].xt().z() << std::endl;
  }

  std::cout << "IONPs created: " << ionpInstances.size() << std::endl;

}

//Function to read Aggregate position files
void Data::readAgg ( std::string& path_to_ionp_data, std::vector<Aggregate>& ionpInstances, SimMatrix& simMatrix, double r_ionp )
{
  std::cout << "Beginning to read Agg position file" << std::endl;

  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		std::cerr << "Unable to open input file: " << path_to_ionp_data << std::endl;
		exit(1);
	}

  while(!file_in.eof())
  {
    double x, y, z;
    file_in >> x >> y >> z;

    ionpInstances.push_back(Aggregate(x,y,z, 0., r_ionp, 0.));

  }

  for (unsigned int i = 0; i < SC_UI(ionpInstances.size()); ++i) {
		for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
      if (i != k){
        if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (2. * r_ionp)) {
					std::cerr << "Aggregates are placed too close to each other: " << i << " and " << k << std::endl;
		      exit(1);
				}
			}
    }
    if (abs(ionpInstances[i].xt().x()) > simMatrix.xLength() / 2.) {
      std::cerr << "Aggregate " << i  <<" is outside of Simulation Volume in x direction: " << ionpInstances[i].xt().x() << " vs. " << simMatrix.xLength() / 2. << std::endl;
		  exit(1);
		}
		if (abs(ionpInstances[i].xt().y()) > simMatrix.yLength() / 2.) {
			std::cerr << "Aggregate " << i  <<" is outside of Simulation Volume in y direction: " << ionpInstances[i].xt().y() << " vs. " << simMatrix.yLength() / 2. << std::endl;
		  exit(1);
		}
		if (abs(ionpInstances[i].xt().z()) > simMatrix.zLength() / 2.) {
			std::cerr << "Aggregate " << i  <<" is outside of Simulation Volume in z direction: " << ionpInstances[i].xt().z() << " vs. " << simMatrix.zLength() / 2. << std::endl;
		  exit(1);
		}
  }

}

//Function to write signal output file
void Data::writeData ( std::vector<std::complex<double>>& Signal_vec, std::vector<double>& timeV )
{

  std::string path_to_data = "../../../data/tmp/";
  int id = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  std::string magnetization_out_name;

  if ( sequence_num == 0 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_signal_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_signal_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_signal_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

  std::cout << magnetization_out_name << std::endl;

	std::ofstream outfile(magnetization_out_name, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		std::cerr << "Unable to open output file: signal" << std::endl;
		exit(1);
	}
	for (unsigned int i = 0; i < Signal_vec.size(); ++i) {
		outfile << real(Signal_vec[i]) << "\t" << timeV[i] << std::endl;
	}
	outfile.close();

}

//Function to write optional phase output file (phase of every proton at every time)
void Data::writeData ( std::vector<double>& Phase_vec )
{

  std::string path_to_data = "../../../data/tmp/";
  int id = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  std::string magnetization_out_name;

  if ( sequence_num == 0 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_Phase_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_Phase_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    magnetization_out_name = path_to_data +  std::to_string( id ) + "_Phase_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

  std::cout << magnetization_out_name << std::endl;

	std::ofstream outfile(magnetization_out_name, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		std::cerr << "Unable to open output file: phase" << std::endl;
		exit(1);
	}
	for (unsigned int i = 0; i < Phase_vec.size(); ++i) {
		outfile << Phase_vec[i] << std::endl;
	}
	outfile.close();

}

//Function to write IONP position output file
void Data::writeData ( std::vector<Ionp>& IONP_vec )
{

  std::string path_to_data = "../../../data/tmp/";
  int id = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  std::string ionps_out_name;

  if ( sequence_num == 0 ){
    ionps_out_name = path_to_data +  std::to_string( id ) + "_IONPs_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    ionps_out_name = path_to_data +  std::to_string( id ) + "_IONPs_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    ionps_out_name = path_to_data +  std::to_string( id ) + "_IONPs_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

  std::cout << ionps_out_name << std::endl;

	std::ofstream outfile(ionps_out_name, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		std::cerr << "Unable to open output file: IONPs" << std::endl;
		exit(1);
	}
  for (unsigned int i = 0; i < IONP_vec.size(); ++i) {
    outfile << IONP_vec[i].xt().x()<<"\t"<<IONP_vec[i].xt().y()<<"\t"<<IONP_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

//Function to write Aggregate position output file
void Data::writeData ( std::vector<Aggregate>& Agg_vec )
{

  std::string path_to_data = "../../../data/tmp/";
  int id = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  std::string agg_out_name;

  if ( sequence_num == 0 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Aggs_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Aggs_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Aggs_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

  std::cout << agg_out_name << std::endl;

	std::ofstream outfile(agg_out_name, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		std::cerr << "Unable to open output file: Aggregates" << std::endl;
		exit(1);
	}
  for (unsigned int i = 0; i < Agg_vec.size(); ++i) {
    outfile << Agg_vec[i].xt().x()<<"\t"<<Agg_vec[i].xt().y()<<"\t"<<Agg_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

//Function to write optional Proton position output file
void Data::writeData ( std::vector<Proton>& Proton_vec )
{

  std::string path_to_data = "../../../data/tmp/";
  int id = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  std::string agg_out_name;

  if ( sequence_num == 0 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Protons_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Protons_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    agg_out_name = path_to_data +  std::to_string( id ) + "_Protons_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

  std::cout << agg_out_name << std::endl;

	std::ofstream outfile(agg_out_name, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		std::cerr << "Unable to open output file: Aggregates" << std::endl;
		exit(1);
	}
  for (unsigned int i = 0; i < Proton_vec.size(); ++i) {
    outfile << Proton_vec[i].xt().x()<<"\t"<<Proton_vec[i].xt().y()<<"\t"<<Proton_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

//Function to write log file
void Data::log ( bool reset )
{

  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf

  if (reset == false){
    std::string path_to_data = "../../../data/tmp/";
  std::string path_to_log = "../../../data/log/";
  int id = 0;
  int logid = 0;

	for (const auto & entry : std::filesystem::directory_iterator(path_to_data)){
    (void) entry;
    id++;
  }

  for (const auto & entry : std::filesystem::directory_iterator(path_to_log)){
    (void) entry;
    logid++;
  }

  std::string log_out_name = path_to_log + std::to_string( logid ) + "_" + std::to_string( id ) + "_log.txt";

  std::ofstream out(log_out_name);
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to log.txt!
  }
  else{
    std::cout.rdbuf(coutbuf); //reset to standard output again
  }


}


//Get functions of all defined private variables
double Data::get_matrix_x0(){
  return matrix_x0;
}
double Data::get_matrix_y0(){
  return matrix_y0;
}
double Data::get_matrix_z0(){
  return matrix_z0;
}
double Data::get_matrix_xlength(){
  return matrix_xlength;
}
double Data::get_matrix_ylength(){
  return matrix_ylength;
}
double Data::get_matrix_zlength(){
  return matrix_zlength;
}
int Data::get_matrix_nxnodes(){
  return matrix_nxnodes;
}
int Data::get_matrix_nynodes(){
  return matrix_nynodes;
}
int Data::get_matrix_nznodes(){
  return matrix_nznodes;
}
double Data::get_radius_ionp(){
  return radius_ionp;
}
double Data::get_radius_agg(){
  return radius_agg;
}
int Data::get_ionp_number(){
  return ionp_number;
}
int Data::get_agg_number(){
  return agg_number;
}
int Data::get_proton_number(){
  return proton_number;
}
double Data::get_echo_time(){
  return echo_time;
}
double Data::get_time_step(){
  return time_step;
}
double Data::get_echo_spacing(){
  return echo_spacing;
}
int Data::get_numThreads(){
  return numThreads;
}
int Data::get_sequence_num(){
  return sequence_num;
}
double Data::get_velocity_IONP(){
  return velocityIONP;
}
Point3D Data::get_direction_IONP(){
  return directionIONP;
}
double Data::get_thickness_IONPcoating(){
  return thickness_IONPcoating;
}
double Data::get_thickness_IONPcoating_std(){
  return thickness_IONPcoating_std;
}
int Data::get_bound_proton_number(){
  return bound_proton_number;
}
double Data::get_background_T2(){
  return background_T2;
}
double Data::get_radius_ionp_1(){
  return radius_ionp_1;
}
double Data::get_radius_ionp_std(){
  return radius_ionp_std;
}
double Data::get_radius_ionp_std_1(){
  return radius_ionp_std_1;
}
double Data::get_radius_agg_std(){
  return radius_agg_std;
}
double Data::get_concentration(){
  return concentration;
}
double Data::get_mag_moment(){
  return mag_moment;
}
double Data::get_ionp_percentage(){
  return ionp_percentage;
}
Point3D Data::get_Cube_Origin(){
  return Cube_Origin;
}
Point3D Data::get_Cube_size(){
  return Cube_size;
}
bool Data::get_DLS_file(){
  return DLS_file;
}

void Data::set_time_step(double dt){
  time_step = dt;
}
