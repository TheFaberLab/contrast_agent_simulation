/*****************************************************************************
*Class to read configuration files and write output files

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef DATA_HPP_
#define DATA_HPP_

class Data {

public:

                     Data                (std::string filename_in);

    void             readData            (  );
    void             readIONP            ( std::string& path_to_ionp_data, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, double r_ionp );
    void             readAgg             ( std::string& path_to_ionp_data, std::vector<Aggregate>& ionpInstances, SimMatrix& simMatrix, double r_ionp );
    void             readHist            ( std::string& path_to_ionp_data, std::vector<double>& randomHist );

    void             writeData           ( std::vector<std::complex<double>>& Signal_vec, std::vector<double>& timeV );
    void             writeData           ( std::vector<double>& Phase_vec );
    void             writeData           ( std::vector<Ionp>& IONP_vec );
    void             writeData           ( std::vector<Aggregate>& Agg_vec );
    void             writeData           ( std::vector<Proton>& Proton_vec );
    void             log                 ( bool reset );

    void             init_logger         ( void );


    double get_matrix_x0();
    double get_matrix_y0();
    double get_matrix_z0();
    double get_matrix_xlength();
    double get_matrix_ylength();
    double get_matrix_zlength();
    int get_matrix_nxnodes();
    int get_matrix_nynodes();
    int get_matrix_nznodes();
    double get_radius_ionp();
    double get_radius_agg();
    int get_ionp_number();
    int get_agg_number();
    int get_proton_number();
    double get_echo_time();
    double get_time_step();
    double get_echo_spacing();
    int get_numThreads();
    int get_sequence_num();
    double get_velocity_IONP();
    Point3D get_direction_IONP();
    double get_thickness_IONPcoating();
    double get_thickness_IONPcoating_std();
    int get_bound_proton_number();
    double get_background_T2();
    double get_radius_ionp_1();
    double get_radius_ionp_std();
    double get_radius_ionp_std_1();
    double get_radius_agg_std();
    double get_concentration();
    double get_mag_moment();
    double get_ionp_percentage();
    Point3D get_Cube_Origin();
    Point3D get_Cube_size();
    bool get_DLS_file();

    void set_time_step(double dt);


private:
  std::string filename;
  double matrix_x0, matrix_y0, matrix_z0; //matrix
  double matrix_xlength, matrix_ylength, matrix_zlength;
  int matrix_nxnodes, matrix_nynodes, matrix_nznodes;
  double radius_ionp, radius_agg; //particles
  int ionp_number, agg_number;
  int proton_number, logID, numThreads;
  double echo_time, time_step, echo_spacing;
  int sequence_num;
  double velocityIONP;
  Point3D directionIONP;
  double thickness_IONPcoating, thickness_IONPcoating_std;
  int bound_proton_number;
  double background_T2;
  double radius_ionp_1, radius_ionp_std, radius_ionp_std_1, radius_agg_std, concentration, mag_moment, ionp_percentage;
  Point3D Cube_Origin, Cube_size;
  bool DLS_file;

};

inline Data::Data(std::string filename_in) : filename(filename_in), matrix_x0(), matrix_y0(), matrix_z0(), matrix_xlength(), matrix_ylength(), matrix_zlength(), matrix_nxnodes(), matrix_nynodes(), matrix_nznodes(), radius_ionp(), radius_agg(), ionp_number(), agg_number(), proton_number(), logID(0), numThreads(), echo_time(), time_step(), echo_spacing(), sequence_num(0), velocityIONP(), directionIONP(), thickness_IONPcoating(), thickness_IONPcoating_std(), bound_proton_number(), background_T2(), radius_ionp_1(), radius_ionp_std(), radius_ionp_std_1(), radius_agg_std(), concentration(), mag_moment(), ionp_percentage(), Cube_Origin(), Cube_size(), DLS_file() {}


#endif
