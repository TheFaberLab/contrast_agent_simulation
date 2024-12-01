/*****************************************************************************
*Class to apply different MR sequences

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef MRSEQUENCE_HPP_
#define MRSEQUENCE_HPP_

class MRsequence
{
  public:
                     MRsequence          ( SimMatrix& simM, int sequence_n, int prot_n, int bound_prot_n, double r_i, double e_time, double time_s, double e_spacing, int numThread, int m_nxnodes, double m_xlength, Point3D& dirIONP, double velIONP, double thickness_IONPcoat, double back_T2, double magM, Point3D& Cubeorig, Point3D& Cubesize, bool BGrid, bool var_Step );

    void             replaceProton       ( std::vector<Proton>& protonInstances, std::vector<Ionp>& ionpInstances, double dt, uint i, uint j, RandomList& RandomNormalX, RandomList& RandomNormalY, RandomList& RandomNormalZ, RandomList& RandomDiff );
    void             replaceIONP         ( std::vector<Ionp>&   ionpInstances );

    void             Sequence            ( std::vector<Proton>& protonInstances, std::vector<Ionp>& ionpInstances_orig, std::vector<std::complex<double>>& signalVector, std::vector<double>& timeV, std::vector<std::complex<double>>& signalVectorMSE, std::vector<double>& timeVMSE);


 private:
    SimMatrix simMatrix;
    int sequence_num, protons_n, bound_protons_n;
    double r_ionp, echo_time, time_step, echo_spacing;
    int numThreads, matrix_nxnodes;
    double matrix_xlength;
    Point3D directionIONP;
    double vIONP, thickness_IONPcoating, background_T2, magnetic_moment;
    Point3D cube_o, cube_s;
    bool Grid, var_Steps;


};

inline MRsequence::MRsequence(SimMatrix& simM, int sequence_n, int prot_n, int bound_prot_n, double r_i, double e_time, double time_s, double e_spacing, int numThread, int m_nxnodes, double m_xlength, Point3D& dirIONP, double velIONP, double thickness_IONPcoat, double back_T2, double magM, Point3D& Cubeorig, Point3D& Cubesize, bool BGrid, bool var_Step):
  simMatrix(simM),
  sequence_num(sequence_n),
  protons_n(prot_n),
  bound_protons_n(bound_prot_n),
  r_ionp(r_i),
  echo_time(e_time),
  time_step(time_s),
  echo_spacing(e_spacing),
  numThreads(numThread),
  matrix_nxnodes(m_nxnodes),
  matrix_xlength(m_xlength),
  directionIONP(dirIONP),
  vIONP(velIONP),
  thickness_IONPcoating(thickness_IONPcoat),
  background_T2(back_T2),
  magnetic_moment(magM),
  cube_o(Cubeorig),
  cube_s(Cubesize),
  Grid(BGrid),
  var_Steps(var_Step){
}

#endif
