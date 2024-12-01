/*****************************************************************************
*Class to place IONP in SimVolume in different spacial distributions

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef SIMSPACE_HPP_
#define SIMSPACE_HPP_

namespace SimSpace {

    void             simulationSpace         (std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int ionp_n, double conc, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double thickness_coating, double thickness_ionp_std, double percentage_ionp, bool DLS_file );
    void             surroundingSpace        ( std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int ionp_n, double r_ionp );
    void             simulationSpace         ( std::vector<Aggregate>& aggregateInstances, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int agg_n, int ionp_n, double conc, double r_agg, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double r_agg_std, double thickness_coating, double thickness_ionp_std, double percentage_ionp );
    void             simulationSpace         ( std::vector<Aggregate>& aggregateInstances, std::vector<Ionp>& ionpInstances, SimMatrix& simMatrix, int agg_n, int ionp_n, double conc, double r_ionp, double r_ionp_1, double r_ionp_std, double r_ionp_std_1, double thickness_coating, double percentage_ionp );

}


#endif
