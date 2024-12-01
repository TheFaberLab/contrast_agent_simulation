/*****************************************************************************
*Class to define the simulation matrix (header .hpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef SIMMATRIX_HPP_
#define SIMMATRIX_HPP_

#include <vector>
#include "Point3D.hpp"

class SimMatrix {
public:
    SimMatrix(const Point3D& p, double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double x0, double y0, double z0, double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double l, int n);

    SimMatrix(const SimMatrix& m); //copy
    ~SimMatrix();

    inline Point3D x0() const;
    inline double xLength() const;
    inline double yLength() const;
    inline double zLength() const;
    inline int nxNodes() const;
    inline int nyNodes() const;
    inline int nzNodes() const;


   inline void addBValue(double b);
   inline void addPosValue(const Point3D& p);

   double calculateBfield(const Point3D& p1, const Point3D& p2, const double r); //calculate B between 2 points
   double calculateBfield(const Point3D& p1, const Point3D& p2, const double r, const double magM); //calculate B between 2 points

   void calculateBfield(std::vector< std::vector< std::vector<double> > >& BfieldGrid, std::vector<Ionp>& ionpInstances, int numThreads, double magnetic_moment);

   double interpolate(const Point3D& p1, std::vector< std::vector< std::vector<double> > >& BfieldGrid);


private:
    Point3D m_x0; //origin
    const double m_x_length; //fixed
    const double m_y_length;
    const double m_z_length;
    const int m_nx_nodes;
    const int m_ny_nodes;
    const int m_nz_nodes;

public:
    std::vector<double> b_field;
    std::vector<int>::iterator b_field_iter;
    std::vector<Point3D> pos_nodes;
    std::vector<int>::iterator Iter;
};

//inline definitions
inline SimMatrix::SimMatrix(const Point3D& p, double x_l, double y_l, double z_l, int nx, int ny, int nz) // most general (point)
    : m_x0(p), m_x_length(x_l), m_y_length(y_l), m_z_length(z_l), m_nx_nodes(nx), m_ny_nodes(ny), m_nz_nodes(nz) {}
inline SimMatrix::SimMatrix(double x0, double y0, double z0, double x_l, double y_l, double z_l, int nx, int ny, int nz) // most general (exactly the numbers)
    : m_x0(x0,y0,z0), m_x_length(x_l), m_y_length(y_l), m_z_length(z_l), m_nx_nodes(nx), m_ny_nodes(ny), m_nz_nodes(nz) {}
inline SimMatrix::SimMatrix(double x_l, double y_l, double z_l, int nx, int ny, int nz) //origin(0,0,0), different lengths, nNodes
    : m_x0(0., 0., 0.), m_x_length(x_l), m_y_length(y_l), m_z_length(z_l), m_nx_nodes(nx), m_ny_nodes(ny), m_nz_nodes(nz) {}
inline SimMatrix::SimMatrix(double l, int n) // origin(0,0,0), same length, nNodes
    : m_x0(0., 0., 0.), m_x_length(l), m_y_length(l), m_z_length(l), m_nx_nodes(n), m_ny_nodes(n), m_nz_nodes(n) {}

inline SimMatrix::SimMatrix(const SimMatrix& m)
    : m_x0(m.m_x0), m_x_length(m.m_x_length), m_y_length(m.m_y_length), m_z_length(m.m_z_length), m_nx_nodes(m.m_nx_nodes), m_ny_nodes(m.m_ny_nodes), m_nz_nodes(m.m_nz_nodes) {}
inline SimMatrix::~SimMatrix() {}


inline Point3D SimMatrix::x0() const {
    return m_x0;
}
inline double SimMatrix::xLength() const {
    return m_x_length;
}
inline double SimMatrix::yLength() const {
    return m_y_length;
}
inline double SimMatrix::zLength() const {
    return m_z_length;
}
inline int SimMatrix::nxNodes() const {
    return m_nx_nodes;
}
inline int SimMatrix::nyNodes() const {
    return m_ny_nodes;
}
inline int SimMatrix::nzNodes() const {
    return m_nz_nodes;
}


inline void SimMatrix::addBValue(double b) {
   b_field.push_back(b);
}
inline void SimMatrix::addPosValue(const Point3D& p) {
   pos_nodes.push_back(p);
}


#endif
