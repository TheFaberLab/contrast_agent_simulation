/*****************************************************************************
*Class to define an aggregate (header .hpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef AGGREGATE_HPP_
#define AGGREGATE_HPP_

#include "Point3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"

class Aggregate : public Particle{
public:
    Aggregate(double x = 0., double y = 0., double z = 0., double t=0.,  double r = 0, const int n = 1);
    Aggregate(const Point3D& p, double t = 0., double r = 0., int n = 1);
    Aggregate(double r = 0., int n = 1);
    Aggregate(int n = 1);

    inline double diffTime() const;
    inline double radius() const;
    inline int nparticles() const;

private:
    inline void setDiffTime(double t);
    void diffuse(double d, double dt);

    double m_diff_time = 0;
    double m_radius;
    int m_nparticles;
};

//inline definitions
inline Aggregate::Aggregate(double x, double y, double z, double t, double r, int n) : Particle(x,y,z,t), m_radius(r), m_nparticles(n)  {}
inline Aggregate::Aggregate(const Point3D& p, double t, double r, int n) : Particle(p, t), m_radius(r), m_nparticles(n) {}
inline Aggregate::Aggregate(double r, int n) : Particle(0.,0.,0.,0.), m_radius(r), m_nparticles(n) {}
inline Aggregate::Aggregate(int n) : Particle(0.,0.,0.,0.), m_radius(), m_nparticles(n) {}

inline double Aggregate::diffTime() const {
    return m_diff_time;
}
inline double Aggregate::radius() const {
    return m_radius;
}
inline int Aggregate::nparticles() const {
    return m_nparticles;
}

inline void Aggregate::setDiffTime(double t) {
    m_diff_time = t;
}


#endif
