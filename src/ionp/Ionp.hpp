/*****************************************************************************
*Class to define an IONP (header .hpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef IONP_HPP_
#define IONP_HPP_

#include "Point3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"

class Ionp : public Particle {
public:
    Ionp(double x = 0., double y = 0., double z = 0., double t = 0., double r = 0., double c = 0.);
    Ionp(const Point3D& p, double t = 0., double r = 0., double c = 0.);

    inline double diffTime() const;
    inline double radius() const;
    inline double coating() const;
    inline void setRadius(double r);
    inline void setCoating(double c);

private:
    inline void setDiffTime(double t);
    void diffuse(double d, double dt);

    double m_diff_time = 0.;
    double m_radius;
    double m_coating;
};

//inline definitions
inline Ionp::Ionp(double x, double y, double z, double t, double r, double c) : Particle(x, y, z, t), m_radius(r), m_coating(c) {} 
inline Ionp::Ionp(const Point3D& p, double t, double r, double c) : Particle(p, t), m_radius(r), m_coating(c) {}



inline double Ionp::diffTime() const {
    return m_diff_time;
}
inline double Ionp::radius() const {
    return m_radius;
}
inline double Ionp::coating() const {
    return m_coating;
}

inline void Ionp::setDiffTime(double t) {
    m_diff_time = t;
}

inline void Ionp::setRadius(double r) {
    m_radius = r;
}

inline void Ionp::setCoating(double c) {
    m_coating = c;
}

#endif
