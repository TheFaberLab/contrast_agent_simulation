/*****************************************************************************
*Class to define a particle (header .hpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <cmath>
#include <vector>
#include "Point3D.hpp"
#include "RandNumberBetween.hpp"

class Particle {
public:
    Particle(double x = 0., double y = 0., double z = 0., double t = 0);
    Particle(const Point3D& p, double t = 0);

    Particle(const Particle& p); //copy
    ~Particle();

    inline Point3D x0() const;
    inline Point3D xt() const;
    inline double diffTime() const;

    inline void set(const Point3D& p0, const Point3D& pt, double t);
    inline void set(double x0, double y0, double z0, double xt, double yt, double zt, double t);
    inline void setX0(const Point3D& p);
    inline void setX0(double x, double y, double z);
    inline void setXt(const Point3D& p);
    inline void setXt(double x, double y, double z);
    inline void setDiffTime(double t);

    void placeRandom(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high, bool init = false);
    void placeRandom(double low, double high, bool init = false);
    void placeRandom(double x, double y, double z, bool init = false);
    void placeRandom(Point3D& center, double rmin, double rmax, bool init = false);
    void placeRandomSphere(double x0, double y0, double z0, double r_min, double r_max, bool init = false);
    void diffuse(double d, double dt, double theta_rand, double phi_rand);
    void diffuse(double x, double y, double z);
    void move(Point3D& direction, double velocity, double dt);

    Particle& assign(const Particle& p); //self-act operators
    Particle& operator=(const Particle& p);

private:
    std::vector<Point3D> m_all;
    Point3D m_x0;
    Point3D m_xt;
    double m_diff_time;
};

//inline definitions
inline Particle::Particle(double x, double y, double z, double t) : m_x0(x,y,z), m_xt(x,y,z), m_diff_time(t) {}
inline Particle::Particle(const Point3D& p, double t) : m_x0(p), m_xt(p), m_diff_time(t) {}

inline Particle::~Particle() {}
inline Particle::Particle(const Particle& p) : m_x0(p.x0()), m_xt(p.xt()), m_diff_time(p.diffTime()) {}

inline void Particle::setX0(const Point3D& p) {
    m_x0 = p;
}
inline void Particle::setX0(double x, double y, double z) {
    m_x0.set(x, y, z);
}
inline void Particle::setXt(const Point3D& p) {
    m_xt = p;
}
inline void Particle::setXt(double x, double y, double z) {
    m_xt.set(x, y, z);
}
inline void Particle::setDiffTime(double t) {
    m_diff_time = t;
}

inline Point3D Particle::x0() const {
    return m_x0;
}
inline Point3D Particle::xt() const {
    return m_xt;
}
inline double Particle::diffTime() const {
    return m_diff_time;
}

inline Particle& Particle::operator=(const Particle& p) { //return the method result
    return assign(p);
}

#endif
