/*****************************************************************************
*Class to define a particle (cpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#include "Point3D.hpp"
#include "Particle.hpp"
#include "RandNumberBetween.hpp"
#include "mptmacros.h"

void Particle::set(const Point3D& p0, const Point3D& pt, double t) {
	m_x0.set(p0.x(), p0.y(), p0.z());
	m_xt.set(pt.x(), pt.y(), pt.z());
	m_diff_time = t;
}

void Particle::set(double x0, double y0, double z0, double xt, double yt, double zt, double t) {
	m_x0.set(x0,y0,z0);
	m_xt.set(xt,yt,zt);
	m_diff_time = t;
}


void Particle::placeRandom(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high, bool init) { //different limits

	RandNumberBetween x_rand_numb(x_low, x_high);
	RandNumberBetween y_rand_numb(y_low, y_high);
	RandNumberBetween z_rand_numb(z_low, z_high);

	if (init) { //initial position of particle (generation)
		m_x0.set(x_rand_numb(), y_rand_numb(), z_rand_numb());
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
		m_diff_time = 0.;
	}
	else { //reposition
		m_xt.set(x_rand_numb(), y_rand_numb(), z_rand_numb());
	}

}

void Particle::placeRandom(double x, double y, double z, bool init) {

	if (init) { //initial position of particle (generation)
		m_x0.set(x, y, z);
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
		m_diff_time = 0.;
	}
	else  { //reposition
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
	}

}

void Particle::placeRandom(double low, double high, bool init) {

	RandNumberBetween rand_numb(low, high);

	double x = rand_numb();
	double y = rand_numb();
	double z = rand_numb();

	if (init) { //initial position of particle (generation)
		m_x0.set(x, y, z);
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
		m_diff_time = 0.;
	}
	else  { //reposition
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
	}

}

//For Particles only in Cell use this
void Particle::placeRandom(Point3D& center, double rmin, double rmax, bool init) {

	RandNumberBetween rand_numb(0., 1.);
	RandNumberBetween rand_angle(0., 1.);

	double theta_p = TWO_PI * rand_angle();
	double phi_p = acos(1. - 2. * rand_angle());
	double radial_p;

	while (true){
		radial_p = rmax * std::cbrt(rand_numb());

		if ( radial_p >= rmin ){
			break;
		}
	}

	if (init) { //initial position of particle (generation)
		m_x0.set((center.x() + sin(phi_p) * cos(theta_p) * radial_p), (center.y() + sin(phi_p) * sin(theta_p) * radial_p), (center.z() + cos(phi_p) * radial_p));
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
		m_diff_time = 0.;
	}
	else  { //reposition
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
	}

}

//random placement inside of a sphere
void Particle::placeRandomSphere(double x, double y, double z, double r_min, double r_max, bool init) {
	RandNumberBetween rand_numb(0, 1);
	double theta_p = TWO_PI * rand_numb();
	double phi_p = acos(1. - 2. * rand_numb());
	double radial_p;

	if ( r_min != r_max){
		while (true){
			radial_p = r_max * std::cbrt(rand_numb());

			if ( radial_p >= r_min ){
				break;
			}
		}
	}
	else{
		radial_p = r_max;
	}


	if (init) { //initial position of particle (generation)
		m_x0.set((x + sin(phi_p) * cos(theta_p) * radial_p), (y + sin(phi_p) * sin(theta_p) * radial_p), (z + cos(phi_p) * radial_p));
		m_xt.set(m_x0.x(), m_x0.y(), m_x0.z());
		m_diff_time = 0.;
	}
	else  { //reposition
		m_xt.set((x + sin(phi_p) * cos(theta_p) * radial_p), (y + sin(phi_p) * sin(theta_p) * radial_p), (z + cos(phi_p) * radial_p));
	}

}

void Particle::diffuse(double d, double dt, double theta_rand, double phi_rand) { 

	double theta_p = TWO_PI * theta_rand;
	double phi_p = acos(1. - 2. * phi_rand);
	double radial_p = sqrt(6 * d * dt);

	m_xt.set((m_xt.x() + sin(phi_p) * cos(theta_p) * radial_p), (m_xt.y() + sin(phi_p) * sin(theta_p) * radial_p), (m_xt.z() + cos(phi_p) * radial_p));

}

void Particle::diffuse(double x, double y, double z) {

	m_xt.set((m_xt.x() + x), (m_xt.y() + y), (m_xt.z() + z));

}

void Particle::move(Point3D& direction, double velocity, double dt) {

	double norm = sqrt(square(direction.x()) + square(direction.y()) + square(direction.z()));

	if (fabs(norm - 1.) > 1.e-10){
		double normFactor = 1./norm;
		direction.set((normFactor * direction.x()), (normFactor * direction.y()), (normFactor * direction.z()));
	}

	m_xt.set((m_xt.x() + direction.x()*velocity*dt), (m_xt.y() + direction.y()*velocity*dt), (m_xt.z() + direction.z()*velocity*dt));

}

Particle &Particle::assign(const Particle& p) {
	if (this != &p) {//if they are not the same
		m_x0 = p.x0();
		m_xt = p.xt();
		m_diff_time = p.diffTime();
	}
	return *this;
}
