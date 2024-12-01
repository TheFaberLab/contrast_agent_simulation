/*****************************************************************************
*Class to define a proton (header .hpp)
* derived from particle

*16/05/2022 Lauritz Klünder
******************************************************************************/

#ifndef PROTON_HPP_
#define PROTON_HPP_

#include "Point3D.hpp"
#include "Particle.hpp"

class Proton : public Particle {
public:

	inline Proton(double x, double y, double z, double t, double m);
	inline Proton(const Point3D& p, double t, double m);
	inline Proton(double t, double m);

	inline void setPhase(double p) {
		m_phase = p;
	}

	inline double setPhase() {
		return m_phase;
	}


	inline void inCell(bool in) {
		m_inCell = in;
	}

	inline bool inCell() {
		return m_inCell;
	}

private:
	double m_phase;
	bool m_inCell;

};

inline Proton::Proton(double x, double y, double z, double t, double m) : Particle(x,y,z,t), m_phase(m), m_inCell() {}
inline Proton::Proton(const Point3D& p, double t, double m) : Particle(p,t), m_phase(m), m_inCell() {}
inline Proton::Proton(double t, double m) : Particle(0.,0.,0.,t), m_phase(m), m_inCell() {}

#endif
