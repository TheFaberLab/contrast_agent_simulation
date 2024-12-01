/*****************************************************************************
*Class to define a 3d point (cpp)

*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <iostream>
#include <math.h>
#include "mptmacros.h"
#include "Point3D.hpp"

Point3D::~Point3D(){
}

void Point3D::set(double x, double y, double z){
    m_x = x;
    m_y = y;
    m_z = z;
}

void Point3D::setX(double x){
    m_x = x;
}

void Point3D::setY(double x){
    m_y = x;
}

void Point3D::setZ(double x){
    m_z = x;
}

void Point3D::print_stdout() const{
    std::cout<<"Coordinates (x,y,z) = (" <<m_x<<","<<m_y<<","<<m_z<<")"<<std::endl;
}

Point3D &Point3D::assign(const Point3D &p){
    if(this!=&p){//if they are not the same
        m_x = p.x();
        m_y = p.y();
        m_z = p.z();
    }
    return *this;
}

Point3D &Point3D::add(const Point3D &p){
    m_x += p.x();
    m_y += p.y();
    m_z += p.z();

    return *this;
}

Point3D& Point3D::subtract(const Point3D& p) {
    m_x -= p.x();
    m_y -= p.y();
    m_z -= p.z();

    return *this;
}

Point3D &Point3D::stretch(double d){
    m_x *= d;
    m_y *= d;
    m_z *= d;

    return *this;
}

bool operator==(const Point3D &p_left,const Point3D &p_right){
    return p_left.x()==p_right.x() && p_left.y()==p_right.y() && p_left.z()==p_right.z();
}

bool operator!=(const Point3D &p_left,const Point3D &p_right){
    return p_left.x()!=p_right.x() || p_left.y()!=p_right.y() || p_left.z()!=p_right.z();
}

Point3D operator+(const Point3D &p_left,const Point3D &p_right){ //left to right addition by using class methods
    return Point3D(p_left).add(p_right);
}

Point3D operator-(const Point3D& p_left, const Point3D& p_right) { //left to right addition by using class methods
    return Point3D(p_left).subtract(p_right);
}

Point3D operator*(const Point3D &p_left, double d){
    return Point3D(p_left).stretch(d);
}

Point3D operator*(double d, const Point3D &p_right){
    return Point3D(p_right).stretch(d);
}

double distance(const Point3D& p_left, const Point3D& p_right) { //distance between 2 points
    double dist = std::sqrt(square(p_left.x() - p_right.x()) + square(p_left.y() - p_right.y()) + square(p_left.z() - p_right.z()) );

    return dist;
}
