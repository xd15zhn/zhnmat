#include <iostream>
#include <cmath>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

Vector3d::Vector3d(): _x(0), _y(0), _z(0) {};
Vector3d::Vector3d(double x, double y, double z) :_x(x), _y(y), _z(z) {};
Vector3d::Vector3d(const Vector3d &vec) { _x=vec._x, _y=vec._y, _z=vec._z; }
Vector3d& Vector3d::operator=(const Vector3d &vec) { _x=vec._x, _y=vec._y, _z=vec._z; return *this; }
Vector3d& Vector3d::operator=(const Mat& m) {
    if (m.row()!=3 || m.col()!=1) TRACELOG(LOG_FATAL, "Size mismatch in equation constructor!");
    _x=m.at(0,0), _y=m.at(1,0), _z=m.at(2,0); return *this;
}
Vector3d::Vector3d(const Mat& m) {
    if (m.row()!=3 || m.col()!=1) TRACELOG(LOG_FATAL, "Size mismatch in equation constructor!");
    _x=m.at(0,0), _y=m.at(1,0), _z=m.at(2,0);
}
/**********************
operator addition and subtraction
**********************/
Vector3d Vector3d::operator+(const Vector3d &vec) const { return Vector3d(_x + vec._x, _y + vec._y, _z + vec._z); }
Vector3d Vector3d::operator-(const Vector3d &vec) const { return Vector3d(_x - vec._x, _y - vec._y, _z - vec._z); }
Vector3d& Vector3d::operator+=(const Vector3d &vec) { _x+=vec._x, _y+=vec._y, _z+=vec._z; return *this; }
Vector3d& Vector3d::operator-=(const Vector3d &vec) { _x-=vec._x, _y-=vec._y, _z-=vec._z; return *this; }
Vector3d Vector3d::operator+(const Mat& m) const
{
    if (m.row()!=3 || m.col()!=1) TRACELOG(LOG_FATAL, "Size mismatch! Addition between vector and matrix.");
    return Vector3d(_x+m.at(0,0), _y+m.at(1,0), _z+m.at(2,0));
};
Vector3d Vector3d::operator-(const Mat& m) const
{
    if (m.row()!=3 || m.col()!=1) TRACELOG(LOG_FATAL, "Size mismatch! Subtraction between vector and matrix.");
    return Vector3d(_x-m.at(0,0), _y-m.at(1,0), _z-m.at(2,0));
};

/**********************
operator multiplication
**********************/
Vector3d Vector3d::operator*(double x) const { return Vector3d(_x * x, _y * x, _z * x); }
Vector3d Vector3d::operator/(double x) const { return Vector3d(_x / x, _y / x, _z / x); }
Vector3d operator*(double n, const Vector3d& m) { return m*n; }
double Vector3d::operator*(const Vector3d& vec) { return _x*vec._x + _y*vec._y + _z*vec._z; }
Vector3d& Vector3d::operator*=(double x) { _x*=x, _y*=x, _z*=x; return *this; }
Vector3d& Vector3d::operator/=(double x) { _x/=x, _y/=x, _z/=x; return *this; }
double Vector3d::operator*(const Mat& m)
{
    if (m.row()!=3 || m.col()!=1) TRACELOG(LOG_FATAL, "Size mismatch! multiplication between vector and matrix.");
    return _x*m.at(0,0) + _y*m.at(1,0) + _z*m.at(2,0);
}

/**********************
other member functions
**********************/
Vector3d Vector3d::operator&(const Vector3d &vec)
{
    double vx = _y*vec._z - _z*vec._y;
    double vy = _z*vec._x - _x*vec._z;
    double vz = _x*vec._y - _y*vec._x;
    return Vector3d(vx, vy, vz);
};
std::ostream& operator<<(std::ostream &os, const Vector3d& vec)
{
    os<<"{"<<vec._x<<", "<<vec._y<<", "<<vec._z<<"}";
    return os;
}
Vector3d &Vector3d::Reset() { _x=0;_y=0;_z=0;return *this; }
Vector3d &Vector3d::Reverse() { _x=-_x;_y=-_y;_z=-_z;return *this; }
double Vector3d::norm2() const{ return sqrt(_x * _x + _y * _y + _z * _z); }
Vector3d Vector3d::Normalvector(void) const {
	double len = norm2();
    if (len==0) return Vector3d(0, 0, 0);
    return Vector3d(_x/len, _y/len, _z/len);
}
Vector3d& Vector3d::Normalize(void)
{
	double len = norm2();
    if (len==0) return *this;
    _x /= len; _y /= len; _z /= len;
    return *this;
}

NAMESPACE_ZHNMAT_R
