#include "zhnmat.hpp"
NAMESPACE_ZHNMAT_L

Vector3d::Vector3d(const Vector3d &vec) { _x=vec._x, _y=vec._y, _z=vec._z; }
Vector3d Vector3d::operator+(const Vector3d &vec) const { return Vector3d(_x + vec._x, _y + vec._y, _z + vec._z); }
Vector3d Vector3d::operator-(const Vector3d &vec) const { return Vector3d(_x - vec._x, _y - vec._y, _z - vec._z); }
Vector3d Vector3d::operator*(double x) const { return Vector3d(_x * x, _y * x, _z * x); }
Vector3d operator*(double n, const Vector3d& m) { return m*n; }
Vector3d& Vector3d::operator=(const Vector3d &vec) { _x=vec._x, _y=vec._y, _z=vec._z; return *this; }
Vector3d& Vector3d::operator+=(const Vector3d &vec) { _x+=vec._x, _y+=vec._y, _z+=vec._z; return *this; }
Vector3d& Vector3d::operator-=(const Vector3d &vec) { _x-=vec._x, _y-=vec._y, _z-=vec._z; return *this; }
Vector3d& Vector3d::operator*=(double x) { _x*=x, _y*=x, _z*=x; return *this; }
double Vector3d::operator*(const Vector3d& vec) { return _x*vec._x + _y*vec._y + _z*vec._z; }
Vector3d &Vector3d::Reset() { _x=0;_y=0;_z=0;return *this; }
Vector3d &Vector3d::Reverse() { _x=-_x;_y=-_y;_z=-_z;return *this; }
double Vector3d::norm2() const{ return sqrt(_x * _x + _y * _y + _z * _z); }

Vector3d& Vector3d::Normalize(void)
{
	double len = norm2();
	_x /= len;
	_y /= len;
	_z /= len;
    return *this;
}

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

Vector3d Vector3d::operator+(const Mat& m) const
{
    MAT_ASSERT_ERROR(m.col()==3, "Size mismatch!");
    MAT_ASSERT_ERROR(m.row()==1, "Size mismatch!");
    return Vector3d(_x+m.at(0,0), _y+m.at(1,0), _z+m.at(2,0));
};

double Vec_angle(Vector3d &v1, Vector3d &v2)
{
    return (v1*v2) / v1.norm2() / v2.norm2();
}

Vector3d Vec_vertical(const Vector3d &v1)
{
    if (ABS(v1._y) <= 1e-6)
        return Vector3d(0, 1, 0);
	else
        return Vector3d(-1, v1._x/v1._y, 0);
}

NAMESPACE_ZHNMAT_R
