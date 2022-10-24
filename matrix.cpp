#include <iostream>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

unsigned char Mat::OutputFormat = USE_BRACKET | USE_SEMICOLON;
double Mat::precision = 16;
Mat::Mat() :_r(0), _c(0), _p(nullptr) {}
int Mat::row() const { return _r; }
int Mat::col() const { return _c; }

/**********************
Construction and Destruction
**********************/
void Mat::initialize() {
    if (_r<=0 || _c<=0) TRACELOG(LOG_FATAL, "Rows and columns must be greater than 0!");
    _p = new double* [_r];
    for (int i = 0; i < _r; ++i)
        _p[i] = new double[_c];
}
Mat::Mat(const Mat& m) {
    _r = m._r; _c = m._c;
    initialize();
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] = m._p[i][j];
}
Mat::Mat(std::vector<double> data) {
    if (data.size()<1) return;
    _r = data.size(); _c = 1;
    initialize();
    for (int i = 0; i < _r; ++i)
        _p[i][0] = data[i];
}
Mat::Mat(int r, int c, double value) {
    _r = r; _c = c;
    initialize();
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] = value;
}
Mat::Mat(int r, int c, std::vector<double> data) {
    if (r*c!=(int)data.size()) TRACELOG(LOG_WARNING, "Input data mismatch.");
    _r = r; _c = c;
    initialize();
    int cnt;
    for (int i = 0; i < _r; ++i) {
        for (int j = 0; j < _c; ++j) {
            cnt = i*_c + j;
            _p[i][j] = (cnt < (int)data.size()) ? data[cnt] : 0;
        }
    }
}
Mat::~Mat() {
    for (int i = 0; i < _r; ++i)
        delete[] _p[i];
    delete[] _p;
}

/**********************
Frequently used operation
**********************/
double Mat::at(int r, int c) const
{
    if (r>=_r || c>=_c)
        TRACELOG(LOG_FATAL, "Read indexes out of range!"
        " Input rows and cols:%d,%d;  Max rows and cols:%d,%d.", r,c,_r,_c);
    return _p[r][c];
}
void Mat::set(int r, int c, double value)
{
    if (r>=_r || c>=_c)
        TRACELOG(LOG_FATAL, "Write indexes out of range!"
        " Input rows and cols:%d,%d;  Max rows and cols:%d,%d.", r,c,_r,_c);
    _p[r][c]=value;
}
Mat Mat::T()
{
    Mat ans(_c, _r);
    for (int i = 0; i < _c; ++i)
        for (int j = 0; j < _r; ++j)
            ans._p[i][j] = _p[j][i];
    return ans;
}

Mat Mat::inv()
{
    Mat ans;
    return ans;
}

/**********************
operator multiplication
**********************/
Mat Mat::operator*(double n) const
{
    Mat ans(_r, _c);
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            ans._p[i][j] = _p[i][j] * n;
    return ans;
}
Mat Mat::operator/(double n) const
{
    Mat ans(_r, _c);
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            ans._p[i][j] = _p[i][j] / n;
    return ans;
}
Mat Mat::operator*(const Mat& m) const
{
    if (_c!=m._r) TRACELOG(LOG_FATAL, "Size mismatch! Multiplication between matrix.");
    Mat ans(_r, m._c);
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < m._c; ++j)
            for (int k = 0; k < _c; k++)
                ans._p[i][j] += (_p[i][k] * m._p[k][j]);
    return ans;
}
Vector3d Mat::operator*(const Vector3d &vec) const
{
    if (_r!=3 || _c!=3) TRACELOG(LOG_FATAL, "Size mismatch! Multiplication between matrix and vector.");
    Vector3d ans;
    ans._x = _p[0][0]*vec._x + _p[0][1]*vec._y + _p[0][2]*vec._z;
    ans._y = _p[1][0]*vec._x + _p[1][1]*vec._y + _p[1][2]*vec._z;
    ans._z = _p[2][0]*vec._x + _p[2][1]*vec._y + _p[2][2]*vec._z;
    return ans;
};
Mat operator*(double n, const Mat& m) { return m*n; }
Mat& Mat::operator*=(double n)
{
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
                _p[i][j] *= n;
    return *this;
}
Mat& Mat::operator*=(const Mat& m)
{
    if (_c!=m._r) TRACELOG(LOG_FATAL, "Size mismatch! Multiplication between matrix.");
    Mat temp(_r, m._c);
    for (int i = 0; i < temp._r; ++i)
        for (int j = 0; j < temp._c; ++j)
            for (int k = 0; k < _c; ++k)
                temp._p[i][j] += (_p[i][k] * m._p[k][j]);
    *this = temp;
    return *this;
}

/**********************
operator addition and subtraction
**********************/
Mat Mat::operator+(const Mat& m) const
{
    if (_r!=m._r || _c!=m._c) TRACELOG(LOG_FATAL, "Size mismatch! Addition between matrix.");
    Mat ans(_r, _c, 0.0);
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < m._c; ++j)
            ans._p[i][j] = _p[i][j] + m._p[i][j];
    return ans;
}
Mat Mat::operator-(const Mat& m) const
{
    if (_r!=m._r || _c!=m._c) TRACELOG(LOG_FATAL, "Size mismatch! Subtraction between matrix.");
    Mat ans(_r, _c, 0.0);
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < m._c; ++j)
            ans._p[i][j] = _p[i][j] - m._p[i][j];
    return ans;
}
Mat& Mat::operator+=(const Mat& m)
{
    if (_r!=m._r || _c!=m._c) TRACELOG(LOG_FATAL, "Size mismatch! Addition between matrix.");
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] += m._p[i][j];
    return *this;
}
Mat& Mat::operator-=(const Mat& m)
{
    if (_r!=m._r || _c!=m._c) TRACELOG(LOG_FATAL, "Size mismatch! Subtraction between matrix.");
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] -= m._p[i][j];
    return *this;
}
Vector3d Mat::operator+(const Vector3d& vec) const
{
    if (_r!=3 || _c!=1) TRACELOG(LOG_FATAL, "Size mismatch! Addition between matrix and vector.");
    return Vector3d(_p[0][0] + vec._x, _p[1][0] + vec._y, _p[2][0] + vec._z);
};
Vector3d Mat::operator-(const Vector3d& vec) const
{
    if (_r!=3 || _c!=1) TRACELOG(LOG_FATAL, "Size mismatch! Subtraction between matrix and vector.");
    return Vector3d(_p[0][0] - vec._x, _p[1][0] - vec._y, _p[2][0] - vec._z);
};

/**********************
assignment
**********************/
Mat& Mat::operator=(const Mat& m)
{
    if (this == &m) return *this;
    if (_r != m._r || _c != m._c) {
        if (_p != nullptr){
            for (int i = 0; i < _r; ++i)
                delete[] _p[i];
            delete[] _p;
        }
        _r = m._r; _c = m._c;
        initialize();
    }
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] = m._p[i][j];
    return *this;
}
Mat& Mat::operator=(const Vector3d& vec)
{
    if (_r != 3 || _c != 1) {
        if (_p != nullptr){
            for (int i = 0; i < _r; ++i)
                delete[] _p[i];
            delete[] _p;
        }
        _r = 3; _c = 1;
        initialize();
    }
    _p[0][0] = vec._x;
    _p[1][0] = vec._y;
    _p[2][0] = vec._z;
    return *this;
}
Mat Mat::operator()(const Rect &rect) const
{
    if (_r<=0 || _c<=0) TRACELOG(LOG_FATAL, "Matrix haven't been initialized!");
    if (rect._w<=0 && rect._h<=0) TRACELOG(LOG_FATAL, "Rect error!");
    if (_r<rect._w || _c<rect._h) TRACELOG(LOG_FATAL, "Rect out of range!");
    Mat ans(rect._h, rect._w);
    for (int i = 0; i < ans._r; ++i)
        for (int j = 0; j < ans._c; ++j)
            ans._p[i][j] = this->_p[i+rect._y][j+rect._x];
    return ans;
}

/**********************
Print
**********************/
std::ostream& operator<<(std::ostream& os, const Mat& m)
{
    if ((Mat::OutputFormat&0xF8)!=0) TRACELOG(LOG_WARNING, "Wrong print format were given.");
    unsigned char opf = Mat::OutputFormat & 0x07;
    if (Mat::precision<=1) TRACELOG(LOG_WARNING, "Wrong print precision were given.");
    
    os.precision(Mat::precision);
    if (opf & USE_BRACKET) os << "[";
    if (opf & WRAP_AROUND) os << std::endl;
    for (int i=0; i<m._r; ++i){
        for (int j=0; j<m._c-1; ++j)
            os << m._p[i][j] << ", ";
        os << m._p[i][m._c-1] << ((opf & USE_SEMICOLON)? "; " : ", ");
        if (opf & WRAP_AROUND) os << std::endl;
    }
    if (opf & USE_BRACKET) os << "]";
    return os;
}

NAMESPACE_ZHNMAT_R
