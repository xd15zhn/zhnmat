#include <iostream>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L
#define DELTA_ZERO                      1e-12
#ifndef ABS
#define ABS(x)                          ((x)>=0?(x):-(x))
#endif
#ifndef LIMIT
#define LIMIT(x, min, max)           (((x)<=(min) ? (min) : ((x)>=(max) ? (max) : (x))))
#endif

Mat::Mat() :_r(0), _c(0), _p(nullptr) {}
int Mat::row() const { return _r; }
int Mat::col() const { return _c; }
int Mat::_precision = 15;
OutputFormat Mat::_format;

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
Mat::Mat(const Vector3d& vec) {
    _r = 3; _c = 1;
    initialize();
    _p[0][0] = vec._x;
    _p[1][0] = vec._y;
    _p[2][0] = vec._z;
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
Mat Mat::atr(int r) const {
    Mat ans(1, _c);
    for (int i = 0; i < _c; i++)
        ans._p[0][i] = _p[r][i];
    return ans;    
}
Mat Mat::atc(int c) const {
    Mat ans(_r, 1);
    for (int i = 0; i < _r; i++)
        ans._p[i][0] = _p[i][c];
    return ans;    
}
void Mat::setr(int r, Mat m) const {
    if (m._r !=1 || m._c != _c)
        TRACELOG(LOG_FATAL, "setr indexes out of range!"
        " Input rows and cols:%d,%d;  Max rows and cols:%d,%d.", m._r,m._c,_r,_c);
    for (int i = 0; i < _c; i++)
        _p[r][i] = m._p[0][i];
}
void Mat::setc(int c, Mat m) const {
    if (m._c !=1 || m._r != _r)
        TRACELOG(LOG_FATAL, "setr indexes out of range!"
        " Input rows and cols:%d,%d;  Max rows and cols:%d,%d.", m._r,m._c,_r,_c);
    for (int i = 0; i < _c; i++)
        _p[i][c] = m._p[i][0];
}
std::vector<double> Mat::To_Vector() const {
    std::vector<double> ans;
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            ans.push_back(_p[i][j]);
    return ans;
}

Mat Mat::T()
{
    Mat ans(_c, _r);
    for (int i = 0; i < _c; ++i)
        for (int j = 0; j < _r; ++j)
            ans._p[i][j] = _p[j][i];
    return ans;
}

/*初等变换法求逆矩阵*/
Mat Mat::inv()
{
    if (_r!=_c) TRACELOG(LOG_FATAL, "can't reverse non-square matrix!");
    Mat copy(*this), ans(_r, _r);
    for (int i=0; i<_r; ++i) ans.set(i, i, 1);
    double max;  // 每一列的最大值
    double absvalue;  // 绝对值暂存
    short row;  // 每一列最大值的行号
    double temp;  // 交换两行时暂存
    double k1;
    /* 化为上三角矩阵 */
    for (short j = 0; j < _r - 1; j++) {  // 第j列
        max = ABS(copy._p[j][j]); row = j;
        for (short i = j+1; i < _r; i++) {
            absvalue = ABS(copy._p[i][j]);
            if (absvalue > max) {  // 寻找每一列的最大值和对应的行号
                max = absvalue; row = i;
            }
        }
        // 将第j列最大值所在的行与第j行交换
        if (row != j) {
            for (short i = j; i < _r; i++) {
                temp = copy._p[row][i];
                copy._p[row][i] = copy._p[j][i];
                copy._p[j][i] = temp;
            }
            for (short i = 0; i < _r; i++) {
                temp = ans._p[row][i];
                ans._p[row][i] = ans._p[j][i];
                ans._p[j][i] = temp;
            }
        }
        // 开始按列消元,即每次将除最大值行以外的一列清零
        for (short i = j + 1; i < _r; i++) {
            k1 = copy._p[i][j] / copy._p[j][j];
            for (short k = j; k < _r; k++)
                copy._p[i][k] -= k1 * copy._p[j][k];
            for (short k = 0; k < _r; k++)
                ans._p[i][k] -= k1 * ans._p[j][k];
        }
    }
    if (ABS(copy._p[_r-1][_c-1]) < DELTA_ZERO) TRACELOG(LOG_FATAL, "Singular matrix!");
    /* 将上三角矩阵化为单位阵 */
    for (short j=_r-1; j>=0; --j) {
        k1 = 1 / copy._p[j][j];
        for (short k=0; k<_r; ++k)
            ans._p[j][k] *= k1;
        if (j == 0) break;
        for (short i=j-1; i>=0; --i) {
            k1 = copy._p[i][j];
            for (short k=0; k<_r; ++k)
                ans._p[i][k] -= k1 * ans._p[j][k];
        }
    }
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
    if (_r<rect._h || _c<rect._w) TRACELOG(LOG_FATAL, "Rect out of range!");
    Mat ans(rect._h, rect._w);
    for (int i = 0; i < ans._r; ++i)
        for (int j = 0; j < ans._c; ++j)
            ans._p[i][j] = this->_p[i+rect._y][j+rect._x];
    return ans;
}

/**********************
Print
**********************/
void Mat::Set_Precision(int precision) {
    _precision = LIMIT(precision, 1, 32);
}
std::ostream& operator<<(std::ostream& os, const Mat& m) {
    os.precision(Mat::_precision);
    os << Mat::_format.prefix;
    for (int i=0; i<m._r-1; ++i){
        for (int j=0; j<m._c-1; ++j)
            os << m._p[i][j] << Mat::_format.rowin;
        os << m._p[i][m._c-1] << Mat::_format.rowout;
    }
    for (int j=0; j<m._c-1; ++j)
        os << m._p[m._r-1][j] << Mat::_format.rowin;
    os << m._p[m._r-1][m._c-1] << Mat::_format.suffix;
    return os;
}

NAMESPACE_ZHNMAT_R
