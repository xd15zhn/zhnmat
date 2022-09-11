#include <iostream>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

unsigned char Mat::OutputFormat = USE_BRACKET | USE_SEMICOLON;
double Mat::precision = 16;
int Mat::row() const { return _r; }
int Mat::col() const { return _c; }

void Mat::initialize()
{
    if (_r<=0 || _c<=0) TRACELOG(LOG_FATAL, "Rows and columns must be greater than 0!");
    _p = new double* [_r];
    for (int i = 0; i < _r; ++i)
        _p[i] = new double[_c];
}

Mat::Mat(const Mat& m)
{
    _r = m._r; _c = m._c;
    initialize();
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] = m._p[i][j];
}
Mat::Mat(std::vector<double> data)
{
    if (data.size()<1) return;
    _r = data.size(); _c = 1;
    initialize();
    for (int i = 0; i < _r; ++i)
        _p[i][0] = data[i];
}
Mat::Mat(int r, int c, double value)
{
    _r = r; _c = c;
    initialize();
    for (int i = 0; i < _r; ++i)
        for (int j = 0; j < _c; ++j)
            _p[i][j] = value;
}
Mat::Mat(int r, int c, std::vector<double> data)
{
    if (r*c!=(int)data.size()) TRACELOG(LOG_WARNING, "Input data mismatch.");
    _r = r; _c = c;
    initialize();
    int cnt;
    for (int i = 0; i < _r; ++i){
        for (int j = 0; j < _c; ++j){
            cnt = i*_c + j;
            _p[i][j] = (cnt < (int)data.size()) ? data[cnt] : 0;
        }
    }
}
Mat::~Mat()
{
    for (int i = 0; i < _r; ++i)
        delete[] _p[i];
    delete[] _p;
}

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

Mat Mat::T(int method)
{
    Mat ans(_c, _r);
    if (method==0) {
        for (int i = 0; i < _c; ++i)
            for (int j = 0; j < _r; ++j)
                ans._p[i][j] = _p[j][i];
        return ans;
    }
    else if (method==1) {
        return ans;
    }
    return ans;
}
Mat Mat::M()
{
    if ((_r!=3) || (_c!=1)) TRACELOG(LOG_FATAL, "can't transform to antisymmetric matrix!");
    Mat ans(3, 3);
    ans._p[1][2] = -_p[0][0];
    ans._p[0][2] = +_p[1][0];
    ans._p[0][1] = -_p[2][0];
    ans._p[2][1] = +_p[0][0];
    ans._p[2][0] = -_p[1][0];
    ans._p[1][0] = +_p[2][0];
    return ans;
}
Mat Mat::inv(int method)
{
    if (_r!=_c) TRACELOG(LOG_FATAL, "can't reverse non-square matrix!");
    Mat copy(*this), ans=eye(_r);
    if (method==0) {
        double max;  // Maximum value per column
        short row;  // Row number of maximum value per column
        double *temp = new double[_r];
        double k1;
        for (short j = 0; j < _r - 1; j++) {  // j is the reference column
            // Find maximum value of column and its row numper
            max = ABS(copy._p[j][j]);
            row = j;
            for (short i = j+1; i < _r; i++) {
                if (ABS(_p[i][j]) > max) {
                    max = ABS(copy._p[i][j]);
                    row = i;
                }
            }
            // Change row of maximum value with first row
            if (row != j) {
                for (short i = j; i < _r; i++)
                    temp[i] = copy._p[row][i];
                for (short i = j; i < _r; i++)
                    copy._p[row][i] = copy._p[j][i];
                for (short i = j; i < _r; i++)
                    copy._p[j][i] = temp[i];
                for (short i = j; i < _r; i++)
                    temp[i] = ans._p[row][i];
                for (short i = j; i < _r; i++)
                    ans._p[row][i] = ans._p[j][i];
                for (short i = j; i < _r; i++)
                    ans._p[j][i] = temp[i];
            }
            //Start column elimination, that is, clear one column except the maximum row at a time
            for (short i = j + 1; i < _r; i++) {
                k1 = copy._p[i][j] / copy._p[j][j];
                copy._p[i][j] = 0;
                for (short k = 0; k < _r; k++) {
                    copy._p[i][k] -= k1 * copy._p[j][k];
                    ans._p[i][k] -= k1 * ans._p[j][k];
                }
            }
        }
        delete[] temp;
        for (short j=_r-1; j>=0; --j) {
            k1 = 1 / copy._p[j][j];
            copy._p[j][j] = 1;
            for (short k=0; k<_r; ++k)
                ans._p[j][k] *= k1;
            if (j == 0) break;
            for (short i=j-1; i>=0; --i) {
                k1 = copy._p[i][j] / copy._p[j][j];
                for (short k=0; k<_r; ++k)
                    ans._p[i][k] -= k1 * ans._p[j][k];
            }
        }
        return ans;
    }
    else if (method==1) {
        std::vector<double> y(_r);
        std::vector<double> x{1, 0, 0};
        y = Solve_LinearEqution(x);
        return ans;
    }
    return ans;
}

std::vector<double> Mat::Solve_LinearEqution(const std::vector<double> m)
{
    if ((int)m.size()!=_c) TRACELOG(LOG_FATAL, "Size mismatch!");
    if (_r!=_c) TRACELOG(LOG_FATAL, "Square matrix required, or use Solve_LeastSquare() instead.");
    if (_r<2) TRACELOG(LOG_FATAL, "Matrix dimention must be at least 2!");
    double max;  // Maximum value per column
    short row;  // Row number of maximum value per column
    double *temp = new double[_r];
    std::vector<double> ans(_r);
    for (short j = 0; j < _r - 1; j++) {  // j is the reference column
        // Find maximum value and its row numper
        max = ABS(_p[j][j]);
        row = j;
        for (short i = j+1; i < _r; i++) {
            if (ABS(_p[i][j]) > max) {
                max = ABS(_p[i][j]);
                row = i;
            }
        }
        // Change row of maximum value with first row
        if (row != j) {
            for (short i = j; i < _r; i++)
                temp[i] = _p[row][i];
            for (short i = j; i < _r; i++)
                _p[row][i] = _p[j][i];
            for (short i = j; i < _r; i++)
                _p[j][i] = temp[i];
        }
        //Start column elimination, that is, clear one column except the maximum row at a time
        for (short i = j + 1; i < _r; i++)
            for (short k = j + 1; k < _r; k++)
                _p[i][k] -= _p[i][j] / _p[j][j] * _p[j][k];
    }
    delete[] temp;
    //特征矩阵一定是降秩的,设特征向量的最后一个元素为1,从特征矩阵的倒数第二行开始回代
    ans[_r - 1] = 1;
    for (short i = _r - 2; i >= 0; i--) {
        ans[i] = 0;
        for (short j = _r - 1; j > i; j--)
            ans[i] -= _p[i][j] * ans[j];
        ans[i] /= _p[i][i];
    }
    return ans;
}

std::vector<double> Mat::Solve_LeastSquare(const std::vector<double> m)
{
    if ((int)m.size()!=_c) TRACELOG(LOG_FATAL, "Size mismatch!");
    std::vector<double> ans=m;
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
        _r = m._r;
        _c = m._c;
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
        _r = 3;
        _c = 1;
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
input and output
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
std::istream& operator>>(std::istream& is, Mat& m)
{
    for (int i=0; i<m._r; i++)
        for (int j=0; j<m._r; j++)
            is >> m._p[i][j];
    return is;
}

NAMESPACE_ZHNMAT_R
