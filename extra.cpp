#include <iostream>
#include <cmath>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

double Vec_angle(Vector3d &v1, Vector3d &v2)
{
    double v1norm = v1.norm2();
    double v2norm = v2.norm2();
    if (v1norm==0 || v2norm==0)
        return 0;
    double ans = (v1*v2) / v1norm / v2norm;
    return (ans<-1) ? -1 : ((ans>1) ? 1 : ans);
}

Mat eye(int n)
{
    Mat ans(n, n);
    for (int i=0; i<n; ++i)
        ans.set(i, i, 1);
    return ans;
}

Mat Antisymmetric(const Mat& m)
{
    if ((m.row()!=3) || (m.col()!=1)) TRACELOG(LOG_FATAL, "can't transform to antisymmetric matrix!");
    Mat ans(3, 3);
    ans.set(1, 2, -m.at(0, 0));
    ans.set(0, 2, +m.at(1, 0));
    ans.set(0, 1, -m.at(2, 0));
    ans.set(2, 1, +m.at(0, 0));
    ans.set(2, 0, -m.at(1, 0));
    ans.set(1, 0, +m.at(2, 0));
    return ans;
}

Mat HConcat(const Mat& m1, const Mat& m2)
{
    if (m1.row()!=m2.row()) TRACELOG(LOG_FATAL, "Horizontal concatenate size mismatch!");
    Mat ans(m1.row(), m1.col()+m2.col());
    for (int i=0; i<m1.row(); ++i)
        for (int j=0; j<m1.col(); ++j)
            ans.set(i, j, m1.at(i, j));
    for (int i=0; i<m2.row(); ++i)
        for (int j=0; j<m2.col(); ++j)
            ans.set(i, j+m1.col(), m2.at(i, j));
    return ans;
}
Mat VConcat(const Mat& m1, const Mat& m2)
{
    if (m1.col()!=m2.col()) TRACELOG(LOG_FATAL, "Vertical concatenate size mismatch!");
    Mat ans(m1.row()+m2.row(), m1.col());
    for (int i=0; i<m1.row(); ++i)
        for (int j=0; j<m1.col(); ++j)
            ans.set(i, j, m1.at(i, j));
    for (int i=0; i<m2.row(); ++i)
        for (int j=0; j<m2.col(); ++j)
            ans.set(i+m1.row(), j, m2.at(i, j));
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

NAMESPACE_ZHNMAT_R
