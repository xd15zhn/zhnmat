#include <iostream>
#include <cmath>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

#ifdef USE_EXTRA
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
#endif


/**********************
QR分解
**********************/
// 列向量的2-范数
double Vector_Norm2(const double *v, const short len) {
	double ans = 0;
	for (short i = 0; i < len; i++)
		ans += v[i] * v[i];
	return sqrt(ans);
}
//求使向量v变换到基向量的法向量u,返回u的1/2模平方,len为两个向量的维度
double Vector_Normal(double *u, const double *v, const short len) {
	for (short i = 0; i < len; i++)
		u[i] = v[i];
	double sigma = Vector_Norm2(u, len);
	if (u[0] < 0)
		sigma = -sigma;
	double rho = sigma*(sigma + u[0]);
	u[0] += sigma;
	return rho;
}
//豪斯霍尔德矩阵法求矩阵m的拟上三角矩阵
Mat Matrix_Hessenberg(Mat m) {
    Mat ans(m);
	double rho, dotsum;
	double u[m.row()-1];  // Householder变换阵的法向量
	double v[m.row()-1];  // 组成下三角阵的列向量
	for (short h = 0; h < m.row()-2; h++) {  // H矩阵编号
		for (short i = h; i < m.row()-1; i++)
			v[i] = m.at(i+1, h);  // 提取列向量
		rho = Vector_Normal(u+h, v + h, m.row() - h-1);  // 确定出对应阶数的H矩阵
		for (short j = 0; j < m.row(); j++) {  // 每一列分别左乘H矩阵
			dotsum = 0;
			for (short i = h; i < m.row()-1; i++)
				dotsum += u[i] * m.at(i+1, j);  // 点乘
			dotsum /= rho;
			for (short i = h; i < m.row()-1; i++)
                m.set(i+1, j, m.at(i+1, j)-dotsum*u[i]);
		}
		for (short i = 0; i < m.row(); i++) {  //每一行分别右乘H矩阵
			dotsum = 0;
			for (short j = h; j < m.row()-1; j++)
				dotsum += u[j] * m.at(i, j+1);
			dotsum /= rho;
			for (short j = h; j < m.row()-1; j++)
                m.set(i, j+1, m.at(i, j+1)-dotsum*u[j]);
		}
	}
    return ans;
}
//将拟上三角矩阵M分别左右乘其豪斯霍尔德矩阵H,并覆盖
Mat Matrix_Householder(Mat m)
{
    Mat ans(m);
	double rho, dotsum;
	double u[2];
	double v[2];
	double U[3][m.row()-1];
	for (short h = 0; h < m.row()-1; h++)
	{
		v[0] = m.at(h, h);
		v[1] = m.at(h+1, h);
		rho = Vector_Normal(u, v, 2);
		U[0][h] = u[0];
		U[1][h] = u[1];
		U[2][h] = rho;
		for (short j = 0; j < m.row(); j++)
		{
			dotsum = u[0] * m.at(h, j);
			dotsum += u[1] * m.at(h+1, j);
			dotsum /= rho;
            ans.set(h, j, m.at(h, j)-dotsum*u[0]);
            ans.set(h+1, j, m.at(h+1, j)-dotsum*u[1]);
		}
	}
	for (short h = 0; h < m.row()-1; h++)
	{
		u[0] = U[0][h];
		u[1] = U[1][h];
		rho = U[2][h];
		for (short i = 0; i < m.row(); i++)
		{
			dotsum = u[0] * m.at(i, h);
			dotsum += u[1] * m.at(i, h+1);
			dotsum /= rho;
            ans.set(i, h, m.at(i, h)-dotsum*u[0]);
            ans.set(i, h+1, m.at(i, h+1)-dotsum*u[1]);
		}
	}
    return ans;
}

NAMESPACE_ZHNMAT_R
