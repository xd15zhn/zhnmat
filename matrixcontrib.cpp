#include <iostream>
#include <cmath>
#include "zhnmat.hpp"
#include "utils.hpp"
NAMESPACE_ZHNMAT_L

Mat eye(int n)
{
    Mat ans(n, n);
    for (int i=0; i<n; ++i)
        ans.set(i, i, 1);
    return ans;
}

double AbsMat(const Mat& m)
{
    double square, ans = 0;
    for (int i=0; i<m.row(); ++i)
        for (int j=0; j<m.col(); ++j)
            ans += m.at(i, j)*m.at(i, j);
    return sqrt(ans);
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

Mat Gaussian_Kernel(double sigma, int sidelen)
{
    constexpr double PI = 3.14159265358979323846;
    if (sidelen%2==0) TRACELOG(LOG_WARNING, "Side length of Gaussian kernel should be an odd number.");
    Mat ans(sidelen, sidelen);
    sidelen >>= 1;
    double div = 0.5/(sigma*sigma);
    double divpi = div/PI;
    for (int i=-sidelen; i<=sidelen; ++i)
        for (int j=-sidelen; j<=sidelen; ++j)
            ans.set(i+sidelen, j+sidelen, divpi*exp(-(i*i+j*j)*div));
    return ans;
}

Mat Convolution(const Mat& m, const Mat& kernel, bool padding)
{
    if (m.row()<kernel.row()) TRACELOG(LOG_WARNING, "Size of kernel nust be equal or less than matrix.");
    if (m.col()<kernel.col()) TRACELOG(LOG_WARNING, "Size of kernel nust be equal or less than matrix.");
    int r = m.row(), c = m.col();
    int pad = kernel.row()>>1;
    double pixel;
    Mat ans;
    if (!padding) return ans;
    ans=Mat(m);
    int x, y;
    for (int i=0; i<r; ++i){
        for (int j=0; j<c; ++j){
            pixel = 0;
            for (int k=0; k<kernel.row(); ++k){
                for (int l=0; l<kernel.col(); ++l){
                    x = i+k-pad; y = j+l-pad;
                    if (x<0 || y<0) continue;
                    if (x>=r || y>=c) continue;
                    //std::cout<<i<<" "<<j<<" "<<k<<" "<<l<<"\n";
                    pixel += m.at(x, y) * kernel.at(k, l);
                }
            }
            ans.set(i, j, pixel);
        }
    }
    return ans;
}

NAMESPACE_ZHNMAT_R
