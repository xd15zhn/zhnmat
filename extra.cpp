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

Mat Insertr(const Mat& m, uint r, const Mat& data) {
    if (r > m.row()) TRACELOG(LOG_FATAL, "Insert row error! %d and %d.", r, m.row());
    Mat ans1 = m(Rect(0, 0, r, m.col()));
    Mat ans2 = m(Rect(r, 0, m.row()-r, m.col()));
    return VConcat(VConcat(ans1, data), ans2);
}
Mat Deleter(const Mat& m, uint r) {
    if (r >= m.row()) TRACELOG(LOG_FATAL, "Delete row error! %d and %d.", r, m.row());
    Mat ans1 = m(Rect(0, 0, r, m.col()));
    Mat ans2 = m(Rect(r+1, 0, m.row()-r-1, m.col()));
    return VConcat(ans1, ans2);
}
Mat Insertc(const Mat& m, uint c, const Mat& data) {
    if (c > m.col()) TRACELOG(LOG_FATAL, "Insert row error! %d and %d.", c, m.col());
    Mat ans1 = m(Rect(0, 0, m.row(), c));
    Mat ans2 = m(Rect(0, c, m.row(), m.col()-c));
    return HConcat(HConcat(ans1, data), ans2);
}
Mat Deletec(const Mat& m, uint c) {
    if (c >= m.col()) TRACELOG(LOG_FATAL, "Delete row error! %d and %d.", c, m.col());
    Mat ans1 = m(Rect(0, 0, m.row(), c));
    Mat ans2 = m(Rect(0, c+1, m.row(), m.col()-c-1));
    return HConcat(ans1, ans2);
}

Mat HConcat(const Mat& m1, const Mat& m2)
{
    if (m1.col()==0) return m2;
    if (m2.col()==0) return m1;
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
    if (m1.row()==0) return m2;
    if (m2.row()==0) return m1;
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
NAMESPACE_ZHNMAT_R
