#ifndef __MAT_H
#define __MAT_H
#include <vector>
#define NAMESPACE_ZHNMAT_L              namespace zhnmat {
#define NAMESPACE_ZHNMAT_R              }
NAMESPACE_ZHNMAT_L

constexpr double EPSILON = 1e-12;

class Mat;

/**********************
Output format. Default:
[1, 2, 3; 4, 5, 6; 7, 8, 9]
**********************/
struct OutputFormat {
    std::string prefix = "[";
    std::string rowin = ", ";
    std::string rowout = "; ";
    std::string suffix = "]";
};

struct Rect
{
    Rect(int y, int x, int h, int w):
        _x(x), _y(y), _w(w), _h(h) {};
    int _x, _y, _w, _h;
};

struct Vector3d
{
    // constructor function
    Vector3d();
    Vector3d(double x, double y, double z);
    Vector3d(const Vector3d& vec);
    Vector3d(const Mat& m);
    Vector3d& operator=(const Vector3d& vec);
    Vector3d& operator=(const Mat& m);

    // operator addition and subtraction
    Vector3d operator+(const Vector3d& vec) const;
    Vector3d operator-(const Vector3d& vec) const;
    Vector3d& operator+=(const Vector3d& vec);
    Vector3d& operator-=(const Vector3d& vec);
    Vector3d operator+(const Mat& m) const;
    Vector3d operator-(const Mat& m) const;

    // operator multiplication
    Vector3d operator*(double x) const;  // vector scalar multiplication
    Vector3d operator/(double x) const;  // vector scalar division
    friend Vector3d operator*(double n, const Vector3d& m);  // vector scalar multiplication
    double operator*(const Vector3d& vec);  // vector dot product
    double operator*(const Mat& m);  // vector dot product
    Vector3d& operator*=(double x);  // vector scalar multiplication
    Vector3d& operator/=(double x);  // vector scalar division

    Vector3d operator&(const Vector3d& vec);  // vector cross product
    friend std::ostream& operator<<(std::ostream &os, const Vector3d& vec);
    Vector3d& Reset();
    Vector3d& Reverse();
    double norm2() const;
    Vector3d Normalvector() const;  // Return the normal vector of this vector.
    Vector3d& Normalize();  // Normalize this vector and return itself.
    double _x, _y, _z;
};

class Mat
{
public:
    Mat();
    Mat(const Mat& m);
    Mat(const Vector3d& vec);
    Mat(std::vector<double> data);
    Mat(int r, int c, double value=0);
    Mat(int r, int c, std::vector<double> data);
    ~Mat();

    // Return how many rows and columns.
    int row() const;
    int col() const;

    // Get and set specified element.
    double at(int r, int c) const;
    void set(int r, int c, double value=0);
    // Get and set elements in a specified row or col.
    Mat atr(int r) const;
    Mat atc(int r) const;
    void setr(int r, Mat m) const;
    void setc(int c, Mat m) const;

    // Matrix transpose.
    Mat T();
    // Matrix inverse.
    Mat inv();

    // operator
    Mat operator*(double n) const;
    Mat operator/(double n) const;
    Mat operator*(const Mat& m) const;
    Vector3d operator*(const Vector3d& vec) const;
    friend Mat operator*(double n, const Mat& m);
    Mat& operator*=(double n);
    Mat& operator*=(const Mat& m);

    Mat operator+(const Mat& m) const;
    Mat operator-(const Mat& m) const;
    Mat& operator+=(const Mat& m);
    Mat& operator-=(const Mat& m);
    Vector3d operator+(const Vector3d& vec) const;
    Vector3d operator-(const Vector3d& vec) const;

    Mat& operator=(const Mat& m);
    Mat& operator=(const Vector3d& vec);
    Mat operator()(const Rect& rect) const;
    friend std::ostream& operator<<(std::ostream& os, const Mat& m);
    static void Set_Precision(int precision=15);
    static OutputFormat _format;

private:
    int _r, _c;
    double** _p;
    void initialize();
    static int _precision;
};

/**********************
functions below can be removed to save space when building the project.
**********************/
#pragma region extra
// Return angle between two vectors.
double Vec_angle(Vector3d &v1, Vector3d &v2);

// Return Identity matrix.
Mat eye(int n);

// Transform to antisymmetric matrix. It requires the matrix to be size of 3-by-1.
Mat Antisymmetric(const Mat& m);

// Horizontal and vertical concatenate of two matrices.
Mat HConcat(const Mat& m1, const Mat& m2);
Mat VConcat(const Mat& m1, const Mat& m2);

#pragma endregion extra

NAMESPACE_ZHNMAT_R
#endif
