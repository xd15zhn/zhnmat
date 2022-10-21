#ifndef __MAT_H
#define __MAT_H
#include <vector>
#define NAMESPACE_ZHNMAT_L              namespace zhnmat {
#define NAMESPACE_ZHNMAT_R              }
#ifndef ABS
#define ABS(x)                          ((x)>=0?(x):-(x))
#endif
NAMESPACE_ZHNMAT_L

constexpr double EPSILON = 1e-12;
enum OUTPUT_FORMAT {
    USE_BRACKET = 1 << 0,
    WRAP_AROUND = 1 << 1,
    USE_SEMICOLON = 1 << 2
};

class Mat;

struct Vector3d
{
    // constructor function
    Vector3d();
    Vector3d(double x, double y, double z);
    Vector3d(const Vector3d& vec);
    Vector3d& operator=(const Vector3d& vec);

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

struct Rect
{
    Rect(int x, int y, int w, int h):
        _x(x), _y(y), _w(w), _h(h) {};
    int _x, _y, _w, _h;
};

class Mat
{
public:
    Mat();
    Mat(const Mat& m);
    Mat(std::vector<double> data);
    Mat(int r, int c, double value=0);
    Mat(int r, int c, std::vector<double> data);
    ~Mat();

    // Return how many rows and columns.
    int row() const;
    int col() const;

    // Get and set specified pixels.
    double at(int r, int c) const;
    void set(int r, int c, double value=0);

    // Transpose. method: 0 for traditional method; 1 for fast method.
    // Method 1 unfinished.
    Mat T(int method=0);

    // Inverse. method: 0 for Gaussian Elimination; 1 for LU decompose.
    // Method 1 unfinished.
    Mat inv(int method=0);

    //Unfinished.
    std::vector<double> Solve_LinearEqution(const std::vector<double>);
    std::vector<double> Solve_LeastSquare(const std::vector<double>);

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
    friend std::istream& operator>>(std::istream& is, Mat& m);  // Unfinished.
    static unsigned char OutputFormat;
    static double precision;

private:
    int _r, _c;
    double** _p;
    void initialize();
    Mat L(), U();
};

/**********************
functions below can be removed when building the project to save space.
**********************/
#pragma region extra
// Return angle between two vectors.
double Vec_angle(Vector3d &v1, Vector3d &v2);

// Return Identity matrix.
Mat eye(int n);

// Transform to antisymmetric matrix. It requires the matrix to be size of 3-by-1.
Mat Antisymmetric(const Mat& m);

// Return the Euclid length of a matrix.
double AbsMat(const Mat& m);

// Horizontal and vertical concatenate of two matrices.
Mat HConcat(const Mat& m1, const Mat& m2);
Mat VConcat(const Mat& m1, const Mat& m2);

// Generate a Gaussian kernel with standard deviation and side length.
Mat Gaussian_Kernel(double sigma, int sidelen);

// Situation for padding=false is not finished yet.
Mat Convolution(const Mat& m, const Mat& kernel, bool padding=true);
#pragma endregion extra

NAMESPACE_ZHNMAT_R
#endif
