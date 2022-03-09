#ifndef __MAT_H
#define __MAT_H
#include <vector>
#define ZHNMAT_VERSION                  "1.1.6"
#define NAMESPACE_ZHNMAT_L              namespace zhnmat {
#define NAMESPACE_ZHNMAT_R              }
NAMESPACE_ZHNMAT_L
#define MAT_ASSERT_ERROR(e, s)          if(!(e)){std::cout<<"Matrix Error: "<<s<<std::endl;abort();}
#define MAT_ASSERT_WARNING(e, s)        if(!(e)){std::cout<<"Matrix Warning: "<<s<<std::endl;}
#ifndef ABS
#define ABS(x)                          ((x)>=0?(x):-(x))
#endif

constexpr double EPSILON = 1e-12;
enum GENERATE_TYPE {
    NORMAL,
    DOWN_SAMPLE  // Generate downsample image pyramid.
};
enum OUTPUT_FORMAT {
    USE_BRACKET = 1 << 0,
    WRAP_AROUND = 1 << 1,
    USE_SEMICOLON = 1 << 2
};

class Mat;

struct Vector3d
{
    Vector3d(): _x(0), _y(0), _z(0) {};
    Vector3d(double x, double y, double z) :_x(x), _y(y), _z(z) {};
    Vector3d(const Vector3d& vec);

    // operator
    Vector3d operator+(const Vector3d& vec) const;
    Vector3d operator-(const Vector3d& vec) const;
    Vector3d& operator+=(const Vector3d& vec);
    Vector3d& operator-=(const Vector3d& vec);
    Vector3d operator+(const Mat& m) const;
    Vector3d operator-(const Mat& m) const;

    Vector3d operator*(double x) const;
    friend Vector3d operator*(double n, const Vector3d& m);
    double operator*(const Vector3d& vec);
    double operator*(const Mat& m);
    Vector3d& operator*=(double x);

    Vector3d operator&(const Vector3d& vec);
    Vector3d& operator=(const Vector3d& vec);
    friend std::ostream& operator<<(std::ostream &os, const Vector3d& vec);
    double norm2() const;
    Vector3d& Reset();
    Vector3d& Normalize();
    Vector3d& Reverse();
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
    Mat() :_r(0), _c(0), _p(nullptr) {}
    Mat(const Mat& m, GENERATE_TYPE type=NORMAL);
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

    // Transform to antisymmetric matrix. It requires the matrix to be size of 3-by-1.
    Mat M();

    // Inverse. method: 0 for Gaussian Elimination; 1 for LU decompose.
    // Method 1 unfinished.
    Mat inv(int method=0);

    //Unfinished.
    std::vector<double> Solve_LinearEqution(const std::vector<double>);
    std::vector<double> Solve_LeastSquare(const std::vector<double>);

    // operator
    Mat operator*(double n) const;
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
    Mat operator()(const Rect& rect) const;
    friend std::ostream& operator<<(std::ostream& os, const Mat& m);
    friend std::istream& operator>>(std::istream& is, Mat& m);  // Unfinished.
    static unsigned char OutputFormat;

private:
    int _r, _c;
    double** _p;
    void initialize();
    Mat L(), U();
};

// Return angle between two vectors.
double Vec_angle(Vector3d &v1, Vector3d &v2);

// Returns a vector perpendicular to the current vector and parallel to the ground.
Vector3d Vec_vertical(const Vector3d &v1);

// Return Identity matrix.
Mat eye(int n);

// Horizontal and vertical concatenate of two matrices.
Mat HConcat(const Mat& m1, const Mat& m2);
Mat VConcat(const Mat& m1, const Mat& m2);

// Generate a Gaussian kernel with standard deviation and side length.
Mat Gaussian_Kernel(double sigma, int sidelen);

// Situation for padding=false is not finished yet.
Mat Convolution(const Mat& m, const Mat& kernel, bool padding=true);

NAMESPACE_ZHNMAT_R
#endif
