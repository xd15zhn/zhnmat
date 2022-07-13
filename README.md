# zhnmat
自用的矩阵运算库。

## 一个简单的例子
```cpp
//main.cpp
#include <iostream>
#include "zhnmat/zhnmat.hpp"
using namespace std;
using namespace zhnmat;
int main()
{
    Mat A(2, 2, vector<double>{1, 2, 3, 4});
    Mat B(2, 2, vector<double>{1, 2, 3, 4});
    cout << A*B << endl;
}
```
```
# CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(untitled)
set(CMAKE_BUILD_TYPE release)
add_executable(${CMAKE_PROJECT_NAME}_exe main.cpp)
find_package(zhnmat REQUIRED)
message(STATUS "zhnmat config: ${zhnmat_DIR}")
target_link_libraries(${CMAKE_PROJECT_NAME}_exe PRIVATE zhnmat)
```
