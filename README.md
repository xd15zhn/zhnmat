# zhnmat
自用的矩阵运算库。

## 一个简单的例子
```cpp
//main.cpp
#include <iostream>
#include "zhnmat.hpp"
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
cmake_minimum_required(VERSION 3.12)
project(untitled)
set(CMAKE_BUILD_TYPE release)
add_executable(${CMAKE_PROJECT_NAME} main.cpp)
find_package(zhnmat REQUIRED)
message(STATUS "zhnmat_VERSION: ${zhnmat_VERSION}")
message(STATUS "zhnmat_DIR: ${zhnmat_DIR}")
message(STATUS "zhnmat_LIBS: ${zhnmat_LIBS}")
message(STATUS "zhnmat_INCLUDE_DIRS: ${zhnmat_INCLUDE_DIRS}")
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${zhnmat_LIBS})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${zhnmat_INCLUDE_DIRS})
```
