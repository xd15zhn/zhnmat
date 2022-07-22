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
add_executable(${CMAKE_PROJECT_NAME}_exe main.cpp)
# list(APPEND CMAKE_PREFIX_PATH "E:/cpplibraries/")
find_package(zhnmat REQUIRED)
# message(STATUS "zhnmat_DIR: ${zhnmat_DIR}")
# message(STATUS "zhnmat_VERSION: ${zhnmat_VERSION}")
target_link_libraries(${CMAKE_PROJECT_NAME}_exe PUBLIC ${zhnmat_LIBS})
target_include_directories(${CMAKE_PROJECT_NAME}_exe PUBLIC ${zhnmat_INCLUDE_DIRS})
```
