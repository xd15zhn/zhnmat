# 自用的矩阵运算库`zhnmat`使用说明
代码仓库  
<https://gitee.com/xd15zhn/zhnmat>  
<https://github.com/xd15zhn/zhnmat>  
安装教程建议参考  
[（二）EGE安装与配置 -CSDN博客](https://blog.csdn.net/qq_39151563/article/details/100161986)  

# 一个简单的例子
```cpp
//main.cpp
#include <iostream>
#include "zhnmat.hpp"
using namespace std;
using namespace zhnmat;
int main() {
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
set(CMAKE_PREFIX_PATH "C:/Users/xd15zhn/Documents/cpplibraries")
add_executable(${CMAKE_PROJECT_NAME} main.cpp)
find_package(zhnmat REQUIRED)
message(STATUS "zhnmat_VERSION: ${zhnmat_VERSION}")
message(STATUS "zhnmat_DIR: ${zhnmat_DIR}")
message(STATUS "zhnmat_LIBS: ${zhnmat_LIBS}")
message(STATUS "zhnmat_INCLUDE_DIRS: ${zhnmat_INCLUDE_DIRS}")
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${zhnmat_LIBS})
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${zhnmat_INCLUDE_DIRS})
```

# 其它功能

## 与MATLAB对应的初始化方法
MATLAB 写法：
```matlab
A = zeros(1, 5);
B = ones(1, 5);
C = 10 * ones(1, 5);
D = [-14,-6,-4; 0,13,0; -19,-20,-21]
```
zhnmat 写法：
```cpp
Mat A(1, 5);
Mat B(1, 5, 1);
Mat C(1, 5, 10);
Mat D(3, 3, vector<double>{-14,-6,-4, 0,13,0, -19,-20,-21});
```

## 矩阵操作
取出部分元素(切片)
```cpp
Mat B(4, 4);
Mat A = B.atr(0);  // 取出第0行元素
Mat A = B.atc(0);  // 取出第0列元素
Mat A = B(Rect(0, 1, 3, 2));  // 取出从第0行第1列开始的3行2列元素
```
串联
```cpp
Mat C1(1, 3, vecdble{1, 1, 0});
Mat C2 = eye(3);  // 单位矩阵
Mat C = VConcat(C1, C2);
```
插入与删除
```cpp
    vecdble data;
    for (size_t i = 0; i < 9; i++)
        data.push_back(rand());
    Mat A(3, 3, data);
    Mat ins(1, 3, vecdble{1, 1, 1});
    Mat B = Insertr(A, 1, ins);  // 将给定矩阵插入到第1行(下标从0开始)
    Mat C = Deletec(B, 2);  // 删除第2列(下标从0开始)
    Mat::_format.rowout = ";\n";
    cout << A << endl;
    cout << B << endl;
    cout << C << endl;
```

## 修改打印格式
下面的代码可以将矩阵输出成可用于 markdown 或 $\LaTeX$ 的格式。
```cpp
    Mat A(4, 3, vector<double>{
        0, 0, 0,
        -3.91895953790310, -5.23722226983912, -0.000310564114977309,
        0.0736634979806905, 0.0984424841255557, 0.00603621863958453,
        1.89333312742797, 0.364831438915809, 8.57639360143337e-05,
    });
    Mat::Set_Precision(6);  // 保留6位有效数字
    Mat::_format.prefix = "[";  // 输出矩阵之前先输出的字符串
    Mat::_format.rowin = " & ";  // 一行内每个元素之间插入的字符串
    Mat::_format.rowout = " \\\\\n";  // 两行之间输出的字符串
    Mat::_format.suffix = "]";  // 输出矩阵之后继续输出的字符串
    cout << A << endl;
```
输出为
```tex
[0 & 0 & 0 \\
-3.91896 & -5.23722 & -0.000310564 \\
0.0736635 & 0.0984425 & 0.00603622 \\
1.89333 & 0.364831 & 8.57639e-05]
```
