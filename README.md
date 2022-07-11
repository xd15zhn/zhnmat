# zhnmat
自用的矩阵运算库。

## 一个简单的例子
```cpp
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

