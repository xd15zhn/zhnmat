#include "zhnmat.hpp"
#include <iostream>
#include <random>
#include <ctime>
#define PRINT_NAME_VALUE(x)         std::cout<<#x<<": "<<x<<std::endl;

int main()
{
    using namespace std;
    using namespace zhnmat;
    Mat P1(3, 3, vector<double>{395.996938447599, 0, 355.3537673950195, 0, 395.996938447599, 237.5834827423096, 0, 0, 1});
    Mat p1inv = P1.inv();
    PRINT_NAME_VALUE(p1inv);
    cout << p1inv*P1 << endl;
    p1inv.set(0, 0, 1);
    Mat test1(eye(4));
    test1.Set_OutputFormat(USE_BRACKET | USE_SEMICOLON | WRAP_AROUND);
    PRINT_NAME_VALUE(test1);
    Mat test2(Gaussian_Kernel(1.52, 5));
    Mat test3(Convolution(test2, test1));
    test2.Set_OutputFormat(USE_BRACKET | USE_SEMICOLON | WRAP_AROUND);
    cout << test2 << endl;
    default_random_engine gen((unsigned int)time(0));  // 生成初始化种子
    normal_distribution<double> NormDis(0, 1);  // 正态分布
    uniform_real_distribution<double> UniFloatDis;  // [0,1]均匀分布
    Mat img(640, 480);
    for (int i=0; i<640; i++)
        for (int j=0; j<480; j++)
            img.set(i, j, NormDis(gen));
    Mat down(img, GENERATE_TYPE::DOWN_SAMPLE);
    cout << img.at(2, 4) << " " <<  down.at(1, 2) << endl;
    cout << img.at(638, 478) << " " << down.at(319, 239) << endl;
    return 0;
}
