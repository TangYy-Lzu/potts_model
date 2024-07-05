#include <iostream>
#include <sstream>

int main()
{
    const int B = 100;
    int array[B];

    // 初始化数组
    for (int i = 0; i < B; ++i)
    {
        array[i] = 0;
    }

// 并行区域
#pragma omp parallel num_threads(2)
    {
#pragma omp parallel for
        for (int b = 0; b < B; b++)
        {
            std::cout << b << std::endl;
        }
    }

    return 0;
}
