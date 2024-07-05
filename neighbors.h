#ifndef neighbors_h
#define neighbors_h

#include "constant.h" // 一些常数

int neighs[SIZE][nei]; // 根据字符串长度定义数组

using namespace std; // 默认用库是std，这样可能造成冲突，比如用std里的函数当变量名

void get_neighbors(int neighs[SIZE][nei]);
void get_neighbors_square(int neighs[SIZE][nei]);
void get_neighbors_diced(int neighs[SIZE][nei]);
void get_neighbors_union_jack(int neighs[SIZE][nei]);

// int main()
// {
//     get_neighbors(neighs);
//     for (int i = 0; i < SIZE; i++)
//     {
//         for (int j = 0; j < nei; j++)
//         {
//             cout << neighs[i][j] << " ";
//         }
//         cout << endl;
//     }
//     return 0;
// }

void get_neighbors(int neighs[SIZE][nei])
{
    if (lattice == square)
    {
        get_neighbors_square(neighs);
    }
    else if (lattice == diced)
    {
        get_neighbors_diced(neighs);
    }
    else if (lattice == union_jack)
    {
        get_neighbors_union_jack(neighs);
    }
}

// Fills the neigbour tablet，填满相邻格子
void get_neighbors_square(int neighs[SIZE][nei]) // 相邻的有4个，这是个二维数组size*nei个元素
{
    constant();

    int UP = 0;
    int DOWN = 1;
    int LEFT = 2;
    int RIGHT = 3;

    int i, j;
    int u, d, r, l;

    for (i = 0; i < L; i++) // i列
    {
        for (j = 0; j < L; j++) // j行
        {
            // Get the (x,y) with periodic boundaries得到有周期边界条件的（x,y）
            u = j + 1 == L ? 0 : j + 1; // if j+1=L，u=0，else j=j+1
            d = j - 1 == -1 ? L - 1 : j - 1;
            r = i + 1 == L ? 0 : i + 1;
            l = i - 1 == -1 ? L - 1 : i - 1;

            //(x,y) to index notation and store in table，得到六个近邻格子的坐标
            neighs[i + j * L][UP] = i + u * L;
            neighs[i + j * L][DOWN] = i + d * L;
            neighs[i + j * L][LEFT] = l + j * L;
            neighs[i + j * L][RIGHT] = r + j * L;
        }
    }

    return;
}

// Fills the neigbour tablet，填满相邻格子
void get_neighbors_diced(int neighs[SIZE][nei]) // 相邻的有nei个，这是个二维数组size*nei个元素
{
    int UP = 0;
    int DOWN = 1;
    int LEFT = 2;
    int RIGHT = 3;
    int UP_LEFT = 4;
    int DOWN_RIGHT = 5;

    int i, j;
    int u, d, r, l;

    for (i = 0; i < L; i++) // i列
    {
        for (j = 0; j < L; j++) // j行
        {
            // Get the (x,y) with periodic boundaries得到有周期边界条件的（x,y）
            u = j + 1 == L ? 0 : j + 1; // if j+1=L，u=0，else j=j+1
            d = j - 1 == -1 ? L - 1 : j - 1;
            r = i + 1 == L ? 0 : i + 1;
            l = i - 1 == -1 ? L - 1 : i - 1;
            if ((i % 3) == (j % 3))
            {
                //(x,y) to index notation and store in table，得到六个近邻格子的坐标
                neighs[i + j * L][UP] = i + u * L;
                neighs[i + j * L][DOWN] = i + d * L;
                neighs[i + j * L][LEFT] = l + j * L;
                neighs[i + j * L][RIGHT] = r + j * L;
                neighs[i + j * L][UP_LEFT] = l + u * L;
                neighs[i + j * L][DOWN_RIGHT] = r + d * L;
            }
            else if (((i + 2) % 3) == ((j) % 3))
            {
                //(x,y) to index notation and store in table，得到三个近邻格子的坐标
                neighs[i + j * L][UP] = i + u * L;
                neighs[i + j * L][DOWN] = SIZE;
                neighs[i + j * L][LEFT] = l + j * L;
                neighs[i + j * L][RIGHT] = SIZE;
                neighs[i + j * L][UP_LEFT] = SIZE;
                neighs[i + j * L][DOWN_RIGHT] = r + d * L;
            }
            else
            {
                //(x,y) to index notation and store in table，得到三个近邻格子的坐标
                neighs[i + j * L][UP] = SIZE;
                neighs[i + j * L][DOWN] = i + d * L;
                neighs[i + j * L][LEFT] = SIZE;
                neighs[i + j * L][RIGHT] = r + j * L;
                neighs[i + j * L][UP_LEFT] = l + u * L;
                neighs[i + j * L][DOWN_RIGHT] = SIZE;
            }
            //     cout << neighs[i + j * L][UP] << " "
            //          << neighs[i + j * L][DOWN] << " "
            //          << neighs[i + j * L][LEFT] << " "
            //          << neighs[i + j * L][RIGHT] << " "
            //          << neighs[i + j * L][UP_LEFT] << " "
            //          << neighs[i + j * L][DOWN_RIGHT] << endl;
        }
    }

    return;
}

// Fills the neigbour tablet，填满相邻格子
void get_neighbors_union_jack(int neighs[SIZE][nei]) // 相邻的有四个，这是个二维数组size*4个元素
{
    int UP = 0;
    int DOWN = 1;
    int LEFT = 2;
    int RIGHT = 3;
    int UP_LEFT = 4;
    int UP_RIGHT = 5;
    int DOWN_RIGHT = 6;
    int DOWN_LEFT = 7;

    int i, j;
    int u, d, r, l;

    for (i = 0; i < L; i++) // i列
    {
        for (j = 0; j < L; j++) // j行
        {
            // Get the (x,y) with periodic boundaries得到有周期边界条件的（x,y）
            u = j + 1 == L ? 0 : j + 1;      // if j+1=L，u=0，else j=j+1
            d = j - 1 == -1 ? L - 1 : j - 1; // L=16，而格子等于15
            r = i + 1 == L ? 0 : i + 1;
            l = i - 1 == -1 ? L - 1 : i - 1;
            if (((i + j) % 2) == 0)
            {
                //(x,y) to index notation and store in table，得到四个近邻格子的坐标
                neighs[i + j * L][UP] = i + u * L;
                neighs[i + j * L][DOWN] = i + d * L;
                neighs[i + j * L][LEFT] = l + j * L;
                neighs[i + j * L][RIGHT] = r + j * L;
                neighs[i + j * L][UP_LEFT] = SIZE;
                neighs[i + j * L][UP_RIGHT] = SIZE;
                neighs[i + j * L][DOWN_RIGHT] = SIZE;
                neighs[i + j * L][DOWN_LEFT] = SIZE;
            }
            else
            {
                //(x,y) to index notation and store in table，得到八个近邻格子的坐标
                neighs[i + j * L][UP] = i + u * L;
                neighs[i + j * L][DOWN] = i + d * L;
                neighs[i + j * L][LEFT] = l + j * L;
                neighs[i + j * L][RIGHT] = r + j * L;
                neighs[i + j * L][UP_LEFT] = l + u * L;
                neighs[i + j * L][UP_RIGHT] = r + u * L;
                neighs[i + j * L][DOWN_RIGHT] = r + d * L;
                neighs[i + j * L][DOWN_LEFT] = l + d * L;
            }
        }
    }

    return;
}

// void get_neighbors_union_jack(int neighs[SIZE][nei]) // 相邻的有四个，这是个二维数组size*4个元素
// {
//     int UP = 0;
//     int DOWN = 1;
//     int LEFT = 2;
//     int RIGHT = 3;
//     int UP_LEFT = 4;
//     int UP_RIGHT = 5;
//     int DOWN_RIGHT = 6;
//     int DOWN_LEFT = 7;

//     int i, j;
//     int u, d, r, l;

//     for (i = 0; i < 2 * L; i += 2) // i列
//     {
//         for (j = 0; j < 2 * L; j += 2) // j行
//         {
//             // Get the (x,y) with periodic boundaries得到有周期边界条件的（x,y）
//             u = j + 1 == L ? 0 : j + 1;      // if j+1=L，u=0，else j=j+1
//             d = j - 1 == -1 ? L - 1 : j - 1; // L=16，而格子等于15
//             r = i + 1 == L ? 0 : i + 1;
//             l = i - 1 == -1 ? L - 1 : i - 1;
//             if (((i + j) % 2) == 0)
//             {
//                 //(x,y) to index notation and store in table，得到四个近邻格子的坐标
//                 neighs[i + j * L][UP] = i + u * L;
//                 neighs[i + j * L][DOWN] = i + d * L;
//                 neighs[i + j * L][LEFT] = l + j * L;
//                 neighs[i + j * L][RIGHT] = r + j * L;
//                 neighs[i + j * L][UP_LEFT] = SIZE;
//                 neighs[i + j * L][UP_RIGHT] = SIZE;
//                 neighs[i + j * L][DOWN_RIGHT] = SIZE;
//                 neighs[i + j * L][DOWN_LEFT] = SIZE;
//             }
//             else
//             {
//                 //(x,y) to index notation and store in table，得到八个近邻格子的坐标
//                 neighs[i + j * L][UP] = i + u * L;
//                 neighs[i + j * L][DOWN] = i + d * L;
//                 neighs[i + j * L][LEFT] = l + j * L;
//                 neighs[i + j * L][RIGHT] = r + j * L;
//                 neighs[i + j * L][UP_LEFT] = l + u * L;
//                 neighs[i + j * L][UP_RIGHT] = r + u * L;
//                 neighs[i + j * L][DOWN_RIGHT] = r + d * L;
//                 neighs[i + j * L][DOWN_LEFT] = l + d * L;
//             }
//         }
//     }

//     return;
// }

#endif