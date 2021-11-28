#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <ctime>
#include "RANSAC.h"

using namespace std;

void RANSAC::read_data(const std::string &filepath, int data_size)
{
    this->points.clear();

    ifstream fp(filepath);
    string line;
    int line_num = 0;

    while (getline(fp, line) && line_num < data_size)
    {
        line_num++;

        double num[3];
        string str;
        std::stringstream ss(line);
        for (double &i: num)
        {
            getline(ss, str, ',');
            i = stod(str);
        }

        this->points.emplace_back(num);
    }
}

double RANSAC::solve_distance(RANSAC::Point3d M, RANSAC::Point3d P, RANSAC::Point3d N)
{
    double A = N.x;
    double B = N.y;
    double C = N.z;
    double D = -A * P.x - B * P.y - C * P.z;

    return fabs(A * M.x + B * M.y + C * M.z + D) / sqrt(A * A + B * B + C * C);
}

void RANSAC::solve_plane(RANSAC::Point3d A, RANSAC::Point3d B, RANSAC::Point3d C)
{
    //定义两个常量
    const double pi = 3.1415926535;
    RANSAC::Point3d N(0, 0, 1);

    //计算平面的单位法向量，即BC与BA的叉积
    RANSAC::Point3d Nx = (B - C) ^ (B - A);
    Nx.norm();

    //计算单位旋转向量与旋转角(范围0到pi)
    RANSAC::Point3d Nv = Nx ^ N;
    double angle = acos(Nx * N);

    //两个向量的夹角不大于pi/2,这里单独处理一下
    if (angle > pi / 2.0)
    {
        angle = pi - angle;
        Nv = Nv * (-1);
    }

    //利用私有变量把返回值保存出来,不太优美...
    this->plane_P = B;
    this->plane_Q = Quaternion(Nv * sin(angle / 2), cos(angle / 2));
    this->plane_N = Nx;
}

void RANSAC::ransac_core()
{
    //数据规模
    unsigned long size = this->points.size();
    //迭代的最大次数,每次得到更好的估计会优化iters的数值,默认10000
    int iters = 10000;
    //数据和模型之间可接受的差值,默认0.25
    double sigma = 0.15;
    //内点数目
    int pretotal = 0;
    //希望的得到正确模型的概率,默认0.99
    double per = 0.999;
    for (int i = 0; i < iters; i++)
    {
        //随机从数据中算则三个点去求解模型,这里取得点可能是重复的
        int index[3];
        srand((unsigned) time(NULL));
        index[0] = rand() % (size);
        index[1] = rand() % (size);
        index[2] = rand() % (size);
        //TODO 下面这样可能有问题，卡住了就注释掉，放弃保证点不重复
        //如果取的点是重复的,就重新取一次
        if (index[0] == index[1] || index[1] == index[2] || index[2] == index[0])
        {
            i--;
            continue;
        }
        solve_plane(this->points[index[0]], this->points[index[1]], this->points[index[2]]);

        //计算内点的数目
        int total_inlier = 0;
        for (int j = 0; j < size; j++)
        {
            if (solve_distance(points[j], this->plane_P, this->plane_N) < sigma)
            {
                total_inlier = total_inlier + 1;
            }
        }

        //判断当前的模型是否比之前估算的模型更好
        if (total_inlier > pretotal)
        {
            iters = int(log(1.0 - per) / log(1.0 - pow(total_inlier / size, 2)));
            pretotal = total_inlier;
        }

        //判断是否当前模型已经符合超过一半的点
        if (total_inlier > size / 2)
            break;
    }
}

void RANSAC::solve()
{
    ransac_core();
}

void RANSAC::test()
{
    //////////////////////////////////////////
    RANSAC::Point3d A(1, 2, 3);
    RANSAC::Point3d B(2, 3, 4);
    RANSAC::Point3d C = (A + B) ^ ((A - B)) * (A * B);
    cout << C.x << " " << C.y << " " << C.z << endl;
    C.norm();
    cout << C.x << " " << C.y << " " << C.z << endl;

    //////////////////////////////////////////
    RANSAC::Point3d D(1, 2, 3);
    RANSAC::Point3d E(0, 3, 4);
    RANSAC::Point3d F(1, 8, -10);
    solve_plane(D, E, F);
    cout << solve_distance(D, E, F) << endl;
    cout << plane_P.x << " " << plane_P.y << " " << plane_P.z << endl;
    cout << plane_Q.x << " " << plane_Q.y << " " << plane_Q.z << " " << plane_Q.w << endl;
    cout << plane_N.x << " " << plane_N.y << " " << plane_N.z << endl;

    //////////////////////////////////////////
    read_data("../pcltest.csv", 5000);
    ransac_core();
    cout << plane_P.x << " " << plane_P.y << " " << plane_P.z << endl;
    cout << plane_Q.x << " " << plane_Q.y << " " << plane_Q.z << " " << plane_Q.w << endl;
    cout << plane_N.x << " " << plane_N.y << " " << plane_N.z << endl;
}
