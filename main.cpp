#include "RANSAC.h"

using namespace std;

/**
 * TODO 目标
 *  1. 实现ransac算法，对csv中的点云进行平面拟合
 *  2. 找到一套评估标准，保证拟合出来的平面是足够可靠的（如果拟合结果理论上不出问题那就不写这个了）
 */
int main(int argc, char **argv)
{
    RANSAC ransac;
    ransac.test();
}