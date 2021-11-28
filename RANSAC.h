#ifndef RANSAC_RANSAC_H
#define RANSAC_RANSAC_H

#include <string>
#include <vector>
#include <cmath>

//我认为我应该需要一个类来实现RANSAC算法，否则变量的传递会非常麻烦
class RANSAC
{
    struct Point3d
    {
        double x;
        double y;
        double z;

        explicit Point3d(const double *num)
        {
            this->x = num[0];
            this->y = num[1];
            this->z = num[2];
        }

        explicit Point3d(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z)
        {
        }

        //非常简单的归一化
        void norm()
        {
            double sum = this->x * this->x + this->y * this->y + this->z * this->z;
            sum = sqrt(sum);
            if (sum == 0)
                return;
            this->x = this->x / sum;
            this->y = this->y / sum;
            this->z = this->z / sum;
        }

        Point3d operator+(const Point3d M) const
        {
            double num[3];
            num[0] = this->x + M.x;
            num[1] = this->y + M.y;
            num[2] = this->z + M.z;

            Point3d point(num);
            return point;
        }

        Point3d operator-(const Point3d M) const
        {
            double num[3];
            num[0] = this->x - M.x;
            num[1] = this->y - M.y;
            num[2] = this->z - M.z;

            Point3d point(num);
            return point;
        }

        //数乘
        Point3d operator*(const double a) const
        {
            Point3d point;
            point.x = this->x * a;
            point.y = this->y * a;
            point.z = this->z * a;

            return point;
        }

        //内积
        double operator*(const Point3d M) const
        {
            return this->x * M.x + this->y * M.y + this->z * M.z;
        }

        //外积
        Point3d operator^(const Point3d M) const
        {
            double num[3];
            num[0] = this->y * M.z - this->z * M.y;
            num[1] = this->z * M.x - this->x * M.z;
            num[2] = this->x * M.y - this->y * M.x;

            Point3d point(num);
            return point;
        }
    };

    struct Quaternion
    {
        double x;
        double y;
        double z;
        double w;

        explicit Quaternion(Point3d point, double w)
        {
            this->x = point.x;
            this->y = point.y;
            this->z = point.z;
            this->w = w;
        }

        explicit Quaternion(double x = 0, double y = 0, double z = 0, double w = 0) : x(x), y(y), z(z), w(w)
        {
        }
    };

public:
    RANSAC() = default;

    //FIXME 发现一个很严重的问题
    // 点云中的点是一帧一帧连续计算出来的，相关性比较大
    // 如果直接用前1000个点计算，就计算不出正确的结果来，必须要多用一些点（比如5000个），或者随机取1000个点！
    //读取数据
    void read_data(const std::string &filepath, int data_size = 1000);

    //进行求解
    void solve();

    //测试代码
    void test();

private:
    /**
     * 求解点M到平面的距离
     * @param M 点M
     * @param P 平面上一点P
     * @param N 平面的法向量
     */
    static double solve_distance(Point3d M, Point3d P, Point3d N);

    /**
     * 根据三点求解平面方程
     * @param A 点A
     * @param B 点B
     * @param C 点C
     */
    void solve_plane(Point3d A, Point3d B, Point3d C);

    /**
     * RANSAC核心算法
     */
    void ransac_core();

private:
    //保存的点集
    std::vector<Point3d> points;
    //计算得到的平面上的一点
    Point3d plane_P;
    //计算得到的平面的四元数
    Quaternion plane_Q;
    //计算得到的平面的法向量
    Point3d plane_N;
};


#endif //RANSAC_RANSAC_H
