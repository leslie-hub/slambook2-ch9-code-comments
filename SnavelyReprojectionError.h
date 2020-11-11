#ifndef SnavelyReprojection_H
#define SnavelyReprojection_H

#include <iostream>
#include "ceres/ceres.h"
#include "rotation.h"
//这是一个算投影残差的类  这里就像是 p138页的struct CURVE_FITTING_COST

class SnavelyReprojectionError {
public:
    //初始化 用有参构造获取观测数据
    SnavelyReprojectionError(double observation_x, double observation_y) : observed_x(observation_x),
                                                                           observed_y(observation_y) {}
//重载() //残差的计算
    template<typename T>
    bool operator()(const T *const camera,
                    const T *const point,
                    T *residuals) const {
        // camera[0,1,2] are the angle-axis rotation
        T predictions[2];//存放预测的像素点
        CamProjectionWithDistortion(camera, point, predictions);//得到去畸变后的像素，这个像素是预测值，也就是估计值 优化相机和 point
        //构造残差
        residuals[0] = predictions[0] - T(observed_x);
        residuals[1] = predictions[1] - T(observed_y);

        return true;
    }

    // camera : 9 dims array
    // [0-2] : angle-axis rotation 旋转
    // [3-5] : translateion 平移
    //  [6] focal length 焦距
    //  [7-8] second and forth order radial distortion 畸变系数
    //这里是认为 fx=fy  然后不考虑 cx cy
    // point : 3D location.世界坐标下的3d点
    // predictions : 2D predictions with center of the image plane.
    //去畸变的后的3d点在相机下的投影
    template<typename T>
    static inline bool CamProjectionWithDistortion(const T *camera, const T *point, T *predictions) {
        // Rodrigues' formula
        T p[3];
        AngleAxisRotatePoint(camera, point, p);//得到旋转后的点 p
        // camera[3,4,5] are the translation
        p[0] += camera[3];
        p[1] += camera[4];
        p[2] += camera[5];
        //p经过了旋转和平移 就是相机下的坐标
        //见书p102页 畸变部分
        // Compute the center fo distortion 不懂为什么这里有负号
        T xp = -p[0] / p[2];
        T yp = -p[1] / p[2];

        const T &l1 = camera[7];//second order radial distortion
        const T &l2 = camera[8];// forth order radial distortion

        T r2 = xp * xp + yp * yp;
        //下面的式子展开就是p102 式5.12
        // 这里认为 x方向的畸变和 y方向的畸变参数是一样的
        T distortion = T(1.0) + r2 * (l1 + l2 * r2);

        const T &focal = camera[6];
        predictions[0] = focal * distortion * xp;
        predictions[1] = focal * distortion * yp;

        return true;
    }

    //可以参考p139页的代码
    //将cost_function封装在类里面
    //这里是创建cost_function
    static ceres::CostFunction *Create(const double observed_x, const double observed_y)
    {

        //这里面残差是二维 待优化的是相机以及3d点 所以自动求导模板中的数字2,9,3分别是它们的维度
        return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 2, 9, 3>
                (
                        new SnavelyReprojectionError(observed_x, observed_y))
        );
    }

private:
    //观测数据
    double observed_x;
    double observed_y;
};

#endif // SnavelyReprojection.h
