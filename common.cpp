
//BAL数据集格式
// 相机数量 点的数量 观测数据数量 #第一行
// 相机索引编号 点索引编号 #第二行
//  相机_1的内容        #
//      .
//      .
//      .
//  相机_x的内容        #
//  点_1的内容        #
//      .
//      .
//      .
//  点_x的内容        #

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "common.h"
#include "rotation.h"
#include "random.h"

typedef Eigen::Map<Eigen::VectorXd> VectorRef; //向量的维度可变
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;

//将文件流入fptr 并判断文件有没有流入
template<typename T>
void FscanfOrDie(FILE *fptr, const char *format, T *value) {
    int num_scanned = fscanf(fptr, format, value);//返回1才说明文件
    if (num_scanned != 1)
        std::cerr << "Invalid UW data file. ";
}

// 给3d点/三维向量加噪声
void PerturbPoint3(const double sigma, double *point)
{
    for (int i = 0; i < 3; ++i)
        point[i] += RandNormal() * sigma;
}

//找到一组数据的中值
double Median(std::vector<double> *data) {
    int n = data->size();
    std::vector<double>::iterator mid_point = data->begin() + n / 2;
    std::nth_element(data->begin(), mid_point, data->end());
    return *mid_point;
}

//读取BAL数据集
BALProblem::BALProblem(const std::string &filename, bool use_quaternions) {
    FILE *fptr = fopen(filename.c_str(), "r");//以读的形式

    if (fptr == NULL) { //判断文件是否有效
        std::cerr << "Error: unable to open file " << filename;
        return;
    };

    // This wil die horribly on invalid files. Them's the breaks.
    FscanfOrDie(fptr, "%d", &num_cameras_);//得到相机的个数
    FscanfOrDie(fptr, "%d", &num_points_);//得到3d点的个数
    FscanfOrDie(fptr, "%d", &num_observations_);//得到观测数据的个数

    //输出一下 相机的个数  3d点的个数  观测数据的个数
    std::cout << "Header: " << num_cameras_
              << " " << num_points_
              << " " << num_observations_;


    //point_index_是对应3d点的索引编号
    //camera_index_是对应相机的索引编号
    //observations_ 是对应观测值的索引编号

    point_index_ = new int[num_observations_];//给point_index_开辟num_observations_大小的空间
    camera_index_ = new int[num_observations_];//camera_index_开辟num_observations_大小的空间
    observations_ = new double[2 * num_observations_];//给observations_开辟2*num_observations_大小的空间
                                                      // 因为observations_是二维的像素坐标

    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;//算所有参数的数量 相机9维 3d点三维
    parameters_ = new double[num_parameters_];//给 parameters_开辟num_parameters_大小的空间
                                              // parameters_这就是存放相机和3d点的参数
    //按格式传入数据
    for (int i = 0; i < num_observations_; ++i) {
        FscanfOrDie(fptr, "%d", camera_index_ + i);//把相机索引编号传入camera_index_
        FscanfOrDie(fptr, "%d", point_index_ + i);//把3d点索引编号传入point_index_
        for (int j = 0; j < 2; ++j) {
            FscanfOrDie(fptr, "%lf", observations_ + 2 * i + j);//把观测数据传入observations_
        }
    }
    //上面的数据传入完毕了 接下来是 相机和3d点的参数了 parameters_
    for (int i = 0; i < num_parameters_; ++i) {
        FscanfOrDie(fptr, "%lf", parameters_ + i);//
    }

    fclose(fptr);

    use_quaternions_ = use_quaternions;
    //如果数据集使用四元数 把旋转向量转为四元数
    if (use_quaternions) {
        // Switch the angle-axis rotations to quaternions.
        num_parameters_ = 10 * num_cameras_ + 3 * num_points_;//参数的个数 数据集使用四元数的话，相机维度是10维
        double *quaternion_parameters = new double[num_parameters_];
        double *original_cursor = parameters_;//original_cursor指向parameters_
        double *quaternion_cursor = quaternion_parameters;//quaternion_cursor指向和parameters_大小一样的空间
        for (int i = 0; i < num_cameras_; ++i) {
            AngleAxisToQuaternion(original_cursor, quaternion_cursor);//旋转向量转化为四元数
            quaternion_cursor += 4;//移动到下一个四元数
            original_cursor += 3;//移动到下一个旋转向量
            for (int j = 4; j < 10; ++j) {
                *quaternion_cursor++ = *original_cursor++;//旋转后面的数据不改动
            }
        }
        //把上面的相机参数弄完了，直接让二级指针*quaternion_cursor指向剩余的3d点
        // Copy the rest of the points.
        for (int i = 0; i < 3 * num_points_; ++i) {
            *quaternion_cursor++ = *original_cursor++;
        }
        // Swap in the quaternion parameters.
        delete[]parameters_;//释放parameters_
        parameters_ = quaternion_parameters;
    }
}

//好像没用到，就是将数据写进一个文件里
void BALProblem::WriteToFile(const std::string &filename) const {
    FILE *fptr = fopen(filename.c_str(), "w");

    if (fptr == NULL) {
        std::cerr << "Error: unable to open file " << filename;
        return;
    }

    fprintf(fptr, "%d %d %d %d\n", num_cameras_, num_cameras_, num_points_, num_observations_);

    for (int i = 0; i < num_observations_; ++i) {
        fprintf(fptr, "%d %d", camera_index_[i], point_index_[i]);
        for (int j = 0; j < 2; ++j) {
            fprintf(fptr, " %g", observations_[2 * i + j]);
        }
        fprintf(fptr, "\n");
    }

    for (int i = 0; i < num_cameras(); ++i) {
        double angleaxis[9];
        if (use_quaternions_) {
            //OutPut in angle-axis format.
            QuaternionToAngleAxis(parameters_ + 10 * i, angleaxis);
            memcpy(angleaxis + 3, parameters_ + 10 * i + 4, 6 * sizeof(double));
        } else {
            memcpy(angleaxis, parameters_ + 9 * i, 9 * sizeof(double));
        }
        for (int j = 0; j < 9; ++j) {
            fprintf(fptr, "%.16g\n", angleaxis[j]);
        }
    }

    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i) {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j) {
            fprintf(fptr, "%.16g\n", point[j]);
        }
    }

    fclose(fptr);
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
//保存ply点云文件
void BALProblem::WriteToPLYFile(const std::string &filename) const {
    std::ofstream of(filename.c_str());

    of << "ply"
       << '\n' << "format ascii 1.0"
       << '\n' << "element vertex " << num_cameras_ + num_points_
       << '\n' << "property float x"
       << '\n' << "property float y"
       << '\n' << "property float z"
       << '\n' << "property uchar red"
       << '\n' << "property uchar green"
       << '\n' << "property uchar blue"
       << '\n' << "end_header" << std::endl;

    //用绿色点描述相机中心
    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras(); ++i) {
        const double *camera = cameras() + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        of << center[0] << ' ' << center[1] << ' ' << center[2]
           << " 0 255 0" << '\n';
    }

    //3d点用白色的点描述
    // Export the structure (i.e. 3D Points) as white points.
    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i) {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j) {
            of << point[j] << ' ';
        }
        of << " 255 255 255\n";
    }
    of.close();
}

//个人理解是把相机的中心取出来，并用绿色的点描述
//该函数还有一个作用就是，如果相机旋转用四元数描述，它会转为旋转向量
void BALProblem::CameraToAngelAxisAndCenter(const double *camera,
                                            double *angle_axis,
                                            double *center) const {
    VectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_) {
        QuaternionToAngleAxis(camera, angle_axis);
    } else {
        angle_axis_ref = ConstVectorRef(camera, 3);//拿到旋转向量
    }

    // c = -R'*t 把平移的向量进行反向旋转并乘-
    Eigen::VectorXd inverse_rotation = -angle_axis_ref;
    AngleAxisRotatePoint(inverse_rotation.data(),
                         camera + camera_block_size() - 6,//这是【平移
                         center);
    VectorRef(center, 3) *= -1.0;
}


//将相机的中心恢复？个人觉得是把相机中的平移恢复
//这个函数还有一个作用就是，如果数据集用四元数 由于我们使用CameraToAngelAxisAndCenter()函数把四元数变成了旋转向量
//因此利用该函数将旋转向量还原到四元数
void BALProblem::AngleAxisAndCenterToCamera(const double *angle_axis,
                                            const double *center,
                                            double *camera) const {
    ConstVectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_) {
        AngleAxisToQuaternion(angle_axis, camera);
    } else {
        VectorRef(camera, 3) = angle_axis_ref;
    }

    // t = -R * c
    AngleAxisRotatePoint(angle_axis, center, camera + camera_block_size() - 6);
    VectorRef(camera + camera_block_size() - 6, 3) *= -1.0;
}
//对数据进行归一化（3d点 和 相机 【旋转和平移】）
void BALProblem::Normalize() {
    // Compute the marginal median of the geometry 计算几何图形的边缘中值
    std::vector<double> tmp(num_points_);
    Eigen::Vector3d median;
    //给3d点进行归一化
    double *points = mutable_points();//points指向3d点的首地址 3d点是待优化的点
    //这个循环就是找到中间的点
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < num_points_; ++j) {
            tmp[j] = points[3 * j + i];
        }
        median(i) = Median(&tmp);//  第一次找到点的x 第二次是找到中间点的y 第三次是找到中间点的z
    }

    for (int i = 0; i < num_points_; ++i) {
        VectorRef point(points + 3 * i, 3);
        tmp[i] = (point - median).lpNorm<1>();//去中心点后求 1型范数，就是向量的求模
                                              // 每一个tmp都放着去中心后的点的1型范数

    }

    const double median_absolute_deviation = Median(&tmp);//找到这些去了中心点的点的范数中值

    // Scale so that the median absolute deviation of the resulting
    // reconstruction is 100

    const double scale = 100.0 / median_absolute_deviation;//得到缩放的比例

    //这里才是真正实现了3d点的归一化
    // X = scale * (X - median)
    for (int i = 0; i < num_points_; ++i) {
        VectorRef point(points + 3 * i, 3);
        point = scale * (point - median);
    }
    //给相机归一化
    double *cameras = mutable_cameras();
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras_; ++i) {
        double *camera = cameras + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);//此处center为相机的中心
        // center = scale * (center - median)   这个是相机中心去掉3d点的中心并进行尺度缩放得到最后的相机中心
        VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
        AngleAxisAndCenterToCamera(angle_axis, center, camera);//最后再将归一化后的相机中心还原回相机
                                                               // 主要是还原t 因为这时候t是相机经过归一化处理的到的
    }
}

//给数据加上噪声
void BALProblem::Perturb(const double rotation_sigma,//旋转噪声
                         const double translation_sigma,//平移噪声
                         const double point_sigma) {//3d点的噪声
    //assert的作用是现计算表达式 expression ，如果其值为假（即为0），那么它先向stderr打印一条出错信息，然后通过调用 abort 来终止程序运行。
    assert(point_sigma >= 0.0);
    assert(rotation_sigma >= 0.0);
    assert(translation_sigma >= 0.0);

    double *points = mutable_points();
    if (point_sigma > 0) {
        //用循环给每个点加上噪声
        for (int i = 0; i < num_points_; ++i) {
            PerturbPoint3(point_sigma, points + 3 * i);
        }
    }
    //给相机中的平移和旋转加上噪声
    for (int i = 0; i < num_cameras_; ++i) {
        double *camera = mutable_cameras() + camera_block_size() * i;

        double angle_axis[3];
        double center[3];
        // Perturb in the rotation of the camera in the angle-axis
        // representation
        CameraToAngelAxisAndCenter(camera, angle_axis, center);//用这个就是确保用旋转向量描述旋转（如果数据集用的四元数，就会转换为旋转向量）
    //给相机的旋转向量加上噪声
        if (rotation_sigma > 0.0) {
            PerturbPoint3(rotation_sigma, angle_axis);
        }
        AngleAxisAndCenterToCamera(angle_axis, center, camera);//之后把t恢复 (如果数据集用的四元数，也会还原为四元数)
    //给相机的平移加上噪声
        if (translation_sigma > 0.0)
            PerturbPoint3(translation_sigma, camera + camera_block_size() - 6);
    }
}
