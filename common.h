#pragma once

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


//注： const double *t //说明t指向的内容不可修改 而t的指向可以修改
//在归一化处理的时候用的就是mutable_points() 和 mutable_cameras() 因为归一化改变了参数
//同理  加入噪声函数也是用了mutable_points() 和 mutable_cameras()

/// 从文件读入BAL dataset
class BALProblem {
public:
    /// load bal data from text file
    //explict用显示构造
    explicit BALProblem(const std::string &filename, bool use_quaternions = false);

//析构函数 释放内存
    ~BALProblem() {
        delete[] point_index_;
        delete[] camera_index_;
        delete[] observations_;
        delete[] parameters_;
    }

    /// save results to text file
    void WriteToFile(const std::string &filename) const;

    /// save results to ply pointcloud 保存为点云文件
    void WriteToPLYFile(const std::string &filename) const;

    //归一化
    void Normalize();

   //加入给相机和点噪声
    void Perturb(const double rotation_sigma,
                 const double translation_sigma,
                 const double point_sigma);

    //获得相机的维度，本实验用的是9维
    int camera_block_size() const { return use_quaternions_ ? 10 : 9; }

    //获得点的维度  奇怪，3d点当然是三维咯
    int point_block_size() const { return 3; }

    //获得相机个数
    int num_cameras() const { return num_cameras_; }

    //获得点的个数
    int num_points() const { return num_points_; }

    //获得观测的个数
    int num_observations() const { return num_observations_; }

    //获得参数的个数
    int num_parameters() const { return num_parameters_; }

    //获得点的索引号
    const int *point_index() const { return point_index_; }

    //获得相机的索引号
    const int *camera_index() const { return camera_index_; }

    //获得观测的数据
    const double *observations() const { return observations_; }

    //获得参数
    const double *parameters() const { return parameters_; }

    //获得相机参数
    const double *cameras() const { return parameters_; }

    //获得点的参数
    // parameters_ + camera_block_size() * num_cameras_因为上面是相机的参数，因此要隔开
    const double *points() const { return parameters_ + camera_block_size() * num_cameras_; }

    /// 待优化camera参数的起始地址
    double *mutable_cameras() { return parameters_; }

    /// 待优化点的参数的起始地址
    double *mutable_points() { return parameters_ + camera_block_size() * num_cameras_; }

    //得到不同索引号下的相机参数的地址
    double *mutable_camera_for_observation(int i) {
        return mutable_cameras() + camera_index_[i] * camera_block_size();
    }

    //得到不同索引号下的点的参数的地址
    double *mutable_point_for_observation(int i) {
        return mutable_points() + point_index_[i] * point_block_size();
    }

    //得到不同索引号下的相机参数的地址
    const double *camera_for_observation(int i) const {
        return cameras() + camera_index_[i] * camera_block_size();
    }

    //得到不同索引号下的点的参数的地址
    const double *point_for_observation(int i) const {
        return points() + point_index_[i] * point_block_size();
    }

private:
    //个人理解是把相机的中心取出来，并用绿色的点描述
    //该函数还有一个作用就是，如果相机旋转用四元数描述，它会转为旋转向量
    void CameraToAngelAxisAndCenter(const double *camera,
                                    double *angle_axis,
                                    double *center) const;

    //将相机的中心恢复？个人觉得是把相机中的平移恢复
    //这个函数还有一个作用就是，如果数据集用四元数 由于我们使用CameraToAngelAxisAndCenter()函数把四元数变成了旋转向量
    //因此利用该函数将旋转向量还原到四元数
    void AngleAxisAndCenterToCamera(const double *angle_axis,
                                    const double *center,
                                    double *camera) const;

    int num_cameras_;//相机个数
    int num_points_;//3d点个数
    int num_observations_;//观测值个数
    int num_parameters_;//相机加3d点参数的数量
    bool use_quaternions_;//是否用四元数

    int *point_index_;      // 每个observation对应的point index
    int *camera_index_;     // 每个observation对应的camera index
    double *observations_;//指向观测值的指针
    double *parameters_;//指向参数的的指针
};
