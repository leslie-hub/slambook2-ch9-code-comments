//
// Created by nnz on 2020/11/11.
//
#include <iostream>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/robust_kernel_impl.h>
#include "common.h"
#include "sophus/se3.hpp"
using namespace std;
using namespace Eigen;
using namespace Sophus;

//定义一个结构体表示相机 相机9维
struct PoseAndIntrinsics
{
    PoseAndIntrinsics() {}
    //显示构造进行赋值，把数据集的内容赋值过去
    explicit PoseAndIntrinsics(double *data_addr)
    {
        rotation=SO3d::exp(Sophus::Vector3d(data_addr[0],data_addr[1],data_addr[2]));
        translation=Sophus::Vector3d (data_addr[3],data_addr[4],data_addr[5]);
        focal = data_addr[6];
        k1=data_addr[7];
        k2=data_addr[8];
    }

    //set_to函数是将优化的系数放进内存
    void set_to(double *data_addr)
    {
        auto r=rotation.log();
        for (int i=0;i<3;i++)
        {
            data_addr[i]=r[i];
            data_addr[3+i]=translation[i];
        }
        data_addr[6]=focal;
        data_addr[7]=k1;
        data_addr[8]=k2;
    }

    SO3d rotation;//李群 旋转
    Sophus::Vector3d translation=Sophus::Vector3d::Zero();//平移
    double focal=0;//焦距
    double k1=0,k2=0;//畸变系数

};

//定义两个顶点 一个是 相机（PoseAndIntrinsics） 一个是路标点（3d点）
//顶点 ：相机（PoseAndIntrinsics）
class VertexPoseAndIntrinsics : public g2o::BaseVertex<9,PoseAndIntrinsics>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexPoseAndIntrinsics() {}

    //重置
    virtual void setToOriginImpl() override
    {
        _estimate=PoseAndIntrinsics();//给待优化系数赋上初始值
    }

    //更新
    virtual void oplusImpl(const double  *update) override
    {
        _estimate.rotation=SO3d::exp(Vector3d(update[0],update[1],update[2]))*_estimate.rotation;
        _estimate.translation+=Vector3d(update[3],update[4],update[5]);
        _estimate.focal+=update[6];
        _estimate.k1+=update[7];
        _estimate.k2+=update[8];
    }
    //根据路标点进行在相机上的投影
    Vector2d project(const Vector3d point)
    {
        Sophus::Vector3d pc=_estimate.rotation*point+_estimate.translation;
        pc=-pc/pc[2];//把相机坐标归一化
        double r2=pc.squaredNorm();//相机坐标归一化后的模的平方
        double distortion = 1.0+(_estimate.k1+_estimate.k2*r2)*r2;
        return Sophus::Vector2d (_estimate.focal*distortion*pc[0],_estimate.focal*distortion*pc[1]);
    }

    //存盘和读盘 ： 留空
    virtual bool read(istream &in) {}

    virtual bool write(ostream &out) const {}
};

//顶点 ：路标点（3d点）
class VertexPoint : public g2o::BaseVertex<3,Sophus::Vector3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexPoint() {}    
    
    //重置
    virtual void setToOriginImpl() override
    {
        _estimate=Sophus::Vector3d (0,0,0);//待优化系数的初始化
    }

    //更新
    virtual void oplusImpl(const double  *update) override
    {
        _estimate+=Sophus::Vector3d (update[0],update[1],update[2]);
    }
    
    //存盘和读盘 ： 留空
    virtual bool read(istream &in) {}

    virtual bool write(ostream &out) const {}
};

//定义边 这里面边的定义比前几章的要简单很多
class EdgeProjection : public g2o::BaseBinaryEdge<2,Sophus::Vector2d,VertexPoseAndIntrinsics,VertexPoint>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    //计算残差
    virtual void computeError() override
    {
        auto v0 = (VertexPoseAndIntrinsics *) _vertices[0];
        auto v1 = (VertexPoint *) _vertices[1];
        auto proj= v0 -> project(v1->estimate());
        _error=proj-_measurement;
    }
    //存盘和读盘 ： 留空
    virtual bool read(istream &in) {}

    virtual bool write(ostream &out) const {}
};

void SolveBA(BALProblem &bal_problem);
int main(int argc,char **argv)
{
    if(argc!=2)
    {
        cout<<"usge : ./g2o_ba problem-16-22106-pre.txt"<<endl;
        return 1;
    }
    BALProblem bal_problem(argv[1]);//传入数据
    bal_problem.Normalize();//对数据进行归一化
    bal_problem.Perturb(0.1,0.5,0.5);//给数据加上噪声（相机旋转、相机平移、路标点）
    bal_problem.WriteToPLYFile("initial_g20");
    SolveBA(bal_problem);//求解BA
    bal_problem.WriteToPLYFile("final_g20");
    return 0;
}
void SolveBA(BALProblem &bal_problem)
{
    //获得 相机和点的维度
    const int camera_block_size=bal_problem.camera_block_size();
    const int point_block_size=bal_problem.point_block_size();
    //获得相机和点各自参数的首地址
    double *cameras=bal_problem.mutable_cameras();
    double *points=bal_problem.mutable_points();
    //获得观测值的首地址
    const double *observations=bal_problem.observations();
    //构建图优化
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<9,3>> BlockSolverType;//两个顶点的维度
    typedef g2o::LinearSolverCSparse<BlockSolverType::PoseMatrixType> LinearSolverTpe;
    //使用LM方法
    auto solver = new g2o::OptimizationAlgorithmLevenberg(
            g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverTpe>())
            );
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);//设置求解器
    optimizer.setVerbose(true);//打开调试输出

    //加入顶点
    //因为顶点有很多个，所以需要容器
    //容器 vertex_pose_intrinsics 和 vertex_points存放两顶点的地址
    vector<VertexPoseAndIntrinsics *> vertex_pose_intrinsics;
    vector<VertexPoint *> vertex_points;

    for(int i=0;i<bal_problem.num_cameras();i++)
    {
        VertexPoseAndIntrinsics *v = new VertexPoseAndIntrinsics();
        double *camera=cameras+camera_block_size*i;//获得每个相机的首地址
        v->setId(i);//设置编号
        v->setEstimate(PoseAndIntrinsics(camera));//传入待优化的系数 此处为相机
        optimizer.addVertex(v);//加入顶点
        vertex_pose_intrinsics.push_back(v);
    }

    for(int i=0; i<bal_problem.num_points();i++)
    {
        VertexPoint *v=new VertexPoint();//获得每个路标点的首地址
        double *point=points+point_block_size*i;
        v->setId(i+bal_problem.num_cameras());
        v->setEstimate(Sophus::Vector3d(point[0],point[1],point[3]));//传入待优化的系数 此处为路标点
        //g20需要手动设置待边缘化（Marg）的点
        v->setMarginalized(true);//设置边缘化
        optimizer.addVertex(v);//加入顶点
        vertex_points.push_back(v);//将顶点一个一个放回到容器里面
    }

    //加入边
    for(int i=0;i<bal_problem.num_observations();i++)
    {
        EdgeProjection *edge = new EdgeProjection;
        edge->setId(i);//设置编号
        edge->setVertex(0,vertex_pose_intrinsics[bal_problem.point_index()[i]]);//加入顶点
        edge->setVertex(1,vertex_points[bal_problem.camera_index()[i]]);//加入顶点
        edge->setMeasurement(Sophus::Vector2d(observations[2*i+0],observations[2*i+1]));//设置观测数据
        edge->setInformation(Eigen::Matrix2d::Identity());//设置信息矩阵
        edge->setRobustKernel(new g2o::RobustKernelHuber());//设置核函数
        optimizer.addEdge(edge);//加入边
    }
    optimizer.initializeOptimization();
    optimizer.optimize(40);
    //优化后在存到内从中去
    for( int i=0;i<bal_problem.num_cameras();i++)
    {
        double *camera =cameras + camera_block_size*i;
        auto vertex = vertex_pose_intrinsics[i];//把优化后的顶点地址给vertex
        auto estimate = vertex->estimate();//这样estimate就指向了优化后的相机结构体 此时estimate本质上指向了 相机结构体
        estimate.set_to(camera);//这样camera就指向了优化后的相机（利用了相机结构体的set_to()函数）
    }
    for (int i = 0; i < bal_problem.num_points(); ++i) {
        double *point = points + point_block_size * i;
        auto vertex = vertex_points[i];
        for (int k = 0; k < 3; ++k) point[k] = vertex->estimate()[k];
    }
    //经过上面的两个循环后，原来待优化的系数就被优化完毕了，并且优化后的系数 还是存放在 bal_problem.mutable_cameras()和 bal_problem.mutable_points() 所对应的地址中

}