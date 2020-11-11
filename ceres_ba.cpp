//
// Created by nnz on 2020/11/10.
//
#include <iostream>
#include <ceres/ceres.h>
#include "common.h"
#include "SnavelyReprojectionError.h"
using namespace std;
//求解最小二乘的函数
void SolveBA(BALProblem &bal_Problem);

int main(int argc,char **argv)
{
    if(argc!=2)
    {
        cout<<"usge : ./ceres_ba problem-16-22106-pre.txt"<<endl;
        return 1;
    }
    BALProblem bal_Problem(argv[1]);//读入数据
    bal_Problem.Normalize();//进行归一化处理
    bal_Problem.Perturb(0.1, 0.5, 0.5);//利用Perturb加入噪声
    bal_Problem.WriteToPLYFile("initial.ply");//将优化前的数据（相机和3d点） 保存在initial.ply文件中
    SolveBA(bal_Problem);//求解最小二乘问题
    bal_Problem.WriteToPLYFile("final.ply");//将优化后的数据（相机和3d点） 保存在final.ply文件中
    return 0;
}
//重点
void SolveBA(BALProblem &bal_Problem)
{
    const int point_block_size=bal_Problem.point_block_size();//
    const int camera_block_size=bal_Problem.camera_block_size();//
    //注意这里获得待优化系数首地址的时候要用mutable_points()和mutable_cameras()
    // 因为这两个函数指向的地址的内容是允许改变的（优化系数肯定要变的啦）
    double *points =bal_Problem.mutable_points();//获得待优化系数3d点 points指向3d点的首地址
    double *cameras = bal_Problem.mutable_cameras();//获得待优化系数相机 cameras指向相机的首地址
    const double *observations= bal_Problem.observations();//获得观测数据 observations指向观测数据的首地址
    ceres::Problem problem;
    //要用循环
    for(int i=0;i<bal_Problem.num_observations();i++)
    {
        ceres::CostFunction *cost_function;
        cost_function=SnavelyReprojectionError::Create(observations[2*i],observations[2*i+1]);
        ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);//核函数
        //bal_Problem.point_index()这返回的是一个地址指向索引号的首地址
        double *camera = cameras + camera_block_size*bal_Problem.camera_index()[i];
        double *point  = points  + point_block_size*bal_Problem.point_index()[i];
        //构建最小二乘问题
        problem.AddResidualBlock(
                // cons_function    ( new ceres::AutoDiffCostFunction<SnavelyReprojectionError,2,9,3>
                //                (
                //                 new SnavelyReprojectionError( observed_x,observed_y ))
                //                );
                cost_function,//代价函数
                loss_function,//核函数
                camera,//待优化的相机
                point//待优化的3d点
                );
    }
    //显示相关的信息
    // show some information here ...
    std::cout << "bal problem file loaded..." << std::endl;
    std::cout << "bal problem have " << bal_Problem.num_cameras() << " cameras and "
              << bal_Problem.num_points() << " points. " << std::endl;
    std::cout << "Forming " << bal_Problem.num_observations() << " observations. " << std::endl;

    std::cout << "Solving ceres BA ... " << endl;

    //配置求解器
    ceres::Solver::Options options;//这里有很多配置选项可以填
    options.linear_solver_type=ceres::LinearSolverType::SPARSE_SCHUR;//消元
    options.minimizer_progress_to_stdout= true;//输出到cout
    ceres::Solver::Summary summary;
    ceres::Solve(options,&problem,&summary);
    std::cout << summary.FullReport() << "\n";
}