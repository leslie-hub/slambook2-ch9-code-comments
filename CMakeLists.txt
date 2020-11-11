cmake_minimum_required(VERSION 2.8)

project(bundle_adjustment)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-O3 -std=c++11")

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

Find_Package(G2O REQUIRED)
Find_Package(Eigen3 REQUIRED)
Find_Package(Ceres REQUIRED)
Find_Package(Sophus REQUIRED)
Find_Package(CSparse REQUIRED)

SET(G2O_LIBS g2o_csparse_extension g2o_stuff g2o_core cxsparse)

include_directories(${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} ${CSPARSE_INCLUDE_DIR})

add_library(bal_common common.cpp)
add_executable(ceres_ba ceres_ba.cpp)
add_executable(g2o_ba g2o_ba.cpp)

target_link_libraries(ceres_ba ${CERES_LIBRARIES} bal_common)
target_link_libraries(g2o_ba ${G2O_LIBS} bal_common )
