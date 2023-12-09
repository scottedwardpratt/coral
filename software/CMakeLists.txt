cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
set(GSL_INCLUDE_DIR $ENV{GSL_INCLUDE_DIR})
project(coral CXX C)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(SYSTEM include ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include)

file(GLOB coral_src ./src/*/*.cc ./src/*/*/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_coral STATIC ${coral_src})
target_link_libraries(msu_coral ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})



