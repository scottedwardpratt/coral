cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
set(CXX $ENV{CXX})
project(gausstest_project CXX)

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/commonutils/software/include)

link_directories(${GITHOME_MSU}/coral/software/lib ${GITHOME_MSU}/commonutils/software/lib)
message("Libraries: commonutils coral ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

find_library(${GITHOME_MSU}/coral/software/lib msu_coral)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)

set(SOURCE_FILES ./gausstest.cc)
add_executable(gausstest ${SOURCE_FILES})

SET_SOURCE_FILES_PROPERTIES(
	gausstest.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a
 )
SET_SOURCE_FILES_PROPERTIES(
 	gausstest.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a
)

 
add_dependencies(gausstest extern_msu_coral extern_msu_commonutils)
target_link_libraries(gausstest msu_coral msu_commonutils ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
