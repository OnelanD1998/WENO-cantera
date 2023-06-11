cmake_minimum_required(VERSION 3.15)

message(STATUS "The CMAKE_VERSION is ${CMAKE_VERSION}.")

project(singleMPIPrj)

find_package(MPI REQUIRED)

message(STATUS "The PROJECT_NAME is ${PROJECT_NAME}.")

include_directories("${MPI_CXX_INCLUDE_DIRS}" "/home/shuidong2/anaconda3/envs/libcantera/include" "/home/shuidong2/anaconda3/envs/libcantera/include" "/home/shuidong2/anaconda3/envs/libcantera/include/eigen3" "/home/conda/feedstock_root/build_artifacts/cantera-recipe_1667434613990/_build_env/include")

link_directories("/home/shuidong2/anaconda3/envs/libcantera/lib" "/home/shuidong2/anaconda3/envs/libcantera/lib" "/home/conda/feedstock_root/build_artifacts/cantera-recipe_1667434613990/_build_env/lib")

add_executable(${PROJECT_NAME} CtWenoFoam.C)

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES} cantera sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsollapackdense sundials_sunlinsollapackband fmt yaml-cpp lapack blas)