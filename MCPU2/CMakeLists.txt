cmake_minimum_required(VERSION 3.10)
project(MCPU2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# Include directories
include_directories(include)
include_directories(include/xtc)

# Find MPI
set(MPI_C_COMPILER mpicc)
set(MPI_CXX_COMPILER mpicxx)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Auto detect source files
file (GLOB UTILS_SOURCES "src/utils/*.cpp")
# file (GLOB_RECURSE ALL_SOURCES "src/*.cpp")

# Source files
set(SOURCES
    main.cpp
    ${UTILS_SOURCES}
    "src/mpi_utils.cpp"
    "src/logger.cpp"
    "src/read_config.cpp"
    "src/pdb_parser.cpp"
    "src/protein_structure.cpp"
    "src/mcpu_forcefield.cpp"
    "src/system_mcpu.cpp"
    "src/context_mcpu.cpp"
    "src/cell_list.cpp"
)

# Add the executable
add_executable(MCPU2 ${SOURCES})

# Find and link libraries
find_package(fmt CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(MCPU2 PRIVATE ${MPI_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp fmt::fmt)

# For benchmark
add_executable(pairwise_distance_benchmark benchmark/pairwise_distance.cpp)
add_executable(eigen_vector_benchmark benchmark/eigen_vector_1d_vec_benchmark.cpp)
target_link_libraries(pairwise_distance_benchmark PRIVATE Eigen3::Eigen)
target_compile_options(pairwise_distance_benchmark PRIVATE -O3)
target_link_libraries(eigen_vector_benchmark PRIVATE Eigen3::Eigen)
target_compile_options(eigen_vector_benchmark PRIVATE -O3)