cmake_minimum_required(VERSION 3.10)
project(MCPU2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)
include_directories(include/xtc)

# Find MPI
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Source files
set(SOURCES
    main.cpp
    src/read_config.cpp
    src/load_pdb.cpp
    src/mpi_util.cpp
    src/logger.cpp
    src/topology.cpp
    src/context.cpp
    src/cell_list.cpp
    src/system.cpp
    src/mcpu_forcefield.cpp
    src/integrator.cpp
    src/reporter.cpp
    src/xtc/xdrfile.c
    src/xtc/xdrfile_xtc.c
    src/utils/geometry.cpp
    src/utils/hbond.cpp
)

# Add the executable
add_executable(MCPU2 ${SOURCES})

# Find and link libraries
find_package(fmt CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(MCPU2 PRIVATE ${MPI_LIBRARIES} Eigen3::Eigen yaml-cpp::yaml-cpp fmt::fmt)