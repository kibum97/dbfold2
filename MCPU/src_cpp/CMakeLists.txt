cmake_minimum_required(VERSION 3.15) # Bumped for better modern dependency support
project(MCPU2)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wno-template-body)
endif()

# Specify the C++ standard 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# -------------------------------------------------------------------------
# DEPENDENCIES
# -------------------------------------------------------------------------

# Find MPI (Modern Target-based approach)
find_package(MPI REQUIRED)

# Find Eigen (Bumped to 3.4 to enable Matrix3/Vector3 shorthands)
# --- Eigen Dependency Logic ---
include(FetchContent)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
)
FetchContent_MakeAvailable(eigen)

# This creates the 'Eigen3::Eigen' target even if it's not on the system
if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

# Find other dependencies
find_package(fmt CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
# find_package(Eigen3 3.4.0 REQUIRED)

# -------------------------------------------------------------------------
# SOURCE DETECTION
# -------------------------------------------------------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# -------------------------------------------------------------------------
# TARGET CONFIGURATION
# -------------------------------------------------------------------------
add_executable(MCPU2 ${SOURCES})

# Target-specific include directories
# This is cleaner than global include_directories()
target_include_directories(MCPU2 PRIVATE 
    "include"
    "include/xtc"
)

# Link Libraries using Modern CMake Targets
# Note: MPI::MPI_CXX handles both include paths and library linking automatically
target_link_libraries(MCPU2 PRIVATE 
    MPI::MPI_CXX 
    Eigen3::Eigen 
    yaml-cpp::yaml-cpp 
    fmt::fmt
)

# -------------------------------------------------------------------------
# BENCHMARKS (Optional - uncomment as needed)
# -------------------------------------------------------------------------
# if(BUILD_BENCHMARKS)
#     add_executable(pairwise_distance_benchmark benchmark/pairwise_distance.cpp)
#     target_link_libraries(pairwise_distance_benchmark PRIVATE Eigen3::Eigen)
#     target_compile_options(pairwise_distance_benchmark PRIVATE -O3)
# endif()

# When changes made to CMakeLists.txt, cmake -B build -G Ninja -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug
# Regular compilation step: cmake --build build 2>&1 | tee build_errors.txt