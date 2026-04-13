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
# MPI TOGGLE OPTION
# -------------------------------------------------------------------------
option(USE_MPI "Enable MPI support" ON)
if (USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(-DUSE_MPI)
    message(STATUS "MPI support enabled.")
else()
    message(STATUS "MPI support disabled.")
endif()

# -------------------------------------------------------------------------
# DEPENDENCIES
# -------------------------------------------------------------------------
# Find Eigen (Bumped to 3.4 to enable Matrix3/Vector3 shorthands)
# --- Eigen Dependency Logic ---
include(FetchContent)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
)
FetchContent_MakeAvailable(eigen)

FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/archive/refs/tags/v3.12.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(json)

# This creates the 'Eigen3::Eigen' target even if it's not on the system
if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

# Find other dependencies
find_package(fmt CONFIG REQUIRED)
# find_package(yaml-cpp CONFIG REQUIRED)
# find_package(Eigen3 3.4.0 REQUIRED)

# -------------------------------------------------------------------------
# SOURCE DETECTION
# -------------------------------------------------------------------------
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
add_library(MCPU_lib STATIC ${LIB_SOURCES})

target_include_directories(MCPU_lib PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/xtc"
)

target_link_libraries(MCPU_lib PUBLIC 
    Eigen3::Eigen 
    nlohmann_json::nlohmann_json 
    fmt::fmt
)
if(USE_MPI)
    target_link_libraries(MCPU_lib PUBLIC MPI::MPI_CXX)
    target_compile_definitions(MCPU_lib PUBLIC USE_MPI)
endif()

# -------------------------------------------------------------------------
# TARGET CONFIGURATION
# -------------------------------------------------------------------------
add_executable(MCPU2 "apps/main.cpp")
target_link_libraries(MCPU2 PRIVATE MCPU_lib)

# -------------------------------------------------------------------------
# 3. UNIT TESTS (The verification in tests/)
# -------------------------------------------------------------------------
option(BUILD_TESTING "Build the test suite" ON)

if(BUILD_TESTING)
    enable_testing()
    
    file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
    
    if(TEST_SOURCES)
        add_executable(MCPU_tests ${TEST_SOURCES})
        # Link to the core logic + the Catch2 main() runner
        target_link_libraries(MCPU_tests PRIVATE MCPU_lib Catch2::Catch2WithMain)
        
        # This registers the binary with 'ctest'
        add_test(NAME AllTests COMMAND MCPU_tests)
    endif()
endif()

# When changes made to CMakeLists.txt, cmake -B build -G Ninja -DUSE_MPI=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug
# Regular compilation step: cmake --build build 2>&1 | tee build_errors.txt