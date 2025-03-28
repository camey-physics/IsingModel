cmake_minimum_required(VERSION 3.10)
project(IsingModelProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(GSL REQUIRED)

# Define paths
set(SRC_DIR src)
set(INCLUDE_DIR include)
set(TEST_DIR tests)

# Compiler optimizations for g++ and Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
endif()

# Add library
add_library(IsingModel ${SRC_DIR}/IsingModel.cpp)
target_include_directories(IsingModel PUBLIC ${INCLUDE_DIR})

# Add simulated annealing executable
add_executable(simulatedAnnealing ${SRC_DIR}/simulatedAnnealing.cpp)
target_link_libraries(simulatedAnnealing PRIVATE IsingModel GSL::gsl GSL::gslcblas)
target_include_directories(simulatedAnnealing PRIVATE ${INCLUDE_DIR})

# Enable testing
enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

# Add test executable
add_executable(IsingModelTest ${TEST_DIR}/IsingModelTest.cpp)
target_link_libraries(IsingModelTest PRIVATE IsingModel GTest::GTest GTest::Main GSL::gsl GSL::gslcblas)

# Discover and register tests
gtest_discover_tests(IsingModelTest)