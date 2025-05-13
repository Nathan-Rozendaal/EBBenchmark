cmake_minimum_required(VERSION 3.30)
project(EBBenchmark)

set(CMAKE_CXX_STANDARD 20)

add_executable(EBBenchmark main.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

# xtensor
find_package(xtensor REQUIRED)
target_include_directories(EBBenchmark PUBLIC ${xtensor_INCLUDE_DIRS})
target_link_libraries(EBBenchmark PUBLIC xtensor)

# aocommon
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include")

# casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)
target_include_directories(EBBenchmark PUBLIC ${CASACORE_INCLUDE_DIRS})
target_link_libraries(EBBenchmark PUBLIC ${CASACORE_LIBRARIES})

# everybeam
find_package(EveryBeam REQUIRED)
include_directories(${EVERYBEAM_INCLUDE_DIR})
target_link_libraries(EBBenchmark PUBLIC ${EVERYBEAM_LIB})
