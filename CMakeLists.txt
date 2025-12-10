cmake_minimum_required(VERSION 3.12)
project(CEC2013LSGO)

set(CMAKE_CXX_STANDARD 23)
#O3 optimization
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)

set(INPUT_DATA_SRC "${CMAKE_SOURCE_DIR}/cdatafiles")
set(INPUT_DATA_DST "${CMAKE_BINARY_DIR}/cdatafiles")

file(GLOB_RECURSE INPUT_DATA_FILES
        CONFIGURE_DEPENDS
        "${INPUT_DATA_SRC}/*")

set(INPUT_DATA_STAMP "${CMAKE_BINARY_DIR}/.cdatafiles.stamp")

add_custom_command(
        OUTPUT "${INPUT_DATA_STAMP}"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${INPUT_DATA_DST}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${INPUT_DATA_SRC}" "${INPUT_DATA_DST}"
        COMMAND ${CMAKE_COMMAND} -E touch "${INPUT_DATA_STAMP}"
        DEPENDS ${INPUT_DATA_FILES}
        COMMENT "Copy `cdatafiles` to build directory"
        VERBATIM)
add_custom_target(copy_cdatafiles ALL DEPENDS "${INPUT_DATA_STAMP}")

file(GLOB BENCHMARK_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/benchmarks/*.cpp")
add_library(benchmarks Benchmarks.cpp eval_func.cpp ${BENCHMARK_FILES})

add_executable(demo demo.cpp)
add_dependencies(demo copy_cdatafiles)
target_link_libraries(demo benchmarks)

add_executable(demo2 demo2.cpp)
add_dependencies(demo2 copy_cdatafiles)
target_link_libraries(demo2 benchmarks)

add_executable(shade algorithms/shade.cpp algorithms/search_algorithm.cpp shade_main.cpp)
add_dependencies(shade copy_cdatafiles)
target_link_libraries(shade benchmarks OpenMP::OpenMP_CXX)

add_executable(lshade algorithms/lshade.cpp algorithms/search_algorithm.cpp main.cpp)
add_dependencies(lshade copy_cdatafiles)
target_link_libraries(lshade benchmarks OpenMP::OpenMP_CXX)

add_executable(uade algorithms/uade.cpp algorithms/search_algorithm.cpp uade_main.cpp)
add_dependencies(uade copy_cdatafiles)
target_link_libraries(uade benchmarks)

add_executable(uade_mt algorithms/uade_mt.cpp algorithms/search_algorithm.cpp uade_mt_main.cpp)
add_dependencies(uade_mt copy_cdatafiles)
target_link_libraries(uade_mt benchmarks OpenMP::OpenMP_CXX)

add_executable(uade_mt2 algorithms/uade_mt2.cpp algorithms/search_algorithm.cpp uade_mt_main.cpp)
add_dependencies(uade_mt2 copy_cdatafiles)
target_link_libraries(uade_mt2 benchmarks OpenMP::OpenMP_CXX)

add_executable(uade_mt3 algorithms/uade_mt3.cpp algorithms/search_algorithm.cpp uade_mt_main.cpp)
add_dependencies(uade_mt3 copy_cdatafiles)
target_link_libraries(uade_mt3 benchmarks OpenMP::OpenMP_CXX)

add_executable(uade_mt4 algorithms/uade_mt4.cpp algorithms/search_algorithm.cpp uade_mt_main.cpp Benchmarks.cpp eval_func.cpp ${BENCHMARK_FILES})
add_dependencies(uade_mt4 copy_cdatafiles)
target_link_libraries(uade_mt4 PRIVATE OpenMP::OpenMP_CXX)

add_executable(uade_mt_ssd_async_io algorithms/uade_mt_ssd_async_io.cpp algorithms/search_algorithm.cpp uade_mt_main.cpp Benchmarks.cpp eval_func.cpp ${BENCHMARK_FILES})
add_dependencies(uade_mt_ssd_async_io copy_cdatafiles)
target_link_libraries(uade_mt_ssd_async_io PRIVATE OpenMP::OpenMP_CXX)

add_executable(test test_main.cpp)