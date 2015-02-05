cmake_minimum_required(VERSION 2.8)
project(kuramoto_simulation)

if(CMAKE_COMPILER_IS_GNUCC)
    add_definitions(-std=c++11)
endif(CMAKE_COMPILER_IS_GNUCC)

find_package(Boost REQUIRED)
include_directories("${Boost_INCLUDE_DIRS}")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost COMPONENTS thread filesystem REQUIRED)

add_executable(kuramoto_simulation kuramoto_simulation.cpp)
target_link_libraries(kuramoto_simulation ${Boost_LIBRARIES})

