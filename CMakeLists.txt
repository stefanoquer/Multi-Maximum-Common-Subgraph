cmake_minimum_required(VERSION 3.16.0)
project(McSplit_Multigraph VERSION 0.1.0)

include(CTest)
enable_testing()
set(CMAKE_CXX_STANDARD 20)
find_package(Threads REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (WIN32)
  set(CMAKE_CXX_FLAGS "-Wall")
endif (WIN32)
if (!WIN32)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra")
endif (!WIN32)
#set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")

add_executable(McSplit_Multigraph mcsp.cpp graph.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(McSplit_Multigraph Threads::Threads)