cmake_minimum_required(VERSION 3.0)
project(spatioTemporalDecisions)

include(CTest)
set(CMAKE_CXX_FLAGS "--std=c++11")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)

include(GetGitRevisionDescription)

git_describe(GIT_DESCRIBE_RESULT)
add_definitions(-DGIT_DESCRIBE=${GIT_DESCRIBE_RESULT})

add_subdirectory(src/main)
