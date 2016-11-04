find_package(Boost COMPONENTS system filesystem thread REQUIRED)
find_package(GSL REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Glog REQUIRED)
find_package(GFlags REQUIRED)
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/src/main)

file(GLOB test_src "*.cpp")

foreach(src_file ${test_src})
  get_filename_component(src_file_base ${src_file} NAME_WE)
  add_executable(${src_file_base} ${src_file})
  target_link_libraries(${src_file_base} spatioTemporalDecisions
    ${Boost_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${GSL_LIBRARIES}
    ${ARMADILLO_LIBRARIES} ${GTEST_LIBRARIES})
endforeach(src_file)
