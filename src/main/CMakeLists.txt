find_package(Boost COMPONENTS system filesystem thread REQUIRED)
find_package(GSL REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Glog REQUIRED)
find_package(GFlags REQUIRED)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(spatioTemporalDecisions agent.cpp allAgent.cpp
  calcCentrality.cpp data.cpp dataDepth.cpp densityEst.cpp
  features.cpp m1SpOptim.cpp mcmcEdgeToEdge2.cpp mcmcGravity2.cpp
  mcmcGravity.cpp model2EdgeToEdge.cpp model2GPowGDist.cpp
  model2GravityEDist.cpp model2GravityGDist.cpp model.cpp
  modelEDist.cpp modelGDist.cpp modelGravityGDist.cpp
  modelIntercept.cpp myopicAgent.cpp noTrtAgent.cpp paramBeta2.cpp
  paramBeta.cpp param.cpp paramEDist.cpp paramGDist.cpp
  paramGravityEDist.cpp paramGravityGDist.cpp paramGravPowGDist.cpp
  paramIntercept.cpp paramTrt.cpp proximalAgent.cpp rand.cpp
  randomAgent.cpp rankAgent.cpp runner.cpp runStats.cpp settings.cpp
  sortMerge.cpp starts.cpp system.cpp timer.cpp toyFeatures5.cpp
  tuneParam.cpp utilities.cpp wnsFeatures3.cpp)

foreach(exec_var tuneGen tuneGenWNS runM1Mles runM1MlesMiss
    runM1MlesWns runM1MlesWnsMiss sandbox bayesP copyParams)
  add_executable(${exec_var} ${exec_var}.cpp)
  target_link_libraries(${exec_var} spatioTemporalDecisions ${Boost_LIBRARIES}
    ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${GSL_LIBRARIES}
    ${ARMADILLO_LIBRARIES})
endforeach(exec_var)
