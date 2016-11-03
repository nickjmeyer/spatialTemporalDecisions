find_package(Boost COMPONENTS system filesystem thread REQUIRED)
find_package(GSL REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Glog REQUIRED)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(spatioTemporalDecisions agent.cpp allAgent.cpp bayesP.cpp
   calcCentrality.cpp copyParams.cpp data.cpp dataDepth.cpp
   densityEst.cpp features.cpp m1SpOptim.cpp mcmcEdgeToEdge2.cpp
   mcmcGravity2.cpp mcmcGravity.cpp model2EdgeToEdge.cpp
   model2GPowGDist.cpp model2GravityEDist.cpp model2GravityGDist.cpp
   model.cpp modelEDist.cpp modelGDist.cpp modelGravityGDist.cpp
   modelIntercept.cpp myopicAgent.cpp noTrtAgent.cpp paramBeta2.cpp
   paramBeta.cpp param.cpp paramEDist.cpp paramGDist.cpp
   paramGravityEDist.cpp paramGravityGDist.cpp paramGravPowGDist.cpp
   paramIntercept.cpp paramTrt.cpp proximalAgent.cpp rand.cpp
   randomAgent.cpp rankAgent.cpp runM1Mles.cpp runM1MlesMiss.cpp
   runM1MlesWns.cpp runM1MlesWnsMiss.cpp runner.cpp runStats.cpp
   sandbox.cpp settings.cpp sortMerge.cpp starts.cpp system.cpp
   test.cpp timer.cpp toyFeatures5.cpp tuneGen.cpp tuneGenWNS.cpp
   tuneParam.cpp utilities.cpp wnsFeatures3.cpp)

foreach(exec_var tuneGen tuneGenWNS)
  add_executable(${exec_var} ${exec_var}.cpp)
  target_link_libraries(${exec_var} spatioTemporalDecisions)
endforeach(exec_var)
