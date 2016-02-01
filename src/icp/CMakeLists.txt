SET(DREAL_SRCS ${DREAL_SRCS}
  ${CMAKE_CURRENT_LIST_DIR}/icp.cpp
  ${CMAKE_CURRENT_LIST_DIR}/sat_icp.cpp
  ${CMAKE_CURRENT_LIST_DIR}/gsat_icp.cpp
  ${CMAKE_CURRENT_LIST_DIR}/picosat_wrapper.cpp
  ${CMAKE_CURRENT_LIST_DIR}/point_grid.cpp
  PARENT_SCOPE)

# add_library(dsolvers nra_solver.cpp)
# add_dependencies(dsolvers JSON)
# target_link_libraries(dsolvers dreal ${EXTRA_LIBS})
