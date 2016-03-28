add_executable(reduced_box_set test_reduced_box_set.cpp)
target_link_libraries(reduced_box_set dreal ${EXTRA_LIBS})
add_test(reduced_box_set ${CMAKE_CURRENT_BINARY_DIR}/reduced_box_set)

add_executable(gbox test_gbox.cpp)
target_link_libraries(gbox dreal ${EXTRA_LIBS})
add_test(gbox ${CMAKE_CURRENT_BINARY_DIR}/gbox)
