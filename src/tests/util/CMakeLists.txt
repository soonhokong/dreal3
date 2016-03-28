add_executable(box box.cpp)
target_link_libraries(box dreal ${EXTRA_LIBS})
add_test(box ${CMAKE_CURRENT_BINARY_DIR}/box)

add_executable(ibex_helper ibex_helper.cpp)
target_link_libraries(ibex_helper dreal ${EXTRA_LIBS})
add_test(ibex_helper ${CMAKE_CURRENT_BINARY_DIR}/ibex_helper)

add_executable(tls_box tls_box.cpp)
target_link_libraries(tls_box dreal ${EXTRA_LIBS})
add_test(tls_box ${CMAKE_CURRENT_BINARY_DIR}/tls_box)
