add_executable(picosatexample_01 example_01.cpp)
target_link_libraries(picosatexample_01 dreal ${EXTRA_LIBS} picosat)
add_test(picosatexample_01 ${CMAKE_CURRENT_BINARY_DIR}/picosatexample_01)

add_executable(picosatexample_02 example_02.cpp)
target_link_libraries(picosatexample_02 dreal ${EXTRA_LIBS} picosat)
add_test(picosatexample_02 ${CMAKE_CURRENT_BINARY_DIR}/picosatexample_02)

add_executable(picosatexample_03 example_03.cpp)
target_link_libraries(picosatexample_03 dreal ${EXTRA_LIBS} picosat)
add_test(picosatexample_03 ${CMAKE_CURRENT_BINARY_DIR}/picosatexample_03)
