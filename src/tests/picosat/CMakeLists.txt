add_executable(picosatexample_01 example_01.cpp)
target_link_libraries(picosatexample_01 dreal ${EXTRA_LIBS} picosat)
add_test(picosatexample_01 ${CMAKE_CURRENT_BINARY_DIR}/picosatexample_01)
