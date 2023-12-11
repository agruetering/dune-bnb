add_subdirectory(bnb)
add_subdirectory(OCP)

add_executable(bnbtest BnBTest.cc)
target_compile_options(bnbtest PRIVATE -g -DOUTPUTDIR="")
target_link_libraries(bnbtest BB libpoly)