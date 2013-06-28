cmake_minimum_required(VERSION 2.4)
project("fortran_test" C Fortran)
#set(CMAKE_VERBOSE_MAKEFILE "ON")

add_executable(LOWERCASE csub.c fmain.f)
set_property(TARGET LOWERCASE PROPERTY "COMPILE_DEFINITIONS" "LOWERCASE")

add_executable(LOWERCASE_ csub.c fmain.f)
set_property(TARGET LOWERCASE_ PROPERTY COMPILE_DEFINITIONS "LOWERCASE_")

add_executable(LOWERCASE__ csub.c fmain.f)
set_property(TARGET LOWERCASE__ PROPERTY COMPILE_DEFINITIONS "LOWERCASE__")

add_executable(UPPERCASE csub.c fmain.f)
set_property(TARGET UPPERCASE PROPERTY COMPILE_DEFINITIONS "UPPERCASE")

add_executable(UPPERCASE_ csub.c fmain.f)
set_property(TARGET UPPERCASE_ PROPERTY COMPILE_DEFINITIONS "UPPERCASE_")

add_executable(UPPERCASE__ csub.c fmain.f)
set_property(TARGET UPPERCASE__ PROPERTY COMPILE_DEFINITIONS "UPPERCASE__")
