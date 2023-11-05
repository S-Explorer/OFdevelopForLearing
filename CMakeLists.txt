cmake_minimum_required(VERSION 3.20)

if (DEFINED ENV{WM_PROJECT})
    message("*********** hello OpenFOAMer ***********
   OpenFOAM version : $ENV{WM_PROJECT}-$ENV{WM_PROJECT_VERSION}
   OpenFOAM install dir: $ENV{WM_PROJECT_DIR}")
else()
    message(FATAL_ERROR "can't find OpenFOAM envs !")
endif()

project(CustomFOAM)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(OpenFOAM_VERSION $ENV{WM_PROJECT_VERSION})
set(OpenFOAM_DIR $ENV{WM_PROJECT_DIR})
set(OpenFOAM_LIB $ENV{FOAM_LIBBIN})
set(OpenFOAM_SRC $ENV{FOAM_SRC})
set(OpenFOAM_USR_BIN $ENV{FOAM_USER_APPBIN})
set(OpenFOAM_USR_LIB $ENV{FOAM_USER_LIBBIN})

add_compile_options(
        -std=c++11
        -m64
        -pthread
        -DOPENFOAM=2212
        -DWM_DP
        -DWM_LABEL_SIZE=32
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof
        -Wno-attributes
        -Wno-unknown-pragmas
        -O3
        -DNoRepository
        -ftemplate-depth-100)

message(STATUS "solver dir : ${OpenFOAM_USR_BIN}")
message(STATUS "lib dir : ${OpenFOAM_USR_LIB}")

add_subdirectory(MultiFoam)
add_subdirectory(newBC)
add_subdirectory(test)