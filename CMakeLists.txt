cmake_minimum_required(VERSION 2.6)

project(fgls)

option(FGLS_TIMING "Output timings of several parts of the program." OFF)
set(FGLS_BLAS "mkl" CACHE STRING "One of 'mkl' or 'goto'")
set(FGLS_GPU_BLAS "" CACHE STRING "Either empty or 'cublas'. TODO: add 'amdBlas' support.")
if(CMAKE_C_COMPILER MATCHES ".*vtcc")
    set(FGLS_VAMPIR_UNDERLYING_COMPILER "gcc" CACHE STRING "The underlying compiler vampirtrace should use")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(SRC
    src/common.c
    src/driver.c
    src/fgls_chol.c
    src/fgls_eigen.c
    src/io.c
    src/timing.c
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
    add_definitions(-DFGLS_GPU_SERIAL)
endif()

if(CMAKE_C_COMPILER MATCHES ".*vtcc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -vt:cc ${FGLS_VAMPIR_UNDERLYING_COMPILER} -vt:mt -vt:inst manual -vt:verbose -vt:noopari")
    add_definitions(-DVAMPIR -DVTRACE)
endif()

if(FGLS_TIMING)
    add_definitions(-DTIMING)
endif()

if(FGLS_BLAS STREQUAL "mkl")
    add_definitions(-DMKL)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} $ENV{FLAGS_MATH_LINKER}")
elseif(FGLS_BLAS STREQUAL "goto")
    message(FATAL_ERROR "TODO: Goto BLAS unimplemented.")
else()
    message(FATAL_ERROR "Unknown BLAS type: ${FGLS_BLAS}")
endif()

if(FGLS_GPU_BLAS STREQUAL "cublas")
    list(APPEND SRC src/fgls_chol_gpu.c)
    add_definitions(-DFGLS_WITH_GPU)

    include_directories("$ENV{CUDA_ROOT}/include")
    link_directories("$ENV{CUDA_ROOT}/lib64")
    link_libraries("cublas")
endif()

add_executable(driver.x ${SRC})
target_link_libraries(driver.x rt)

