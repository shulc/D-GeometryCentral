cmake_minimum_required(VERSION 3.15)
project(gc_c CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# geometry-central build options — disable everything we don't need
set(GC_BUILD_TESTS        OFF CACHE BOOL "" FORCE)
set(GC_BUILD_EXAMPLES     OFF CACHE BOOL "" FORCE)
set(GC_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS     OFF CACHE BOOL "" FORCE)

# geometry-central uses FetchContent for Eigen and Spectra by default.
# Let it do its thing; we only need the CPU-side library.
add_subdirectory(extern/geometry-central EXCLUDE_FROM_ALL)

# ---- The C shim ------------------------------------------------------------
add_library(gc_c STATIC csrc/gc_c.cpp)
target_include_directories(gc_c PUBLIC include)
# geometry-central target transitively pulls Eigen headers via its own
# include dirs, so linking against it is sufficient.
target_link_libraries(gc_c PUBLIC geometry-central)
set_target_properties(gc_c PROPERTIES POSITION_INDEPENDENT_CODE ON)
