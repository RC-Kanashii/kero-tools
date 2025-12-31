cmake_minimum_required(VERSION 3.19)

project(kero-tools VERSION 0.1)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type: Debug or Release" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

# ===== KERO Ablation Study Configuration =====
# Set KERO storage mode for ablation experiments
# Options: ROW, COLUMNAR_NOCOMP, COLUMNAR_COMP
set(KERO_STORAGE_MODE "COLUMNAR_COMP" CACHE STRING
        "KERO storage mode for ablation study: ROW, COLUMNAR_NOCOMP, or COLUMNAR_COMP")
set_property(CACHE KERO_STORAGE_MODE PROPERTY STRINGS
        "ROW" "COLUMNAR_NOCOMP" "COLUMNAR_COMP")

# Define corresponding macro based on storage mode
if(KERO_STORAGE_MODE STREQUAL "ROW")
    add_compile_definitions(KERO_MODE_ROW)
    message(STATUS "KERO Mode: Row-oriented storage (no compression)")
elseif(KERO_STORAGE_MODE STREQUAL "COLUMNAR_NOCOMP")
    add_compile_definitions(KERO_MODE_COLUMNAR_NOCOMP)
    message(STATUS "KERO Mode: Columnar storage (no compression)")
elseif(KERO_STORAGE_MODE STREQUAL "COLUMNAR_COMP")
    add_compile_definitions(KERO_MODE_COLUMNAR_COMP)
    message(STATUS "KERO Mode: Columnar storage + SIMD compression (default)")
else()
    message(FATAL_ERROR "Invalid KERO_STORAGE_MODE: ${KERO_STORAGE_MODE}. Must be ROW, COLUMNAR_NOCOMP, or COLUMNAR_COMP")
endif()
# ===== End KERO Ablation Study Configuration =====

if(APPLE)
    set(CMAKE_XCODE_GENERATE_SCHEME ON)
    set(CMAKE_XCODE_SCHEME_ADDRESS_SANITIZER ON)
    set(CMAKE_XCODE_SCHEME_ADDRESS_SANITIZER_USE_AFTER_RETURN ON)
    add_compile_options(-Wall -Wextra)

elseif(UNIX AND NOT APPLE)
    add_compile_options(-Wall -Wextra -pedantic)

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-fvar-tracking -fvar-tracking-assignments)
    endif()
endif()

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
if(build_type STREQUAL debug)
    add_compile_definitions(DEBUG_MODE)
    add_compile_options(-g3 -O0 -fno-inline -fno-omit-frame-pointer)

    if(APPLE)
        add_compile_options(-fstandalone-debug)
    else()
        add_compile_options(-fno-eliminate-unused-debug-symbols)
    endif()

elseif(build_type STREQUAL release)
    add_compile_options(-O3 -DNDEBUG)
endif()

add_subdirectory(extern/kero-api/)
add_subdirectory(src/)
add_subdirectory(apps/)