cmake_minimum_required(VERSION 3.19)

project(kero-tools VERSION 0.1)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type: Debug or Release" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

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

    if(APPLE)
        add_compile_options(-flto)
    else()
        add_compile_options(-flto=auto)
    endif()
endif()

add_subdirectory(extern/kero-api/)
add_subdirectory(src/)
add_subdirectory(apps/)