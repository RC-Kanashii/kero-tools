add_library(kerotl
    bucket.cpp
    compact.cpp
    datarm.cpp
    disjoin.cpp
    index.cpp
    instr.cpp
    merge.cpp
    outstr.cpp
    shuffle.cpp
    sort.cpp
    split.cpp
    translate.cpp
    validate.cpp
    query.cpp
)

target_include_directories(kerotl
    PUBLIC
    ./include/
    "${PROJECT_SOURCE_DIR}/extern/CLI11/include/"
    "${PROJECT_SOURCE_DIR}/extern/kero-api/include/"
    "${PROJECT_SOURCE_DIR}/extern/logfault"
)

# link libraries
find_package(OpenMP)
target_link_libraries(kerotl
    PUBLIC
    keroutils
    kero
    OpenMP::OpenMP_CXX
)
