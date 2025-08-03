
add_library(keroutils
    encoding.cpp
    sequences.cpp
    skmers.cpp
    colinear.cpp
)

target_link_libraries(keroutils
    kero
)

target_include_directories(keroutils
    PUBLIC
    .
)
