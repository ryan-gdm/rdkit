rdkit_library(Mordred Mordred.cpp LINK_LIBRARIES DataStructs Descriptors Fingerprints GraphMol PartialCharges RDGeneral SmilesParse Subgraphs SubstructMatch)
target_compile_definitions(Mordred PRIVATE RDKIT_OSMORDRED_BUILD)

rdkit_headers(Mordred.h DEST Mordred)

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
    PATH_SUFFIXES /eigen
    REQUIRED
)
include_directories(Mordred ${EIGEN3_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
include_directories(Mordred ${LAPACK_INCLUDE_DIRS})
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # ${LAPACK_INCLUDE_DIRS} isn't working yet on linux so hack it:
    include_directories(Mordred /usr/include/lapacke)
endif ()

link_libraries(Mordred ${LAPACK_LIBRARIES})

# rdkit_test(testMordred testMordred.cpp LINK_LIBRARIES Mordred )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()