INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF ( LIME_ENABLE_TriKota )
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()
  
IF ( LIME_ENABLE_NOX )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lime_1d_conj_heat 
    SOURCES lime_mpd.cpp ConvDiff_EQ.cpp
    DEPLIBS lime
    NAME 1d_conj_heat 
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_1d_conj_heat_files 
        SOURCE_FILES Region_1_nox.xml
                     Region_2_nox.xml
    )

