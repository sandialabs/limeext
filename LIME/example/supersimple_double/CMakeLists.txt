INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF ( LIME_ENABLE_TriKota )
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lime_mpd_double
    SOURCES lime_mpd_double.cpp ss_con1d.f90 CON1D_ModelEval0.cpp CON1D_ModelEval0_w_Resid.cpp  ss_thinwall.f90 THINWALL_ModelEval0.cpp THINWALL_ModelEval0_w_Resid.cpp
    DEPLIBS lime
    NAME lime_mpd_double
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )
