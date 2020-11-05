INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF ( LIME_ENABLE_TriKota )
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()

TRIBITS_ADD_LIBRARY(
  supersimpleUG_fortran_lib
  SOURCES ss_con1d.f90 ss_neutron0.f90 ss_thinwall.f90
  TESTONLY
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lime_mpd1_jfnk
    SOURCES lime_mpd1_jfnk.cpp ss_con1d.f90 CON1D_ModelEval0.cpp CON1D_ModelEval0_w_Resid.cpp  ss_thinwall.f90 THINWALL_ModelEval0.cpp THINWALL_ModelEval0_w_Resid.cpp
    DEPLIBS lime supersimpleUG_fortran_lib
    NAME lime_mpd1_jfnk
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )


 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd_all_jfnk
     SOURCES lime_mpd_all_jfnk.cpp ss_con1d.f90 CON1D_ModelEval0.cpp CON1D_ModelEval0_w_Resid.cpp  ss_thinwall.f90 THINWALL_ModelEval0.cpp THINWALL_ModelEval0_w_Resid.cpp NEUTRON_ModelEval0.cpp
     DEPLIBS lime supersimpleUG_fortran_lib
     NAME lime_mpd_all_jfnk
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )

