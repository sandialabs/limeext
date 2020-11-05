INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF ( LIME_ENABLE_TriKota )
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()
  
TRIBITS_ADD_LIBRARY(
  supersimple_fortran_lib
  SOURCES ss_con1d.f90 ss_neutron0.f90 ss_neutron.f90 ss_thinwall.f90
  TESTONLY
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lime_mpd0
    SOURCES lime_mpd0.cpp NEUTRON_ModelEval0.cpp CON1D_ModelEval0.cpp
    DEPLIBS lime supersimple_fortran_lib
    NAME lime_mpd0 
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )

 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd 
     SOURCES lime_mpd.cpp NEUTRON_ModelEval.cpp CON1D_ModelEval.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )
 
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd1
     SOURCES lime_mpd1.cpp CON1D_ModelEval0.cpp  THINWALL_ModelEval0.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd1
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )
 
IF (LIME_ENABLE_JFNK)
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd1B
     SOURCES lime_mpd1B.cpp CON1D_ModelEval0.cpp CON1D_ModelEval0_w_Resid.cpp  THINWALL_ModelEval0.cpp THINWALL_ModelEval0_w_Resid.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd1B
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )
ENDIF()

 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd2
     SOURCES lime_mpd2.cpp NEUTRON_ModelEval.cpp CON1D_ModelEval.cpp CON1D_ModelEval_w_Resid.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd2
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )
 
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd3
     SOURCES lime_mpd3.cpp NEUTRON_ModelEval.cpp NEUTRON_ModelEval_w_Elim.cpp CON1D_ModelEval.cpp CON1D_ModelEval_w_Resid.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd3
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )
 
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     lime_mpd_all
     SOURCES lime_mpd_all.cpp NEUTRON_ModelEval0.cpp CON1D_ModelEval0.cpp THINWALL_ModelEval0.cpp
     DEPLIBS lime supersimple_fortran_lib
     NAME lime_mpd_all
     COMM serial mpi
     NUM_MPI_PROCS 1
     PASS_REGULAR_EXPRESSION "Test Passed!"
     )

