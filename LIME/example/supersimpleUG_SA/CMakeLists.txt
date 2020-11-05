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
    con1d_sad
    SOURCES con1d_sad.cpp CON1D_ModelEval.cpp
    DEPLIBS lime supersimpleUG_fortran_lib
    NAME con1d_sad 
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     neutron0_sad
     SOURCES neutron0_sad.cpp NEUTRON_ModelEval0.cpp
     DEPLIBS lime supersimpleUG_fortran_lib
     NAME neutron0_sad
     COMM serial mpi
     NUM_MPI_PROCS 1
     )

 TRIBITS_ADD_EXECUTABLE_AND_TEST(
     thinwall_sad
     SOURCES thinwall_sad.cpp THINWALL_ModelEval.cpp
     DEPLIBS lime supersimpleUG_fortran_lib
     NAME thinwall_sad
     COMM serial mpi
     NUM_MPI_PROCS 1
     )

