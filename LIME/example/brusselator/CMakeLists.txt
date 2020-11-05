
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
IF (LIME_ENABLE_TriKota)
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()

SET(example_brusselator_SOURCES
    Brusselator_Equation.cpp
    Brusselator_Species.cpp
    Brusselator_Temperature.cpp
    Basis.cpp
   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    brusselator
    SOURCES ${example_brusselator_SOURCES}
    Example_Brusselator.cpp
    NAME brusselator_fixed_point
    ARGS "--setupxml=Problem_Manager_setup_for_fp.xml --outputdir=output_fp"
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )

TRIBITS_ADD_TEST(
    brusselator
    NAME brusselator_jfnk
    ARGS "--setupxml=Problem_Manager_setup_for_jfnk.xml --outputdir=output_jfnk"
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )


TRIBITS_ADD_EXECUTABLE(
    brusselator_with_p_and_g
    SOURCES ${example_brusselator_SOURCES}
      Example_Brusselator_p_g.cpp
    )
TRIBITS_ADD_TEST(
    brusselator_with_p_and_g
    NAME brusselator_p_and_g
    ARGS "--setupxml=Problem_Manager_setup_for_jfnk.xml --outputdir=output_pg --paramsxml=ParameterSelection.xml"
    NUM_MPI_PROCS 1
    CATEGORIES NIGHTLY
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )
# Above test takes more than 3 minutes to run on pu241 and therefore
# should not be a *BASIC* test.  It will still run in CI and Nightly
# testing.


IF (LIME_ENABLE_TriKota)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        brusselator_trikota
	SOURCES ${example_brusselator_SOURCES}
	Example_Brusselator_TriKota.cpp
	NAME brusselator_trikota
	ARGS "--setupxml=Problem_Manager_setup_for_jfnk.xml --outputdir=output_trikota --paramsxml=ParameterSelection.xml"
	COMM serial mpi
	NUM_MPI_PROCS 1
	PASS_REGULAR_EXPRESSION "Test Passed!"
	)
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_brusselator_files
    SOURCE_FILES Problem_Manager_setup_fp.xml Problem_Manager_setup_for_fp.xml Problem_Manager_setup_for_jfnk.xml dakota.in
    Brusselator_Species_nox.xml
    Brusselator_Temperature_nox.xml
    ParameterSelection.xml
    EXEDEPS brusselator
    )

