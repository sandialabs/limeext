
ADD_SUBDIRECTORY(fixtures)

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
IF ( LIME_ENABLE_TriKota )
  LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../../packages/TriKota/Dakota/install/lib)
ENDIF()

SET(HEADERS "")
SET(SOURCES "")

APPEND_SET(SOURCES
  UnitTestMain.cpp
  UnitTestAdjustTimeStep.cpp
  UnitTestAndersonAcceleration.cpp
  UnitTestFixedPoint.cpp
  UnitTestHierarchySolves.cpp
  UnitTestHierarchySupport.cpp
  UnitTestNonlinElim.cpp
  UnitTestProblemManagerOperators.cpp
  UnitTestStandardTransfers.cpp
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lime_unit_tests
    SOURCES ${SOURCES}
    DEPLIBS lime_test_fixtures
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

