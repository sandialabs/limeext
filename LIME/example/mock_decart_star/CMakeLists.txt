
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(example_mock_decart_star_SOURCES
    main.cpp
    MockModelEvals.cpp
   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    destar
    SOURCES ${example_mock_decart_star_SOURCES}
    main.cpp
    NAME mock_decart_star
    ARGS "--setupxml=Problem_Manager_setup_for_fp.xml --outputdir=output_destar"
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test Passed!"
    )


# Above test takes more than 3 minutes to run on pu241 and therefore
# should not be a *BASIC* test.  It will still run in CI and Nightly
# testing.


TRIBITS_COPY_FILES_TO_BINARY_DIR(example_decart_star_files
    SOURCE_FILES Problem_Manager_setup_fp.xml Problem_Manager_setup_for_fp.xml
    EXEDEPS destar
    )

