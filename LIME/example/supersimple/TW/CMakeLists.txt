INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
    thinwall_sad 
    SOURCES thinwall_sad.cpp THINWALL_ModelEval.cpp ss_thinwall.f90
    DEPLIBS lime
    )

