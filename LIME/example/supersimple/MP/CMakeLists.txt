INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
TRIBITS_ADD_EXECUTABLE(
    lime_mpd 
    SOURCES lime_mpd.cpp ss_con1d.f90 NEUTRON_ModelEval.cpp ss_neutron.f90 CON1D_ModelEval.cpp ss_con1d.f90 
    DEPLIBS lime
    )

TRIBITS_ADD_EXECUTABLE(
    con1d_sad 
    SOURCES con1d_sad.cpp CON1D_ModelEval.cpp ss_con1d.f90
    DEPLIBS lime
    )

TRIBITS_ADD_EXECUTABLE(
    thinwall_sad 
    SOURCES thinwall_sad.cpp THINWALL_ModelEval.cpp ss_thinwall.f90
    DEPLIBS lime
    )

TRIBITS_ADD_EXECUTABLE(
    neutron_sad
    SOURCES neutron_sad.cpp NEUTRON_ModelEval.cpp ss_neutron.f90
    DEPLIBS lime
    )

