
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

SET(LIME_DEBUG ${HAVE_LIME_DEBUG})

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

FortranCInterface_HEADER(LIME_fortran_mangling.h
                         MACRO_NAMESPACE "LIME_MANGLING_")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_fortran_mangling.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# FAD
APPEND_SET(HEADERS
  LIME_Broadcast_File.hpp
  LIME_Problem_As_Operator.hpp
  LIME_Model_Evaluator.hpp
  LIME_Strings.hpp
  LIME_Data_Transfer_Operator.hpp
  LIME_Elimination_Module.hpp
  LIME_FixedPoint_Accelerators.hpp
  LIME_Hierarchy_Problem_Manager.hpp
  LIME_Problem_Manager.hpp
  LIME_PrePostOperator.hpp
  ) 

APPEND_SET(SOURCES
  LIME_Broadcast_File.cpp
  LIME_Strings.cpp
  LIME_Hierarchy_Problem_Manager.cpp
  LIME_FixedPoint_Accelerators.cpp
  LIME_PrePostOperator.cpp
  LIME_Problem_Manager.cpp
  LIME_Model_Evaluator.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  lime
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

