
INCLUDE(TribitsLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

APPEND_SET(HEADERS
  FiveLinearEQs.hpp
  ThreeNonlinEQsRes.hpp
  )

APPEND_SET(SOURCES
  FiveLinearEQs.cpp
  ThreeNonlinEQsRes.cpp
  )

TRIBITS_ADD_LIBRARY(
  lime_test_fixtures
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
