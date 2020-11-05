INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(LIME)

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

COMBINED_OPTION( ${PACKAGE_NAME}_ENABLE_JFNK
  DEP_OPTIONS_NAMES LIME_ENABLE_NOX LIME_ENABLE_Ifpack LIME_ENABLE_ML LIME_ENABLE_Amesos LIME_ENABLE_AztecOO
  DOCSTR "Enables Jacobian-Free Newton-Krylov solves via NOX and associated linear solver/preconditioners"
  )

IF(LIME_ENABLE_JFNK)
  SET(HAVE_LIME_JFNK ON)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
