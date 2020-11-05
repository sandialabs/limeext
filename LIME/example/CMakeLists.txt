
IF (LIME_ENABLE_JFNK)
  ADD_SUBDIRECTORY(brusselator)
  ADD_SUBDIRECTORY(supersimple_double)
  ADD_SUBDIRECTORY(supersimple_float)
ENDIF()
ADD_SUBDIRECTORY(supersimple)
#ADD_SUBDIRECTORY(supersimpleUG_SA)
#ADD_SUBDIRECTORY(supersimpleUG_FP)
#ADD_SUBDIRECTORY(supersimpleUG_JFNK)
ADD_SUBDIRECTORY(1D_conj_heat)
#ADD_SUBDIRECTORY(mock_decart_star)

