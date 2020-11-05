// ********************************************************************************************
// NEUTRON_ModelEval0.cpp
// This is a LIME model evaluator for a standalone wrap of ss_neutron0.f90 code
// ********************************************************************************************
//
#include "NEUTRON_ModelEval0.hpp"
#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Problem_Manager;

// Data from fortran routines are referenced using an extern C block as 
// shown below. The "name mangling" conventions are taken care of automatically
// in LIME by #including the file LIME_fortran_mangling.h

extern "C" {
  #define setup_neutron0  LIME_MANGLING_GLOBAL(setup_neutron0, SETUP_NEUTRON0)
  #define solve_neutron0  LIME_MANGLING_GLOBAL(solve_neutron0, SOLVE_NEUTRON0)
  #define finish_neutron0 LIME_MANGLING_GLOBAL(finish_neutron0, FINISH_NEUTRON0)

  void setup_neutron0();
  void solve_neutron0(); 
  void finish_neutron0(); 
}

//------------NEUTRON_ModelEval0 C++ constructor --------------------------------------------

NEUTRON_ModelEval0::NEUTRON_ModelEval0(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name)
{
  setup_neutron0();
}

//-------------NEUTRON_ModelEval0 solve_standalone ------------------------------------------

void NEUTRON_ModelEval0::solve_standalone() 
{
    solve_neutron0();
}

//------------NEUTRON_ModelEval0 C++ destructor ---------------------------------------------

NEUTRON_ModelEval0::~NEUTRON_ModelEval0()
{
   finish_neutron0();
}

