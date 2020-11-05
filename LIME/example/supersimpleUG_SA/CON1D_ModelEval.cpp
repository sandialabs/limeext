// *************************************************************************************** 
// CON1D_ModelEval.cpp 
// This is a LIME model evaluator for a standalone wrap of the ss_con1d.f90 code 
// *************************************************************************************** 
//
#include "CON1D_ModelEval.hpp"
#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Problem_Manager;

// Data from fortran routines are referenced using an extern C block as 
// shown below. The "name mangling" conventions are taken care of automatically
// in LIME by #including the file LIME_fortran_mangling.h

 extern "C" { 
   #define setup_con1d  LIME_MANGLING_GLOBAL(setup_con1d, SETUP_CON1D)
   #define solve_con1d  LIME_MANGLING_GLOBAL(solve_con1d, SOLVE_CON1D)
   #define finish_con1d LIME_MANGLING_GLOBAL(finish_con1d, FINISH_CON1D)

   void setup_con1d(); 
   void solve_con1d(int*); 
   void finish_con1d(); 
 }

//------------CON1D_ModelEval C++ constructor ------------------------------------------

CON1D_ModelEval::CON1D_ModelEval(LIME::Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name)
{
  setup_con1d();           // This calls the con1d code's “setup_con1d” routine
}

//------------CON1D_ModelEval solve_standalone -----------------------------------------

void CON1D_ModelEval::solve_standalone() 
{
    int nstep=25;
    solve_con1d(&nstep);   // This calls the con1d code's “solve_con1d” routine
}

//------------CON1D_ModelEval C++ destructor -------------------------------------------

CON1D_ModelEval::~CON1D_ModelEval()
{
   finish_con1d();         // This calls the con1d code's “finish_con1d” routine
}
