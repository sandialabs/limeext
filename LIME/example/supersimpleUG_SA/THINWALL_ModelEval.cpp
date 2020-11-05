// *************************************************************************************** 
// THINWALL_ModelEval.cpp 
// This is a LIME model evaluator for a standalone wrap of the super-simple thinwall code 
// *************************************************************************************** 
//
#include "THINWALL_ModelEval.hpp"
#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Problem_Manager;

// Data from fortran routines are referenced using an extern C block as 
// shown below. The "name mangling" conventions are taken care of automatically
// in LIME by #including the file LIME_fortran_mangling.h

 extern "C" {
   #define setup_thinwall     LIME_MANGLING_GLOBAL(setup_thinwall, SETUP_THINWALL)
   #define finish_thinwall    LIME_MANGLING_GLOBAL(finish_thinwall, FINISH_THINWALL)
   #define solve_thinwall     LIME_MANGLING_GLOBAL(solve_thinwall, SOLVE_THINWALL)

   void setup_thinwall(); 
   void finish_thinwall(); 
   void solve_thinwall(int*); 
 }

//------------THINWALL ModelEval C++ constructor  ---------------------------------------------

THINWALL_ModelEval::THINWALL_ModelEval(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name)
{
  setup_thinwall();
}

//------------THINWALL ModelEval C++ destructor ----------------------------------------------

THINWALL_ModelEval::~THINWALL_ModelEval()
{
   finish_thinwall();
}

//------------THINWALL_ModelEval Implimentation of solve_standalone --------------------------

void THINWALL_ModelEval::solve_standalone() 
{
   int nstep=25;
   solve_thinwall_(&nstep);
}

