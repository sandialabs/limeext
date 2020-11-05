#include "CON1D_ModelEval.hpp"

#include "LIME_Problem_Manager.hpp"

 extern "C" {
   void setup_con1d_(); 
   void solve_con1d_(int*); 
   void finish_con1d_(); 

 }

//------------CON1D setup ---------------------------------------------------------

CON1D_ModelEval::CON1D_ModelEval(const LIME::Problem_Manager & pm, const string & name) :
  problem_manager_api(pm)
{
  setup_con1d_();
}


//------------CON1D solve_standalone -----------------------------------------------

void CON1D_ModelEval::solve_standalone() 
{
    int nstep=25;
    solve_con1d_(&nstep);
}


//------------CON1D finish ---------------------------------------------------------

CON1D_ModelEval::~CON1D_ModelEval()
{
   finish_con1d_();
}

