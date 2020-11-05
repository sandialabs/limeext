#include "THINWALL_ModelEval.hpp"

#include "LIME_Problem_Manager.hpp"

using LIME::Model_Evaluator;
using LIME::Problem_Manager;

extern "C" {
  void setup_thinwall_(); 
  void solve_thinwall_(int*); 
  void finish_thinwall_(); 
}

THINWALL_ModelEval::THINWALL_ModelEval(const Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm)
{
  setup_thinwall_();
}

//------------THINWALL-----------------------------------------------------------------

THINWALL_ModelEval::~THINWALL_ModelEval()
{
   finish_thinwall_();
}


//----THINWALL-------------------------------------------------------------------------

void THINWALL_ModelEval::solve_standalone() 
{
    int nstep=25;
    solve_thinwall_(&nstep);
}
