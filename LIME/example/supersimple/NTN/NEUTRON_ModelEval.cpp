#include "NEUTRON_ModelEval.hpp"

#include "LIME_Problem_Manager.hpp"

 extern "C" {
   void setup_neutron_(); 
   void solve_neutron_(); 
   void finish_neutron_(); 

 }

NEUTRON_ModelEval::NEUTRON_ModelEval(const LIME::Problem_Manager & pm, const string & name) :
  problem_manager_api(pm)
{
  setup_neutron_();
}

//------------NEUTRON-----------------------------------------------------------------

NEUTRON_ModelEval::~NEUTRON_ModelEval()
{
   finish_neutron_();
}


//----NEUTRON-------------------------------------------------------------------------

void NEUTRON_ModelEval::solve_standalone() 
{
    solve_neutron_();
}
