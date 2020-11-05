#include <algorithm>
#include <iostream>

#include "NEUTRON_ModelEval.hpp"
#include "CON1D_ModelEval.hpp"

#include "LIME_Data_Transfer_Operator.hpp"
#include "LIME_fortran_mangling.h"
#include "LIME_Problem_Manager.hpp"

using LIME::Data_Transfer_Operator;
using LIME::Model_Evaluator;
using LIME::Problem_Manager;

extern "C" {

  // fortran module data

  #define neutron_q LIME_MANGLING_MODULE(neutron_mod, q, NEUTRON_MOD, Q)
  extern float neutron_q[4];

  #define neutron_t LIME_MANGLING_MODULE(neutron_mod, t, NEUTRON_MOD, T)
  extern float neutron_t[4];

  #define neutron_r_max LIME_MANGLING_MODULE(neutron_mod, r_max, NEUTRON_MOD, R_MAX)
  extern float neutron_r_max;

  // fortran subroutines
  #define neutron_setup LIME_MANGLING_GLOBAL(setup_neutron, SETUP_NEUTRON)
  void neutron_setup();

  #define solve_neutron LIME_MANGLING_GLOBAL(solve_neutron, SOLVE_NEUTRON)
  void solve_neutron();

  #define residual_neutron LIME_MANGLING_GLOBAL(residual_neutron, RESIDUAL_NEUTRON)
  void residual_neutron();

  #define finish_neutron LIME_MANGLING_GLOBAL(finish_neutron, FINISH_NEUTRON)
  void finish_neutron();
}

// Set the "correct" gold results
const float NEUTRON_ModelEval::gold_qs[] = {  3558.78, 3807.15, 4412.58 };

//------------NEUTRON setup -----------------------------------------------------------------

NEUTRON_ModelEval::NEUTRON_ModelEval(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name), q(neutron_q), t(neutron_t)
{
  neutron_setup();
}

//-------------NEUTRON solve -----------------------------------------------------------------

void NEUTRON_ModelEval::solve_standalone() 
{
  //std::cout << "NEUTRON_modelEval::solve_standalone() called" << std::endl;
    solve_neutron();
}


//-------------NEUTRON converged query --------------------------------------------------------

bool NEUTRON_ModelEval::is_converged() 
{
  residual_neutron();
  //std::cout << "NEUTRON_ModelEval R_max = " << VAR_R_MAX << std::endl;
    return( neutron_r_max < 1.e-5 );
}

//------------NEUTRON finish -----------------------------------------------------------------

NEUTRON_ModelEval::~NEUTRON_ModelEval()
{
   finish_neutron();
}


//--------------------------------------------------------------------------------------------
conduction_2_neutronics::conduction_2_neutronics(Teuchos::RCP<Model_Evaluator> from, 
                                                 Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{ }

//--------------------------------------------------------------------------------------------
bool conduction_2_neutronics::perform_data_transfer() const
{
   const Teuchos::RCP<CON1D_ModelEval> & con1d_me =  Teuchos::rcp_dynamic_cast<CON1D_ModelEval>(source());
   const Teuchos::RCP<NEUTRON_ModelEval> & neutron_me =  Teuchos::rcp_dynamic_cast<NEUTRON_ModelEval>(target());

  std::copy(con1d_me->t, con1d_me->t+3, neutron_me->t);
  //std::cout << "conduction_2_neutronics transfer" << std::endl;
  //std::cout << neutron_me->t[1] << std::endl;
  
  return false;
}

//--------------------------------------------------------------------------------------------
neutronics_2_conduction::neutronics_2_conduction(Teuchos::RCP<Model_Evaluator> from, Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{
}

//--------------------------------------------------------------------------------------------
bool neutronics_2_conduction::perform_data_transfer() const
{
   const Teuchos::RCP<NEUTRON_ModelEval> & neutron_me = Teuchos::rcp_dynamic_cast<NEUTRON_ModelEval>(source());
   const Teuchos::RCP<CON1D_ModelEval> & con1d_me = Teuchos::rcp_dynamic_cast<CON1D_ModelEval>(target());

   std::copy(neutron_me->q, neutron_me->q+3, con1d_me->qs);
   //std::cout << "neutronics_2_conduction transfer" << std::endl;
  //std::cout << con1d_me->qs[1] << std::endl;
  //std::cout << neutron_me->q[0] << ", " << neutron_me->q[1] << ", " << neutron_me->q[2] << ", " << con1d_me->qs[1] << std::endl;

  return false;
}
