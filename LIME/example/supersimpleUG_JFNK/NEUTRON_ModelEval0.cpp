// ********************************************************************************************
//
// NEUTRON_ModelEval0.cpp
// This is a ss_neutron0 LIME model evaluator used for coupling with ss_con1d
// It includes data fransfer operators needed for this coupling
// 
// ********************************************************************************************
#include <iostream>

#include "NEUTRON_ModelEval0.hpp"
#include "CON1D_ModelEval0.hpp"

#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Data_Transfer_Operator;
using LIME::Model_Evaluator;
using LIME::Problem_Manager;

// Data from fortran modules and fortran routines are referenced using an extern C 
// block as shown below. The "name mangling" conventions are taken care of automatically
// in LIME by #including the file LIME_fortran_mangling.h

extern "C" {

   // fortran90 module data that we need access to

  #define neutron0_q     LIME_MANGLING_MODULE(neutron0_mod,q, NEUTRON0_MOD,Q)
  #define neutron0_t     LIME_MANGLING_MODULE(neutron0_mod,t, NEUTRON0_MOD,T)
  #define neutron0_r_max LIME_MANGLING_MODULE(neutron0_mod,r_max, NEUTRON0_MOD,R_MAX)

  extern float neutron0_q[4];
  extern float neutron0_t[4];
  extern float neutron0_r_max;

  // fortran90 subroutines that we need to call

  #define setup_neutron0  LIME_MANGLING_GLOBAL(setup_neutron0, SETUP_NEUTRON0)
  #define solve_neutron0  LIME_MANGLING_GLOBAL(solve_neutron0, SOLVE_NEUTRON0)
  #define finish_neutron0 LIME_MANGLING_GLOBAL(finish_neutron0, FINISH_NEUTRON0)

  void setup_neutron0();
  void solve_neutron0(); 
  void finish_neutron0(); 
}

// Set the "correct" gold results
// const float NEUTRON_ModelEval0::gold_qs[] = {  3558.78, 3807.15, 4412.58 };

//------------NEUTRON_ModelEval0 C++ constructor --------------------------------------------

NEUTRON_ModelEval0::NEUTRON_ModelEval0(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name), q(neutron0_q), t(neutron0_t)
{
  setup_neutron0();
}

//-------------NEUTRON_ModelEval0 solve_standalone ------------------------------------------

void NEUTRON_ModelEval0::solve_standalone() 
{
  //std::cout << "NEUTRON_ModelEval0::solve_standalone() called" << std::endl;
    solve_neutron0();
}

//------------NEUTRON_ModelEval0 C++ destructor ---------------------------------------------

NEUTRON_ModelEval0::~NEUTRON_ModelEval0()
{
   finish_neutron0();
}

//------------ NEUTRON_ModelEval0 Implimentation of perform_elimination ------------------------

bool NEUTRON_ModelEval0::perform_elimination()
{
   solve_neutron0();
   return true;
}


//-------------------------------------------------------------------------------------------
// The required data transfer operators are included next
//-------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------
conduction_2_neutronics::conduction_2_neutronics(Teuchos::RCP<Model_Evaluator> from, 
                                                 Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{
}

//-------------------------------------------------------------------------------------------
bool conduction_2_neutronics::perform_data_transfer() const
{
   const Teuchos::RCP<CON1D_ModelEval0> & con1d_me = Teuchos::rcp_dynamic_cast<CON1D_ModelEval0>(source());
   const Teuchos::RCP<NEUTRON_ModelEval0> & neutron_me = Teuchos::rcp_dynamic_cast<NEUTRON_ModelEval0>(target());

  std::copy(con1d_me->t, con1d_me->t+3, neutron_me->t);
  //std::cout << "conduction_2_neutronics transfer" << std::endl;
  //std::cout << neutron_me->t[1] << std::endl;
  
  return true;
}

//-------------------------------------------------------------------------------------------
neutronics_2_conduction::neutronics_2_conduction(Teuchos::RCP<Model_Evaluator> from, 
                                                 Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{ }

//--------------------------------------------------------------------------------------------
bool neutronics_2_conduction::perform_data_transfer() const
{
   const Teuchos::RCP<NEUTRON_ModelEval0> & neutron_me = Teuchos::rcp_dynamic_cast<NEUTRON_ModelEval0>(source());
   const Teuchos::RCP<CON1D_ModelEval0> & con1d_me = Teuchos::rcp_dynamic_cast<CON1D_ModelEval0>(target());

   std::copy(neutron_me->q, neutron_me->q+3, con1d_me->qs);
   //std::cout << "neutronics_2_conduction transfer" << std::endl;
   //std::cout << con1d_me->qs[1] << std::endl;
   //std::cout << neutron_me->q[0] << ", " << neutron_me->q[1] << ", " << neutron_me->q[2] << ", " << con1d_me->qs[1] << std::endl;

  return true;
}


//--------------------------------------------------------------------------------------------

