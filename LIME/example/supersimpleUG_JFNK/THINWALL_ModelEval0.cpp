// ***************************************************************************************** 
// 
// THINWALL_ModelEval0.cpp 
// This is the ss_thinwall LIME model evaluator used for a simple coupling with ss_con1d
// 
// *****************************************************************************************
#include "THINWALL_ModelEval0.hpp"
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

   // fortran module data

  #define thinwall_t    LIME_MANGLING_MODULE(thinwall_mod, t   , THINWALL_MOD, T   )
  #define thinwall_trad LIME_MANGLING_MODULE(thinwall_mod, trad, THINWALL_MOD, TRAD)
  #define thinwall_dt   LIME_MANGLING_MODULE(thinwall_mod, dt  , THINWALL_MOD, DT  )
  #define thinwall_time LIME_MANGLING_MODULE(thinwall_mod, time, THINWALL_MOD, TIME)
  #define thinwall_r    LIME_MANGLING_MODULE(thinwall_mod, r   , THINWALL_MOD, R   )
  #define thinwall_tol  LIME_MANGLING_MODULE(thinwall_mod, tol , THINWALL_MOD, TOL )

  extern float thinwall_t[1];   //     thinwall temperature 
  extern float thinwall_trad;   //     thinwall radiation BC temperature 
  extern float thinwall_dt;     //     solution time step
  extern float thinwall_time;   //     current simulation time
  extern float thinwall_r[1];   //     residual value
  extern float thinwall_tol ;   //     convergence tolerance

  #define VAR_T     thinwall_t
  #define VAR_TRAD  thinwall_trad
  #define VAR_DT    thinwall_dt
  #define VAR_TMAX  100000.0               //  hardwired to large value
  #define VAR_TIME  thinwall_time
  #define VAR_R_MAX thinwall_r
  #define VAR_TOL   thinwall_tol

  // fortran subroutines

  #define setup_thinwall     LIME_MANGLING_GLOBAL(setup_thinwall, SETUP_THINWALL)
  #define finish_thinwall    LIME_MANGLING_GLOBAL(finish_thinwall, FINISH_THINWALL)
  #define time_step_thinwall LIME_MANGLING_GLOBAL(time_step_thinwall, TIME_STEP_THINWALL)
  #define update_thinwall    LIME_MANGLING_GLOBAL(update_thinwall, UPDATE_THINWALL)
  #define residual_thinwall  LIME_MANGLING_GLOBAL(residual_thinwall, RESIDUAL_THINWALL)

   void setup_thinwall(); 
   void finish_thinwall(); 
   void time_step_thinwall(); 
   void update_thinwall(); 
   void residual_thinwall(); 
 }

// Set the "correct" gold results
// const float THINWALL_ModelEval0::gold_w_temp = 346.6981;

//------------THINWALL ModelEval0 C++ constructor  ---------------------------------------------

THINWALL_ModelEval0::THINWALL_ModelEval0(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name), t(VAR_T), trad(VAR_TRAD)
{
  setup_thinwall();
}

//------------THINWALL ModelEval0 C++ destructor ----------------------------------------------

THINWALL_ModelEval0::~THINWALL_ModelEval0()
{
   finish_thinwall();
}


//------------THINWALL_ModelEval0 Implimentation of solve_standalone --------------------------

void THINWALL_ModelEval0::solve_standalone() 
{
   //std::cout << "THINWALL_ModelEval0::solve_standalone called" << std::endl;
   time_step_thinwall();
}

//------------THINWALL_ModelEval0 Implimentation of is_converged -----------------------------

bool THINWALL_ModelEval0::is_converged() 
{
  residual_thinwall();
  //std::cout << "DEBUG THINWALL_ModelEval0, VAR_R_MAX=" << VAR_R_MAX[0] << std::endl;
    return( VAR_R_MAX[0] < VAR_TOL );
}


//------------THINWALL_ModelEval0 Implimentation of get_time_step ----------------------------

double THINWALL_ModelEval0::get_time_step() const
{
   return VAR_DT;
}


//------------THINWALL_ModelEval0 Implimentation of get_max_time -----------------------------

double THINWALL_ModelEval0::get_max_time() const
{
   return VAR_TMAX;
}

//------------THINWALL_ModelEval0 Implimentation of get_current_time -------------------------

double THINWALL_ModelEval0::get_current_time() const
{
   return VAR_TIME;
}

//------------THINWALL_ModelEval0 Implimentation of update_time ------------------------------

void THINWALL_ModelEval0::update_time() 
{
   update_thinwall();
}


//-------------------------------------------------------------------------------------------
// The required data transfer operators are included next
//-------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------
conduction_2_thinwall::conduction_2_thinwall(Teuchos::RCP<Model_Evaluator> from, 
                                             Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{
}

//--------------------------------------------------------------------------------------------
bool conduction_2_thinwall::perform_data_transfer() const
{
   const Teuchos::RCP<CON1D_ModelEval0> & con1d_me = Teuchos::rcp_dynamic_cast<CON1D_ModelEval0>(source());
   const Teuchos::RCP<THINWALL_ModelEval0> & thinwall_me = Teuchos::rcp_dynamic_cast<THINWALL_ModelEval0>(target());

  thinwall_me->trad = con1d_me->t[2];
  //std::cout << "conduction_2_thinwall transfer" << std::endl;
  //std::cout << thinwall_me->trad << std::endl;
  
  return true;
}

//--------------------------------------------------------------------------------------------
thinwall_2_conduction::thinwall_2_conduction(Teuchos::RCP<Model_Evaluator> from, 
                                             Teuchos::RCP<Model_Evaluator> to)
      : Data_Transfer_Operator(from, to)
{ }

//--------------------------------------------------------------------------------------------
bool thinwall_2_conduction::perform_data_transfer() const
{
   const Teuchos::RCP<THINWALL_ModelEval0> & thinwall_me = Teuchos::rcp_dynamic_cast<THINWALL_ModelEval0>(source());
   const Teuchos::RCP<CON1D_ModelEval0> & con1d_me = Teuchos::rcp_dynamic_cast<CON1D_ModelEval0>(target());

   con1d_me->trbc = thinwall_me->t[0];
   //std::cout << "thinwall_2_conduction transfer" << std::endl;
   //std::cout << con1d_me->trbc << std::endl;

  return true;
}

