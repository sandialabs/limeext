#include <iostream>
#include "CON1D_ModelEval.hpp"

#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Model_Evaluator;
using LIME::Problem_Manager;

extern "C" {

  // fortran module data

  #define con1d_qs LIME_MANGLING_MODULE(con1d_mod, qs, CON1D_MOD, QS)
  extern float con1d_qs[3];

  #define con1d_t LIME_MANGLING_MODULE(con1d_mod, t, CON1D_MOD, T)
  extern float con1d_t[3];

  #define con1d_tr LIME_MANGLING_MODULE(con1d_mod, tr, CON1D_MOD, TR)
  extern float con1d_tr;

  #define con1d_dt LIME_MANGLING_MODULE(con1d_mod, dt, CON1D_MOD, DT)
  extern float con1d_dt;

  #define con1d_tmax LIME_MANGLING_MODULE(con1d_mod, tmax, CON1D_MOD, TMAX)
  extern float con1d_tmax;

  #define con1d_time LIME_MANGLING_MODULE(con1d_mod, time, CON1D_MOD, TIME)
  extern float con1d_time;

  #define con1d_r_max LIME_MANGLING_MODULE(con1d_mod, r_max, CON1D_MOD, R_MAX)
  extern float con1d_r_max;

  #define con1d_tol LIME_MANGLING_MODULE(con1d_mod, tol, CON1D_MOD, TOL)
  extern float con1d_tol;

  // fortran subroutines

  #define setup_con1d LIME_MANGLING_GLOBAL(setup_con1d, SETUP_CON1D)
  void setup_con1d();

  #define finish_con1d LIME_MANGLING_GLOBAL(finish_con1d, FINISH_CON1D)
  void finish_con1d();

  #define take_time_step_con1d LIME_MANGLING_GLOBAL(take_time_step_con1d, TAKE_TIME_STEP_CON1D)
  void take_time_step_con1d();

  #define update_con1d LIME_MANGLING_GLOBAL(update_con1d, UPDATE_CON1D)
  void update_con1d();

  #define con1d_residual LIME_MANGLING_GLOBAL(residual_con1d, RESIDUAL_CON1D)
  void con1d_residual();
}

// Set the "correct" gold results
const float CON1D_ModelEval::gold_temps[] = { 523.4678, 491.3500, 383.0898 };

//------------CON1D setup ---------------------------------------------------------

CON1D_ModelEval::CON1D_ModelEval(Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name), qs(con1d_qs), t(con1d_t)
{
  setup_con1d();
}


//------------CON1D solve_standalone -----------------------------------------------

void CON1D_ModelEval::solve_standalone() 
{
    //std::cout << "CON1D_ModelEval::solve_standalone called" << std::endl;
    take_time_step_con1d();
}


//------------CON1D is_converged ---------------------------------------------------

bool CON1D_ModelEval::is_converged() 
{
  con1d_residual();
  //std::cout << "CON1D_ModelEval R_max = " << VAR_R_MAX << std::endl;
    return( con1d_r_max < con1d_tol );
}


//------------CON1D finish ---------------------------------------------------------

CON1D_ModelEval::~CON1D_ModelEval()
{
   finish_con1d();
}


//------------CON1D residual -------------------------------------------------------

void CON1D_ModelEval::residual_con1d()
{
   con1d_residual();
}

//------------ get_time_step --------------------------------------------------------

double CON1D_ModelEval::get_time_step() const
{
   return con1d_dt;
}


//------------ get_max_time --------------------------------------------------------

double CON1D_ModelEval::get_max_time() const
{
   return con1d_tmax;
}

//------------ get_current_time ----------------------------------------------------

double CON1D_ModelEval::get_current_time() const
{
   return con1d_time;
}

//------------ update_time --------------------------------------------------------

void CON1D_ModelEval::update_time() 
{
   update_con1d();
}

