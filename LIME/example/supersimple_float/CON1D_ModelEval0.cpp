// ********************************************************************************************
//
// CON1D_ModelEval.cpp
// This is a LIME model evaluator for a stand-alone wrap of the super-simple con1d code
// 
// ********************************************************************************************
#include <iostream>
#include "CON1D_ModelEval0.hpp"

#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

// Data from fortran modules and fortran routines are referenced using an extern C 
// block as shown below. The "name mangling" conventions are taken care of automatically
// in LIME by #including the file LIME_fortran_mangling.h

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
const float CON1D_ModelEval0::gold_temps[] = { 523.4678, 491.3500, 383.0898 };
const float CON1D_ModelEval0::gold_temps_w_thinwall[] = { 536.6940, 507.2687, 406.1250 };

//------------CON1D_ModelEval C++ constructor ------------------------------------------

CON1D_ModelEval0::CON1D_ModelEval0(LIME::Problem_Manager & pm, const string & name) :
  Model_Evaluator(pm, name), qs(con1d_qs), t(con1d_t), trbc(con1d_tr)
{
  setup_con1d();
}

//------------CON1D_ModelEval C++ destructor -------------------------------------------

CON1D_ModelEval0::~CON1D_ModelEval0()
{
   finish_con1d();
}


//------------CON1D_ModelEval solve_standalone -----------------------------------------

void CON1D_ModelEval0::solve_standalone() 
{
    // std::cout << "CON1D_ModelEval0::solve_standalone called" << std::endl;
    take_time_step_con1d();
}


//------------ get_time_step --------------------------------------------------------

double CON1D_ModelEval0::get_time_step() const
{
   return con1d_dt;
}


//------------ get_max_time --------------------------------------------------------

double CON1D_ModelEval0::get_max_time() const
{
   return con1d_tmax;
}

//------------ get_current_time ----------------------------------------------------

double CON1D_ModelEval0::get_current_time() const
{
   return con1d_time;
}

//------------ update_time --------------------------------------------------------

void CON1D_ModelEval0::update_time() 
{
   update_con1d();
}


//------------CON1D is_converged ---------------------------------------------------

bool CON1D_ModelEval0::is_converged() 
{
    con1d_residual();
    // std::cout << "DEBUG CON1D_ModelEval0 R_max = " << con1d_r_max << std::endl;
    return( con1d_r_max < con1d_tol );
}
