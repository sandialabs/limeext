// ********************************************************************************************
//
// CON1D_ModelEval0_w_Resid.cpp
// This is an additional member function that inherits from the main LIME model evaluator 
// for a stand-alone wrap of the super-simple con1d code, and adds the functionality needed
// for passing a residual up to LIME for evaluation
// 
// ********************************************************************************************
#include <iostream>
#include "CON1D_ModelEval0_w_Resid.hpp"

#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

//-------------------------------------------------------------------------
// Additional data we wish to access in our wrapped application
//-------------------------------------------------------------------------

using LIME::Problem_Manager;

extern "C" {

  // fortran module data

  #define con1d_qs LIME_MANGLING_MODULE(con1d_mod, qs, CON1D_MOD, QS)
  extern float con1d_qs[3];

  #define con1d_t LIME_MANGLING_MODULE(con1d_mod, t, CON1D_MOD, T)
  extern float con1d_t[3];

  #define con1d_r LIME_MANGLING_MODULE(con1d_mod, r, CON1D_MOD, R)
  extern float con1d_r[3];

  #define con1d_dt LIME_MANGLING_MODULE(con1d_mod, dt, con1d_mod, dt)
  extern float con1d_dt;

  #define con1d_tmax LIME_MANGLING_MODULE(con1d_mod, tmax, CON1D_MOD, TMAX)
  extern float con1d_tmax;

  #define con1d_time LIME_MANGLING_MODULE(con1d_mod, time, CON1D_MOD, TIME)
  extern float con1d_time;

  #define con1d_r_max LIME_MANGLING_MODULE(con1d_mod, r_max, CON1D_MOD, R_MAX)
  extern float con1d_r_max;

  #define con1d_tol LIME_MANGLING_MODULE(con1d_mod, tol, CON1D_MOD, TOL)
  extern float con1d_tol;

  #define con1d_nvars LIME_MANGLING_MODULE(con1d_mod, n_vars, CON1D_MOD, N_VARS)
  extern int con1d_nvars;

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

//------------CON1D_w_Resid constructor ----------------------------------------------

CON1D_ModelEval0_w_Resid::CON1D_ModelEval0_w_Resid(Problem_Manager & pm, const string & name) :
  CON1D_ModelEval0(pm, name), r(con1d_r)
{
  // Create an Epetra Map decribing the data layout for both the
  // unknowns array and residuals array.  This is trivial in serial
  // but contains distribution info in parallel.
  epetra_map_ = Teuchos::rcp(new Epetra_Map(con1d_nvars, 0, pm.Comm()));

  // Create a state vector conforming to the data layout map
  // This vector gets used to create conformal solver objects by the solvers
  // as well as to convey initial values for time stepping and nonlinear
  // iterations
  me_interface_soln_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
}

// ------------CON1D_w_Resid get_x_init --------------------------------------------------

Teuchos::RCP<const Epetra_Vector>
CON1D_ModelEval0_w_Resid::get_x_init() const
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

// ------------CON1D_w_Resid get_x_state --------------------------------------------------

Teuchos::RCP<Epetra_Vector>
CON1D_ModelEval0_w_Resid::get_x_state()
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

// ------------CON1D_w_Resid get_x_state --------------------------------------------------

Teuchos::RCP<Epetra_Vector>
CON1D_ModelEval0_w_Resid::get_x_state() const
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

//------------CON1D_w_Resid initializeSolution -----------------------------------------------

void
CON1D_ModelEval0_w_Resid::initializeSolution()
{
  // Copy initial values from our wrapped application
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
}

//------------CON1D_w_Resid createInArgs -----------------------------------------------

EpetraExt::ModelEvaluator::InArgs 
CON1D_ModelEval0_w_Resid::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  // Register our identify
  inArgs.setModelEvalDescription(my_name_);

  // Signal that we can do calculations with incoming state x
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}

//------------CON1D_w_Resid createOutArgs -----------------------------------------------

EpetraExt::ModelEvaluator::OutArgs 
CON1D_ModelEval0_w_Resid::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  // Register our identify - consistent with createInArgs
  outArgs.setModelEvalDescription(my_name_);

  // Signal that we can compute a residual vector
  outArgs.setSupports(OUT_ARG_f, true);

  return outArgs;
}

//------------CON1D_w_Resid evalModel -----------------------------------------------

void 
CON1D_ModelEval0_w_Resid::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  // Create a "View" of the incoming solution state, x
  //  This avoids a copy of data and becomes important for real-sized problems
  const Epetra_Vector x(View, *(inArgs.get_x().get()), 0);

  //cout << "CON1D_ModelEval0_w_Resid::evalModel called. state:\n";

  // Copy values into our wrapped application's data array
  for( int i = 0; i < x.MyLength(); ++i )
  {
    //cout << "x[" << i << "] = " << x[i] << endl;
    (this->t)[i] = x[i];
  }

  if( outArgs.get_f().get() ) // A non-NULL f-vector signals a residual fill request
  {
    // Get a reference to the vector we will populate
    Epetra_Vector & f = *(outArgs.get_f().get());

    // compute our application's nonlinear residual
    residual_con1d_();

    // Copy application residual array into f-vector provided by outArgs
    for( int i = 0; i < f.MyLength(); ++i )
    {
      f[i] = (this->r)[i];
      //cout << "f[" << i << "] = " << f[i] << endl;
    }
  }

  // We don't yet support this and respond to such a request with an error
  if( outArgs.get_W().get() ) // Signals either a computeJacobian or a computePreconditioner
    throw std::runtime_error("CON1D_ModelEval0_w_Resid::evalModel : Jacobian matrix support not available.");
}
