// ********************************************************************************************
//
// THINWALL_ModelEval0_w_Resid.cpp
// This is an additional member function that inherits from the main LIME model evaluator 
// for a stand-alone wrap of the super-simple thinwall code, and adds the functionality needed
// for passing a residual up to LIME for evaluation
// 
// ********************************************************************************************
#include <iostream>
#include "THINWALL_ModelEval0_w_Resid.hpp"

#include "LIME_Problem_Manager.hpp"
#include "LIME_fortran_mangling.h"

using LIME::Problem_Manager;

//-------------------------------------------------------------------------
// Additional data we wish to access in our wrapped application
//-------------------------------------------------------------------------

extern "C" {

  // fortran module data

  #define thinwall_r LIME_MANGLING_MODULE(thinwall_mod, r, THINWALL_MOD, R)
  extern float thinwall_r[1];

  #define thinwall_nvars LIME_MANGLING_MODULE(thinwall_mod, n_vars, THINWALL_MOD, N_VARS)
  extern int thinwall_nvars;

  #define residual_thinwall LIME_MANGLING_GLOBAL(residual_thinwall, RESIDUAL_THINWALL)
   void residual_thinwall();
 }

//------------THINWALL_w_Resid constructor ----------------------------------------------

THINWALL_ModelEval0_w_Resid::THINWALL_ModelEval0_w_Resid(Problem_Manager & pm, const string & name) :
  THINWALL_ModelEval0(pm, name), r(thinwall_r)
{
  // Create an Epetra Map decribing the data layout for both the
  // unknowns array and residuals array.  This is trivial in serial
  // but contains distribution info in parallel.
  epetra_map_ = Teuchos::rcp(new Epetra_Map(thinwall_nvars, 0, pm.Comm()));

  // Create a state vector conforming to the data layout map
  // This vector gets used to create conformal solver objects by the solvers
  // as well as to convey initial values for time stepping and nonlinear
  // iterations
  me_interface_soln_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
}

// ------------THINWALL_w_Resid get_x_init --------------------------------------------------

Teuchos::RCP<const Epetra_Vector>
THINWALL_ModelEval0_w_Resid::get_x_init() const
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

// ------------THINWALL_w_Resid get_x_state --------------------------------------------------

Teuchos::RCP<Epetra_Vector>
THINWALL_ModelEval0_w_Resid::get_x_state()
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

// ------------THINWALL_w_Resid get_x_state --------------------------------------------------

Teuchos::RCP<Epetra_Vector>
THINWALL_ModelEval0_w_Resid::get_x_state() const
{
   // Copy current values from our wrapped application into LIME array
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
   return me_interface_soln_;
}

// ------------THINWALL_ModelEval0_w_Resid set_x_state --------------------------------------------------

void
THINWALL_ModelEval0_w_Resid::set_x_state(Teuchos::RCP<Epetra_Vector> x)
{ 
  // Copy incoming state into our wrapped application
  for( int i = 0; i < x->MyLength(); ++i )
    (this->t)[i] = (*x)[i];
}


//------------THINWALL_w_Resid initializeSolution -----------------------------------------------

void
THINWALL_ModelEval0_w_Resid::initializeSolution()
{
  // Copy initial values from our wrapped application
  for( int i = 0; i < me_interface_soln_->MyLength(); ++i )
    (*me_interface_soln_)[i] = (this->t)[i];
}

//------------THINWALL_w_Resid createInArgs -----------------------------------------------

EpetraExt::ModelEvaluator::InArgs 
THINWALL_ModelEval0_w_Resid::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  // Register our identify
  inArgs.setModelEvalDescription(my_name_);

  // Signal that we can do calculations with incoming state x
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}

//------------THINWALL_w_Resid createOutArgs -----------------------------------------------

EpetraExt::ModelEvaluator::OutArgs 
THINWALL_ModelEval0_w_Resid::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  // Register our identify - consistent with createInArgs
  outArgs.setModelEvalDescription(my_name_);

  // Signal that we can compute a residual vector
  outArgs.setSupports(OUT_ARG_f, true);

  return outArgs;
}

//------------THINWALL_w_Resid evalModel -----------------------------------------------

void 
THINWALL_ModelEval0_w_Resid::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  // Create a "View" of the incoming solution state, x
  //  This avoids a copy of data and becomes important for real-sized problems
  const Epetra_Vector x(View, *(inArgs.get_x().get()), 0);

  //cout << "THINWALL_ModelEval0_w_Resid::evalModel called. state:\n";

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
    residual_thinwall();

    // Copy application residual array into f-vector provided by outArgs
    for( int i = 0; i < f.MyLength(); ++i )
    {
      f[i] = (this->r)[i];
      //cout << "f[" << i << "] = " << f[i] << endl;
    }
  }

  // We don't yet support this and respond to such a request with an error
  if( outArgs.get_W().get() ) // Signals either a computeJacobian or a computePreconditioner
    throw std::runtime_error("THINWALL_ModelEval0_w_Resid::evalModel : Jacobian matrix support not available.");
}
