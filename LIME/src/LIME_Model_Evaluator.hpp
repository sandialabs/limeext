// ********************************************************************************************
// LIME 1.0:  Lightweight Integrating Multiphysics Environment for coupling codes, Version 1.0
// Copyright (c) 2012, Sandia National Laboratories
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted 
// provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright notice, this list of 
//      conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright notice, this list of 
//      conditions and the following disclaimer in the documentation and/or other materials  
//      provided with the distribution.
//    * Neither the name of Sandia National Laboratories nor the names of its contributors may 
//      be used to endorse or promote products derived from this software without specific prior 
//      written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// (Note: This is a BSD 3-Clause License.  For info. see www.opensource.org/licenses/BSD-3-Clause )
// --------------------------------------------------------------------------------------------
//
// LIME_Model_Evaluator.hpp
//
// Description: The LIME Model_Evaluator class is the starting point for wrapping a physics
// code so that it can participate in multi-physics coupling.  The LIME Problem_Manager requires
// each coupled physics derive from (or inherit) this as a base class.  Physics are then free
// to implement the Model_Evaluator capabilties that their particular codes are able to support.
//
// The Model_Evaluator constructor requires one argument, the LIME Problem_Manager, that owns
// and orchestrates this model evaluator's participation in a multi-physics coupling.  An optional
// descriptive name can be supplied if desired, this name will be used to output error messages
// originating from this Model_Evaluator.
//
// There are six broad capabilties in the Model_Evaluator class, most physics will only implement
// a few of the available features, the remaining features have reasonable defaults so that users
// can safely ignore all Model_Evaluator interfaces except for those features their physics actually
// support.  The key capabilities are time integration and stepping, solving physics, residual
// support, preconditioning, predictors and variable scaling, convergence notification and testing, 
// and parameter sensitivity response functions.  Many of these features are disabled by default
// and can be enabled by implementing the supports_feature() method and returning boolean true, to
// indicate that your physics is able to support this capability.  An example is the 
// supports_standalone_solve() feature.  If you physics is able to perform a standalone solve, then
// your model evaluator class should notify the Problem_Manager of this ability by returning true
// from this function.
//
// 
// ********************************************************************************************

#ifndef lime_model_evaluator_hpp
#define lime_model_evaluator_hpp

#include <map>
#include <set>
#include <vector>

#include "EpetraExt_ModelEvaluator.h" // base class

//-----------------------------------------------------------------------------
// 

class Epetra_Vector;

namespace LIME {

  class Problem_Manager;
  class Data_Transfer_Operator;

  typedef unsigned ID_TYPE;

class Model_Evaluator : public EpetraExt::ModelEvaluator
{
public:

  // The Problem_Manager notifies each Model_Evaluator in a coupled
  // simulation whether the step just completed converged or not.  This
  // enumeration is passed via the step_converged() method below.  If 
  // your physics cares about whether the step converged or not (e.g.
  // you want to do output at the end of a step and the output differs
  // based on whether the step converged or not), implement the function
  // and test the value of the passed enumeration.
  typedef enum {
      converged = 0
    , not_converged
  } convergence;

  Model_Evaluator(LIME::Problem_Manager & pm, string name = "Unnamed ME");

  virtual ~Model_Evaluator() {}

  void assign_id(ID_TYPE id)
  { my_id_ = id; }

  ID_TYPE get_id() const
  { return my_id_; }

  const string & get_name() const
  { return my_name_; }

  // for xLPR we need to reset model evaluators so they behave if 
  // if they're being invoked the first time.  This means resetting
  // internal time and step counters, and potentially lots of other
  // interesting behaviors (this to handle the UQ interface).
  virtual void reinitialize() {}

  // Set initial condition for solution vector incorporating parameters, etc.
  virtual void initializeSolution() {}

  // Get the owning Problem_Manager
  virtual const LIME::Problem_Manager & get_owner_pm() const
  { return owning_manager_; }
  virtual LIME::Problem_Manager & get_owner_pm()
  { return owning_manager_; }

  // Get the top-level Problem_Manager
  //virtual const LIME::Problem_Manager & get_const_top_pm() const;
  //virtual LIME::Problem_Manager & get_top_pm();

  // Inherited from EpetraExt::ModelEvaluator.
  virtual InArgs createInArgs() const { return InArgs(); }
  virtual OutArgs createOutArgs() const { return OutArgs(); }
  virtual Teuchos::RCP<const Epetra_Map> get_x_map() const
    { return Teuchos::null; }
  virtual Teuchos::RCP<const Epetra_Map> get_f_map() const
    { return Teuchos::null; }
  virtual Teuchos::RCP<const Epetra_Vector> get_x_init() const
    { return Teuchos::null; }
  virtual Teuchos::RCP<const Epetra_Map> get_g_map(int) const
    { return Teuchos::null; }
  virtual Teuchos::RCP<Epetra_Operator> create_W() const
    { return Teuchos::null; }
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
    { throw "*ERR* This is a base problem class which should have been inherited."; }

  // A workaround for setting state prior to doing elimination xfers
  virtual Teuchos::RCP<Epetra_Vector> get_x_state()
    { return Teuchos::null; }
  virtual Teuchos::RCP<Epetra_Vector> get_x_state() const
    { return Teuchos::null; }
  virtual void set_x_state(Teuchos::RCP<Epetra_Vector>)
    { }
  virtual Teuchos::RCP<Epetra_Vector> get_p_state()
    { return Teuchos::null; }
  virtual Teuchos::RCP<Epetra_Vector> get_p_state() const
    { return Teuchos::null; }

  // support for problem response functions
  virtual Teuchos::RefCountPtr<const Teuchos::Array<std::string> > get_p_names(int l) const
    { return Teuchos::null; }
  virtual Teuchos::RCP<Epetra_Vector> get_response_vec() const
    { return Teuchos::null; }

  // residual support
  virtual bool supports_residual() const { return false; }
  bool computed_residual() const { return m_computed_residual; }
  void computed_residual(bool b) const { m_computed_residual = b; }

  // Extensions to the EpetraExt::ModelEvaluator to support our multiphysics coupling strategies

  // Some time-stepping API methods
  virtual bool is_transient() const { return false; }
  virtual double get_max_time() const { return 0; }
  virtual double get_time_step() const { return 0; }
  virtual double get_current_time() const { return 0; }
  virtual unsigned int get_max_steps() const { return 0; }

  // LIME fixed point algorithm interface

  // Notify the Problem_Manager that your code can be solved
  // using a fixed point algorithm.
  virtual bool supports_standalone_solve() const { return false; }

  // You should implement these methods if you've notified the
  // Problem_Manager that you can implement a standalone solve.
  // using a fixed point algorithm.
  virtual void solve_standalone() { }

  // This is an overloaded method that allows for changes to time step during solve
  virtual bool solve_standalone(double & dt) 
  { 
    solve_standalone();
    return true;
  }

  // To determine whether a fixed point step has converged, the 
  // Problem_Manager will interrogate each Model_Evaluator in a
  // coupled simulation to find out whether this particular 
  // physics has converged or not.
  virtual bool is_converged() { return false; }

  // returns true if predictor updated state vector, false ow
  virtual void solve() { }
  virtual bool initialize_time_step() { return false; }

  // pass in solution whether converged or not
  virtual void finalize_solution(const Epetra_Vector &) {}

  // give problems that can't compute residuals opportunity to advance step
  virtual void update_time() {}

  // copy application solution into Epetra vector for PM
  virtual void copy_solution_state_into(Epetra_Vector&) const {}

  // The Problem_Manager will notify each Model_Evaluator whether
  // the step converged or not.  This routine is always called
  // whether a fixed point or non-linear solve was performed.
  virtual void step_converged(convergence) {}

  // Following a completed step and in preparation for the next step,
  // the Problem_Manager will determine the time step to be used by
  // all time-dependent physics.  This method notifies each of the
  // the time step they should use for their next solve.
  virtual void set_time_step(double dt) {}

  // Deprecated in favor of set_time_step       
  virtual void update_time_step(double dt)
  {
    cout << "WARNING: update_time_step(double) is deprecated, use set_time_step(double) instead." << endl;
    return set_time_step(dt);
  }

  // API for allowing predictors
  virtual bool supports_predictor() const { return false; }

  // some interfaces allowing physics-based preconditioning
  virtual bool supports_preconditioning() const { return false; }
  virtual bool compute_preconditioner(const Epetra_Vector& x) const { return false; }
  virtual bool apply_preconditioner(const Epetra_Vector& x, Epetra_Vector& y) const { return false; }

  // some interfaces allowing algebraic scaling
  virtual bool supports_right_scaling() { return false; }
  virtual bool supports_left_scaling() { return false; }
  virtual void compute_right_scaling_vector() { } // right (variable) scaling
  virtual void compute_left_scaling_vector() { } // left (equation) scaling

  // methods to support hierarchical solves
  unsigned num_my_children() const
  { return my_children_.size(); }

  const std::set<ID_TYPE> & get_children() const
  { return my_children_; }

  virtual void copy_incoming_transfers(const Problem_Manager & top_pm);
  virtual std::vector<Teuchos::RCP<Data_Transfer_Operator> > get_incoming_transfers() const;


 protected:

  ID_TYPE my_id_;
  const string my_name_;
  mutable bool m_computed_residual;
  LIME::Problem_Manager & owning_manager_;
  std::set<ID_TYPE> my_children_;
  std::multimap<ID_TYPE, Teuchos::RCP<Data_Transfer_Operator> > incoming_transfers_;
};

} // namespace LIME

#endif

