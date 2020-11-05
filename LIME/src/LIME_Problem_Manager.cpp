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
// LIME_Problem_Manager.cpp
//
// The Problem Manager is the central orchestrating software component of LIME
// 
// ********************************************************************************************

//-----------------------------------------------------------------------------
// Include Directives

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "LIME_Problem_Manager.hpp"
#include "LIME_Model_Evaluator.hpp"
#include "LIME_Problem_As_Operator.hpp"
#include "LIME_FixedPoint_Accelerators.hpp"
#include "LIME_Data_Transfer_Operator.hpp"
#include "LIME_PrePostOperator.hpp"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

#ifdef HAVE_LIME_JFNK
#include "NOX.H"
#include "NOX_Epetra.H"

// Headers needed for Coloring
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_CrsMatrixIn.h"
#endif
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_ConditionNumber.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#endif

// Headers for reading/writing configuration parameter lists
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using std::binary_search;
using std::cout;
using std::endl;
using std::lower_bound;
using std::numeric_limits;
using std::ofstream;
using std::runtime_error;
using std::vector;

using namespace EpetraExt;

using LIME::Model_Evaluator;
#ifdef HAVE_LIME_JFNK
using LIME::PrePostOperator;
#endif

namespace LIME {

//-----------------------------------------------------------------------------
// Problem_Manager C++ constructor

Problem_Manager::Problem_Manager(const Epetra_Comm& comm, bool verbose,
  bool use_stdout) :
  comm_(comm),
  m_clp(false),
  m_verbose(verbose),
  m_setupXMLfilename(""),
  m_paramsXMLfilename(""),
  use_predictor_(true),
  my_hierarchy_(0),
  num_total_dofs_(0),
  su_(fixed_point),
  transfer_type(fixed_point_xfer),
  register_complete_called_(false),
  pl_dir_name_("ParamLists"),
  out_stream_(Teuchos::rcp(new ofstream("solver.log", std::ios_base::out))),
  fp_solve_mode_(jacobi),
  max_fp_iters_(1), // do a single step by default
  fp_abs_tolerance_(0.0),
  max_time_steps_(numeric_limits<unsigned int>::max()),
  residual_active_(false),
  problem_is_transient(false),
  stateful_transfers(false),
  m_use_stdout(use_stdout)
{
  // Go through motions to create a Parameter List subdir, even if it already exists
  mkdir(pl_dir_name_.c_str(), 0777 );

  // Add default command line options
  add_command_line_options();
}

//-----------------------------------------------------------------------------
// Problem_Manager C++ destructor

Problem_Manager::~Problem_Manager()
{
  all_registered_transfers_.clear();
  standard_transfers_.clear();
}

//-----------------------------------------------------------------------------
// notify use of stateful transfers

void Problem_Manager::problem_has_stateful_transfers() 
{
  stateful_transfers = true;
}

//-----------------------------------------------------------------------------
// This query method returns the number of problems at this hierarchy level

unsigned
Problem_Manager::num_problems(unsigned hier) const
{
  if( 0 != my_hierarchy_ )
    throw std::runtime_error("Call to num_problems should only be on top-level Problem_Manager.");

  map<unsigned, unsigned>::const_iterator citer = num_hierarchy_problems_.find(hier);
  if( num_hierarchy_problems_.end() == citer )
    return 0;

  return citer->second;
}

//-----------------------------------------------------------------------------
// Increment the problem count for this problem hierarchy

void
Problem_Manager::increment_num_problems(unsigned hier)
{
  if( 0 != my_hierarchy_ )
    throw std::runtime_error("Call to increment_num_problems should only be on top-level Problem_Manager.");

  num_hierarchy_problems_[hier]++;
}

//-----------------------------------------------------------------------------
// This method calculates the next available unique problem id based on 
// all existing registered problems within a multi-physics problem
// hierarchy.  A return value of true is returned if sucessful, false
// if not.

bool
Problem_Manager::get_next_problem_id(unsigned hier, ID_TYPE & value) const
{
  if( num_problems(hier) >= LIME_MAX_PROBLEMS_PER_HIERARCHY )
    return false;

  value = hier*LIME_MAX_HIERARCHIES + num_problems(hier) + 1;

  return true;
}

//-----------------------------------------------------------------------------
// This method adds a previously created problem to the problem hierarchy.
// The order in which this is called matters for Seidel type fixed point
// in that it determines the order of solves and data trasnfers.

ID_TYPE 
Problem_Manager::add_problem(Teuchos::RCP<Model_Evaluator> problem)
{
  ID_TYPE id = problem->get_id();
  if( 0 == provisional_problem_ids_.count(id) )
    throw std::runtime_error("Invalid problem id.  Did you create the problem using a different Manager?");

  // problem is no longer provisional
  provisional_problem_ids_.erase(id);

  string problemName = problem->createInArgs().modelEvalDescription();
  problems_[id] = problem;
  problem_ids_.push_back(id);
  id_to_name_[id] = problemName;
  name_to_id_[problemName] = id;
  increment_num_problems(my_hierarchy_);

  return id;
}

//-----------------------------------------------------------------------------
// Get a RCP reference to the problem given its id

const Teuchos::RCP<Model_Evaluator> & 
Problem_Manager::get_problem(ID_TYPE id) const
{
  map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::const_iterator iter = problems_.find(id);
  if( iter == problems_.end() )
    throw std::runtime_error("ERROR: No valid problem (id=" + convert_to_string(id) + ") registered with Problem_Manager !!");

  return iter->second;
}

//-----------------------------------------------------------------------------
// Get a RCP reference to the problem given its id

Teuchos::RCP<Model_Evaluator> & 
Problem_Manager::get_problem(ID_TYPE id)
{
  map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::iterator iter = problems_.find(id);
  if( iter == problems_.end() )
    throw std::runtime_error("No valid problem (id=" + convert_to_string(id) + ") registered with Problem_Manager !!");

  return iter->second;
}

//-----------------------------------------------------------------------------
// Get a list of all problem ids that are not Hierarchy_Problem_Managers but
// are Model_Evaluators connected to a wrapped application.  This amounts to 
// getting all the leaves of a problem hierarchy from this level down.

vector<ID_TYPE>
Problem_Manager::get_all_descendants() const
{
  vector<ID_TYPE> descendants;

  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    if( problem.num_my_children() )
      std::copy(problem.get_children().begin(), problem.get_children().end(), std::back_inserter(descendants));
    else
      descendants.push_back(problem_ids_[i]);
  }

  return descendants;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::publish_problem_parameters(Teuchos::RCP<const Epetra_Vector> invec) const
{
  if( !invec->Map().SameAs(m_model_p_vec->Map()) )
    throw std::runtime_error("Attempted to set parameters with incorrectly sized vector.");

  multimap<int, pair<int,int> >::const_iterator aiter;
  for( aiter = m_model_p_mappings.begin(); aiter != m_model_p_mappings.end(); ++aiter )
  {
    const Model_Evaluator & problem  = *get_problem(aiter->second.first);
    Teuchos::RCP<Epetra_Vector> pvec = problem.get_p_state();
    (*pvec)[aiter->second.second] = (*invec)[aiter->first];
    //EpetraExt::ModelEvaluator::InArgs inArgs = problem.createInArgs();
    //EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
    //inArgs.set_p(0, pvec);
    ////(*get_problem(aiter->second.first)->get_p_init(0))[aiter->second.second] = (*m_model_p_vec)[aiter->first];
    //problem.evalModel(inArgs, outArgs);
  }
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::add_command_line_options()
{
  m_clp.setOption( "verbose", "no-verbose", &m_verbose, "Verbosity on or off." );
  m_clp.setOption( "setupxml", &m_setupXMLfilename, "Filename for XML used to setup Problem_Manager. Default is \"Problem_Manager_setup.xml\"" );
}

//-----------------------------------------------------------------------------

double 
Problem_Manager::determine_time_step() const
{
  double time_step = numeric_limits<double>::max();

  // If overall problem is steady-state, setup time loop for single pass
  if( !problem_is_transient )
     time_step = 1.0; // should this change to numeric_limits<double>::max() ?

  // If transient, negotiate a most-conservative time step size
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);

    // Skip to next physics code if this one does not provide temporal support
    if( !problem.is_transient() ) continue;

    // Get this physics code's time and step size related values
    double dt = problem.get_time_step();

    // Error check: time step must be positive
    if (dt <= 0)
      throw runtime_error("Error: dt must be positive (greater than zero)");

    // Reconciliation of time_step
    time_step = std::min(time_step, dt);
  }

  // Do a global reduction on the time step before updating problems
  // with agreed upon time step
  double result = 0;     
  comm_.MinAll(&time_step, &result, 1); // get MPI communicator off Epetra

  return result;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::post_problem_response(const string & g_name, Teuchos::RCP<Model_Evaluator> model, Teuchos::RCP<Epetra_Vector> g_vec)
{
  // Check to ensure that any parameters vectors already registered with the same name
  // have compatible data layout
  if( m_responses.find(g_name) != m_responses.end() )
    throw std::runtime_error("Attempted to register a problem response, \""+g_name+"\" that already exists.");

  pair<Teuchos::RCP<Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > pair1(model, g_vec);
  m_responses.insert( pair<string, pair<Teuchos::RCP<Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > >(g_name, pair1) );
}

//-----------------------------------------------------------------------------

#ifdef HAVE_LIME_JFNK
const Teuchos::ParameterList &
Problem_Manager::get_jfnk_solver_params() const
{
  if( (su_ != jfnk) || Teuchos::is_null(composite_nox_solver_) )
    throw std::runtime_error("Problem_Manager has no valid JFNK solver.");

  return composite_nox_solver_->getList();
}
#endif

//-----------------------------------------------------------------------------

void 
Problem_Manager::compute_problem_responses() const
{
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
    if( 0 < outArgs.Ng() )
    {
      EpetraExt::ModelEvaluator::InArgs inArgs = problem.createInArgs();
      inArgs.set_x(problem.get_x_state());
      const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> g_vec(problem.get_response_vec(), EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT);
      outArgs.set_g(0, g_vec);
      problem.evalModel(inArgs, outArgs);
    }
  }
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::register_complete()
{
  if(problems_.empty())
    throw std::runtime_error("No problems registered with Problem_Manager !!");

  // Do configurations that don't depened on the coupling algorithm
  if( "" == m_setupXMLfilename )
    m_setupXMLfilename = "Problem_Manager_setup.xml";
  Teuchos::RCP<Teuchos::ParameterList> setup_list = readParameterList(m_setupXMLfilename);
  if( !Teuchos::is_null(setup_list) )
  {
    if( setup_list->isParameter("Max Time Steps") ) {
      max_time_steps_ = setup_list->get<int>("Max Time Steps");
    }

    if( setup_list->isParameter("Use Predictor") )
      use_predictor_ = setup_list->get<bool>("Use Predictor");

    if( setup_list->isParameter("Transfer Type") ) {
      string s = setup_list->get<string>("Transfer Type");
      if (s == "Successive") {
        transfer_type = successive_xfer;
      }
      else if (s == "FixedPoint") {
        transfer_type = fixed_point_xfer;
      }
      else {
        string t = "Error: invalid transfer type option " + s;
        t += ", valid choices are FixedPoint or Successive.";
        throw std::runtime_error(t);
      }
    }

    if( setup_list->isParameter("Solver Strategy") ) {
      string s = setup_list->get<string>("Solver Strategy");
      if (s == "JFNK") {
        su_ = jfnk;
      }
      else if (s == "SWITCHING") {
        su_ = switching;
      }
      else if (s != "FIXED_POINT") {
        string t = "Error: invalid solver type option " + s;
        t += ", valid choices are JFNK, FIXED_POINT, and SWITCHING.";
        throw std::runtime_error(t);
      }
    }
  }

  if( "" != m_paramsXMLfilename )
    configure_model_parameters(m_paramsXMLfilename);

  int maxGID = 0;
  int compositeGID = 0;
  vector<int> composite_global_ids;

  // Iterate over each problem and construct the necessary objects
  // The assumption here is that all problems that support IN_ARG_x provide 
  // active degrees of freedom to the composite coupled solution vector.
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    problem.initializeSolution();
    if (problem.supports_residual()) { residual_active_ = true; }
    if (problem.is_transient()) { problem_is_transient = true; }
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
    if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
    {
      // Create index mapping for this problem into the composite problem
      const Epetra_Vector & x_init = *(problem.get_x_init());
      problem_LID_to_composite_GID_map_[problem_ids_[i]] = Teuchos::rcp( new Epetra_IntVector(x_init.Map()) );
      Epetra_IntVector & indices = *(problem_LID_to_composite_GID_map_[problem_ids_[i]]);

      for( int j = 0; j < x_init.MyLength(); ++j ) 
      {
        compositeGID = maxGID + x_init.Map().GID(j);
        composite_global_ids.push_back( compositeGID );
        indices[j] = compositeGID;
      }

      maxGID += x_init.Map().MaxAllGID() + 1;

      problem_residual_data_[problem_ids_[i]]   = Teuchos::rcp( new Epetra_Vector(x_init.Map()) );
    }
  }

  if( !composite_global_ids.empty() )
    composite_map_ = Teuchos::rcp( new Epetra_Map(-1, (int) composite_global_ids.size(), &composite_global_ids[0], 0, comm_) );
  else
  {
    int zero = 0;
    composite_map_ = Teuchos::rcp( new Epetra_Map(-1, 0, &zero, 0, comm_) );
  }
  composite_solution_ = Teuchos::rcp( new Epetra_Vector(*composite_map_) );
  composite_residual_ = Teuchos::rcp( new Epetra_Vector(*composite_map_) );

  refresh_solution_state();

  // For top-level manager, verify that all transfers are between leaves of the problem hierarchy
  if( my_hierarchy_ == 0 )
  {
    typedef multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >::const_iterator MMciter;
    MMciter iter = all_registered_transfers_.begin(),
            iter_end = all_registered_transfers_.end();
    for( ; iter_end != iter; ++iter )
    {
      if( Teuchos::is_null(iter->second->source()) ) // skip old-style registered transfers
        continue;
      else if( (0 != iter->second->source()->num_my_children()) ||
          (0 != iter->second->target()->num_my_children())   )
        throw runtime_error("Transfers can only be between leaves of the coupled problem hierarchy (non-aggregates).");
    }

    // Test leaf problem status
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      //echo_tee("Problem: "+convert_to_string(problem.get_id())+" ("+convert_to_string(problem_ids_[i])+") has "+convert_to_string(problem.get_children().size())+" children.");
      problem.copy_incoming_transfers(*this);
    }
  }

  // initialize solver
  setup_nonlinear_solve();

  register_complete_called_ = true;

  return;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::configure_model_parameters(const string & xml_filename)
{
  Teuchos::RCP<Teuchos::ParameterList> params_list = readParameterList(xml_filename);
  if( !Teuchos::is_null(params_list) )
  {
    vector<ID_TYPE> active_parameter_prob_ids;
    vector<Teuchos::Array<double>::size_type> active_parameter_p_indices;
    vector<string> active_parameter_names;
    map<string, pair<int,size_t> > active_parameter_info;

    const Teuchos::ParameterList & allXMLspecs = params_list->sublist("Expose ME Parameters");
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      // currently only a single p-vector is supported
      Teuchos::RCP<const Teuchos::Array<std::string> > avail_prob_params = problem.get_p_names(0);
      if( !Teuchos::is_null(avail_prob_params) )
      {
        if( allXMLspecs.isSublist(problem.get_name()) )
        {
          vector<string> avail_param_names_vec = avail_prob_params->toVector();
          const Teuchos::ParameterList & meXMLspecs = allXMLspecs.sublist(problem.get_name());
          if( !meXMLspecs.isParameter("p0") )
            throw std::runtime_error("Sublist \""+problem.get_name()+"\" must have an entry for \"p0\".");

          const Teuchos::Array<string>
            readStringArray = Teuchos::getArrayFromStringParameter<string>(meXMLspecs,"p0");
          std::vector<string> vecString = readStringArray.toVector();
          for( size_t j = 0; j < vecString.size(); ++j )
          {
            bool found = false;
            size_t pn;
            for( pn = 0; pn < avail_param_names_vec.size(); ++pn )
              if( avail_param_names_vec[pn] == vecString[j] )
              {
                found = true;
                break;
              }
            if( !found )
              throw std::runtime_error("Problem \"" + problem.get_name() 
                  + "\" does not have parameter \"" + vecString[j] + "\".");

            active_parameter_prob_ids.push_back(problem_ids_[i]);
            active_parameter_p_indices.push_back(pn);
            string param_name = problem.get_name() + "::" + vecString[j];
            active_parameter_names.push_back(param_name);
            active_parameter_info[param_name] = std::make_pair<int,size_t>(problem_ids_[i], pn);
          }
        }
      }
    }

    // Now setup parameter aliases

    // First make a temporary copy of all active parameter names for use in making associations
    set<string> active_param_names;
    for( size_t i = 0; i < active_parameter_names.size(); ++i )
      active_param_names.insert(active_parameter_names[i]);

    m_model_p_names = Teuchos::rcp(new Teuchos::Array<string>);
    int num_unique_names = 0;
    const string alias_block_name("ME Parameters Aliases");
    if( params_list->isSublist(alias_block_name) )
    {
      const Teuchos::ParameterList & aliasXMLspecs = params_list->sublist(alias_block_name);
      Teuchos::ParameterList::ConstIterator pIter = aliasXMLspecs.begin(),
        pEnd = aliasXMLspecs.end();
      for( ; pEnd != pIter; ++pIter )
      {
        const string alias_name = pIter->first;
        if( !aliasXMLspecs.isSublist(alias_name) )
          throw std::runtime_error("Error for Parameter \"" + alias_name + "\". \n" + 
              "\tOnly sublists should exist within \"" + alias_block_name + "\" ParameterList.");

        const Teuchos::ParameterList & alias_set =  aliasXMLspecs.sublist(alias_name);
        Teuchos::ParameterList::ConstIterator aliasIter = alias_set.begin(),
          aliasEnd = alias_set.end();
        for( ; aliasEnd != aliasIter; ++aliasIter )
        {
          if( active_param_names.end() == active_param_names.find(aliasIter->first) )
            throw std::runtime_error("Could not find a parameter entry for \"" + aliasIter->first + "\"");

          //m_model_parameter_aliases.insert( pair<string, string>(alias_name, aliasIter->first) );
          m_model_p_mappings.insert(pair<int,pair<int,int> >(num_unique_names, active_parameter_info[aliasIter->first]));
          active_param_names.erase(aliasIter->first);
        }
        num_unique_names++;
        m_model_p_names->push_back(alias_name);
      }
    }

    set<string>::iterator siter = active_param_names.begin(),
      siter_end = active_param_names.end();
    for( ; siter_end != siter; ++ siter )
    {
      m_model_p_mappings.insert(pair<int,pair<int,int> >(num_unique_names, active_parameter_info[*siter]));
      num_unique_names++;
      m_model_p_names->push_back(*siter);
    }

    // Finalize parameter setup
    m_model_p_map = Teuchos::rcp(new Epetra_LocalMap(num_unique_names, 0, comm_));
    m_model_p_vec = Teuchos::rcp(new Epetra_Vector(*m_model_p_map), false);

    // Fill initial values from appropriate MEs
    // Note: We just loop over all active entries so that wrt aliased parameters, the last assignment wins
    multimap<int, pair<int,int> >::const_iterator aiter;
    for( aiter = m_model_p_mappings.begin(); aiter != m_model_p_mappings.end(); ++aiter )
      (*m_model_p_vec)[aiter->first] = (*get_problem(aiter->second.first)->get_p_init(0))[aiter->second.second];
  }
}

//-----------------------------------------------------------------------------

void
Problem_Manager::refresh_solution_state() const
{
  copy_problems_into_composite(SOLUTION_X);
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_problems_into_composite( VECTOR_TYPE vec_type ) const
{
  if( RESIDUAL == vec_type )
    throw "Problem_Manager ERROR : Problem_Manager::copy_problem_into_composite not yet supported for RESIDUAL.";

  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
      copy_problem_into_composite(problem_ids_[i], vec_type);
  }

  return;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_problem_into_composite( ID_TYPE problem_id, VECTOR_TYPE vec_type ) const
{
  Model_Evaluator & problem  = *(problems_[problem_id]);
  EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();

  Teuchos::RefCountPtr<const Epetra_Vector> problem_vector_ptr   = Teuchos::null;
  Teuchos::RefCountPtr<Epetra_Vector>       composite_vector_ptr = Teuchos::null;

  if( SOLUTION_X  == vec_type )
  {
    problem_vector_ptr = problem.get_x_init();
    if( Teuchos::is_null(problem_vector_ptr) )
      return;

    composite_vector_ptr = composite_solution_;
  }
  else if( RESIDUAL == vec_type )
  {
    problem_vector_ptr = problem_residual_data_[problem_id];
    composite_vector_ptr = composite_residual_;
  }
  else
    throw "Problem_Manager ERROR : Problem_Manager::copy_problem_into_composite, unsupported VECTOR_TYPE.";

  const Epetra_Vector & problem_vector   = *(problem_vector_ptr);
        Epetra_Vector & composite_vector = *(composite_vector_ptr);
  if (!problem_LID_to_composite_GID_map_.empty()) {
    const Epetra_IntVector & indices = *(problem_LID_to_composite_GID_map_[problem_id]);
    for( int i = 0; i < indices.MyLength(); ++i ) 
      composite_vector[indices[i]] = problem_vector[i]; 
  }
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_problem_into_composite( ID_TYPE problem_id, VECTOR_TYPE vec_type, Epetra_Vector * p_composite_vector ) const
{
  if( !composite_map_->SameAs(p_composite_vector->Map()) )
    throw std::runtime_error("Problem_Manager::copy_problem_into_composite, destination vector does not conform to composite data map.");

  Model_Evaluator & problem  = *(get_problem(problem_id));
  EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();

  Teuchos::RefCountPtr<const Epetra_Vector> problem_vector_ptr   = Teuchos::null;

  if( SOLUTION_X  == vec_type )
    problem_vector_ptr = problem.get_x_init();
  else if( RESIDUAL == vec_type )
    problem_vector_ptr = problem_residual_data_[problem_id];
  else
    throw std::runtime_error("Problem_Manager::copy_problem_into_composite, unsupported VECTOR_TYPE.");

  const Epetra_Vector & problem_vector   = *(problem_vector_ptr);
        Epetra_Vector & composite_vector = *(p_composite_vector);
  map<ID_TYPE, Teuchos::RCP<Epetra_IntVector> >::const_iterator iter =  problem_LID_to_composite_GID_map_.find(problem_id);
  const Epetra_IntVector & indices = *(iter->second);

  for( int i = 0; i < indices.MyLength(); ++i ) 
    composite_vector[indices[i]] = problem_vector[i]; 

  return;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_vector_into_composite_solution( ID_TYPE problem_id, const Epetra_Vector & vec )
{
  const Model_Evaluator & problem  = *(problems_[problem_id]);
  if( !problem.get_x_map()->SameAs(vec.Map()) )
    throw std::runtime_error("copy_vector_into_composite_solution has non-conforming vector maps.");

  const Epetra_IntVector & indices = *(problem_LID_to_composite_GID_map_[problem_id]);
  Epetra_Vector & x_composite = *(composite_solution_);

  for( int i = 0; i < indices.MyLength(); ++i ) 
    x_composite[indices[i]] = vec[i]; 

  return;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_composite_solution_into_problem_data() const
{
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
     copy_composite_solution_into_problem_data(problem_ids_[i]);
  }

  return;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::copy_composite_solution_into_problem_data(ID_TYPE probId) const
{
  const Epetra_Vector & x_composite = *(composite_solution_);

  Teuchos::RCP<Epetra_Vector> p_x_state = Teuchos::null;
  if( !Teuchos::is_null(get_problem(probId)->get_x_state()) )
    p_x_state = get_problem(probId)->get_x_state();
  map<ID_TYPE, Teuchos::RCP<Epetra_IntVector> >::const_iterator iter =  problem_LID_to_composite_GID_map_.find(probId);
  const Epetra_IntVector & indices = *(iter->second);
  if( !Teuchos::is_null(p_x_state) )
  {
    for( int i = 0; i < indices.MyLength(); ++i ) 
      (*p_x_state)[i] = x_composite[indices[i]];
    get_problem(probId)->set_x_state(p_x_state);
  }

  return;
}

//-----------------------------------------------------------------------------

const Epetra_Vector * 
Problem_Manager::create_view_into_composite_vector(ID_TYPE probid, VECTOR_TYPE type)
{
  Teuchos::RefCountPtr<Epetra_Vector> composite_vector_ptr = Teuchos::null;

  if( SOLUTION_X  == type )
    composite_vector_ptr = composite_solution_;
  else if( RESIDUAL == type )
    composite_vector_ptr = composite_residual_;
  else
    throw "Problem_Manager ERROR : Problem_Manager::copy_problem_into_composite, unsupported VECTOR_TYPE.";

  map<ID_TYPE, Teuchos::RCP<Epetra_IntVector> >::const_iterator iter =  problem_LID_to_composite_GID_map_.find(probid);
  const Epetra_IntVector & indices = *(iter->second);

  const Epetra_Vector * vec_view = new Epetra_Vector(View, get_problem(probid)->get_x_init()->Map(), &((*composite_vector_ptr)[indices[0]]));

  return vec_view;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::output_status( ostream & os ) 
{ 
#if 0
  map<int, Teuchos::RCP<Model_Evaluator> >::const_iterator 
           iter = problems_.begin(),
      iter_end = problems_.end();

  os << endl << endl << "\t\t********************************"   << endl;
  os                 << "\t\t*******  PROBLEM SUMMARY  ******"   << endl;
  os                 << "\t\t********************************"   << endl;
  os << endl << endl;

  // Loop over each problem being managed and output its dependencies
  for( ; iter_end != iter; ++iter )
  {
    const Model_Evaluator & problem = *(*iter).second;
    EpetraExt::ModelEvaluator::InArgs   inArgs  = problem.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs  outArgs = problem.createOutArgs();

    os << "\tProblem \"" << id_to_name_[iter->first] << "\" (" << (*iter).first
         << ")  supports:" << endl << endl;
    
    os << "\t\tIN_ARG_x_dot --> "       << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot) << endl;
    os << "\t\tIN_ARG_x --> "           << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) << endl;
    os << "\t\tIN_ARG_x_dot_poly  --> " << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot_poly) << endl;
    os << "\t\tIN_ARG_x_poly  --> "     << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_poly) << endl;
    os << "\t\tIN_ARG_t --> "           << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_t) << endl;
    os << "\t\tIN_ARG_alpha --> "       << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha) << endl;
    os << "\t\tIN_ARG_beta --> "        << inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_beta) << endl;
    os << endl;

    os << "\t\tOUT_ARG_f --> "          << outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_f) << endl;
    os << "\t\tOUT_ARG_W --> "          << outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W) << endl;
    os << "\t\tOUT_ARG_f_poly  --> "    << outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_f_poly) << endl;
    os << endl << endl;
  }
#endif
}

//-----------------------------------------------------------------------------

bool 
Problem_Manager::verify_consistent_current_time( double time )
{ 
  bool is_consistent = true;
  double time_eps = 10*numeric_limits<double>::epsilon();

  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);

    // this model evaluator does not implement any temporal support
    if( !problem.is_transient() ) continue;

    // get their time and step size
    double current_time = problem.get_current_time();
    if( fabs(current_time-time) > time_eps ) {
      string msg = "inconsistent time for problem " + convert_to_string(id_to_name_[problem_ids_[i]]) 
                   + " at time " + convert_to_string(current_time)
                   + "\nexpected problem time to be close to " + convert_to_string(time) 
                   + "  Time difference is > " + convert_to_string(time_eps);
      echo_tee(msg);
      is_consistent = false;
    }
  }

  return is_consistent;
}

//-----------------------------------------------------------------------------

void Problem_Manager::transfer_setup(vector<ID_TYPE>& ids, ID_TYPE id) 
{
  // if this id not already involved in a transfer, add it to the
  // list.  we keep the list sorted and add the id only once, no 
  // matter how many transfer operations it's actually involved in

  vector<ID_TYPE>::iterator end = ids.end();
  vector<ID_TYPE>::iterator it = lower_bound(ids.begin(), ids.end(), id);

  if (it == end || id < *it) {
    ids.insert(it, id);
  }
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Data_Transfer_Operator> 
Problem_Manager::find_pre_target_xfer_operator(ID_TYPE id) const
{
  Teuchos::RCP<Data_Transfer_Operator> oper;
  vector< Teuchos::RCP<Data_Transfer_Operator> >::const_iterator
    it = preelimination_transfers_.begin(), end = preelimination_transfers_.end();
  for (; it != end; ++it) {
    if (id == (*it)->target_id()) {
      oper = *it;
      break;
    }
  }
  return oper;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Data_Transfer_Operator> 
Problem_Manager::find_post_source_xfer_operator(ID_TYPE id
  , vector< Teuchos::RCP<Data_Transfer_Operator> >& pets) const
{
  Teuchos::RCP<Data_Transfer_Operator> oper;
  vector< Teuchos::RCP<Data_Transfer_Operator> >::iterator
    it = pets.begin(), end = pets.end();
  for (; it != end; ++it) {
    if (id == (*it)->source_id()) {
      oper = *it;
      pets.erase(it);
      break;
    }
  }
  return oper;
}

//-----------------------------------------------------------------------------

void Problem_Manager::transfer_details(vector<Teuchos::RCP<Data_Transfer_Operator> >& ops
  , vector<ID_TYPE>& sids, vector<ID_TYPE>& tids) 
{
#if 0
  // ensure lists are sorted or logic is broken
  assert( is_sorted(sids.begin(), sids.end()) );
  assert( is_sorted(tids.begin(), tids.end()) );
#endif

  // diagnose and report circular dependencies

  // re-generate the list of elimination modules
  elimination_modules.clear();
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    const Model_Evaluator *problem_api = dynamic_cast<const Model_Evaluator*>(&problem);
    if (problem_api && !problem_api->supports_standalone_solve() 
      && !inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x)
      && inArgs.modelEvalDescription() != "epistemic_model_evaluator"
      && inArgs.modelEvalDescription() != "aleatory_model_evaluator")
    {
      Teuchos::RCP<Elimination_Module> p = Teuchos::rcp_dynamic_cast<Elimination_Module>(get_problem(problem_ids_[i]));
      if (p != Teuchos::null) elimination_modules.push_back(p);
    }
  }

  // n log(m) search for dependent transfers
  // n is number of unique source module evaluators
  // m is number of unique target module evaluators

  vector<ID_TYPE>::iterator s_it = sids.begin(), s_end = sids.end();
  vector<ID_TYPE>::iterator t_it = tids.begin(), t_end = tids.end();

  bool is_pre = &ops == &preelimination_transfers_ ? true : false;

  for (; s_it != s_end; ++s_it) {
    if (binary_search(t_it, t_end, *s_it)) {

      // we have an elimination module that is 
      // both the source and target of either a 
      // pre or post elimination data transfer

      // get the data transfer operator that we have to reschedule

      ID_TYPE lid = 0;
      Teuchos::RCP<Data_Transfer_Operator> oper;
      vector< Teuchos::RCP<Data_Transfer_Operator> >::iterator it, end;
      it = ops.begin(), end = ops.end();
      for (; it != end; ++it) {
        lid = is_pre ? (*it)->source_id() : (*it)->target_id();
        if (*s_it == lid) {
          oper = *it;
          break;
        }
      }
      assert( !Teuchos::is_null(oper) );

      // get the elimination module that we have to reschedule
     
      lid = is_pre ? oper->source_id() : oper->target_id();
      Teuchos::RCP<Elimination_Module> m = Teuchos::rcp_dynamic_cast<Elimination_Module>(get_problem(lid));

      if (!Teuchos::is_null(m)) 
      {
        // then this transfer operator is associated with an elimination module and it must be sequenced

        // (1) remove transfer from elimination list

        elimination_modules.erase(find(elimination_modules.begin(), elimination_modules.end(), m));

        // (2) remove transfer from pre or post list, as appropriate

        ops.erase(find(ops.begin(), ops.end(), oper));

        // (3) add a perform elimination followed by pre elimination data transfer
        //     or post elimination data transfer followed by a perform elimination

        pair_om pom(oper, m);
        if (is_pre) {
          perform_pres.push_back(pom);
        }
        else {
          post_performs.push_back(pom);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

void Problem_Manager::add_transfer(Teuchos::RCP<Data_Transfer_Operator> & xfer ) 
{ 
  verify_not_register_complete();

  if( 0 != my_hierarchy_ )
    throw runtime_error("Transfers can only be added to the top-level Problem_Manager.");

  if( Teuchos::is_null(xfer->source()) ) // old-style using ID_TYPE construction
  {
    pair<ID_TYPE,ID_TYPE> key(xfer->source_id(), xfer->target_id());
    all_registered_transfers_.insert(pair<pair<ID_TYPE,ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >(key,xfer));
  }
  else // new-style using problems in construction
  {
    pair<ID_TYPE,ID_TYPE> key(xfer->source()->get_id(), xfer->target()->get_id());
    // I wanted to do this ...
    //all_registered_transfers_.insert(pair<pair<ID_TYPE,ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >(key,xfer));
    // ... but ended up having to do this contortion to avoid 
    //     non-deletion warnings from Teuchos::RCP at the end of the run.
    //     I don't know why.  - RWH
    all_registered_transfers_.insert(pair<pair<ID_TYPE,ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >(key,Teuchos::rcp(&(*xfer),false)));
  }

  // This will be populated by a static graph analysis during register_complete
  //standard_transfers_.push_back(xfer); 
}

//-----------------------------------------------------------------------------

void Problem_Manager::add_preelimination_transfer(Teuchos::RCP<Data_Transfer_Operator> xfer ) 
{ 
  preelimination_transfers_.push_back(xfer); 

  ID_TYPE sid = (*xfer).source_id();
  ID_TYPE tid = (*xfer).target_id();

  // track which problem ids are involved in transfers
  transfer_setup(pre_sids, sid);
  transfer_setup(pre_tids, tid);

  // determine if pre-elimination dependencies exist
  transfer_details(preelimination_transfers_, pre_sids, pre_tids);
}

//-----------------------------------------------------------------------------

void Problem_Manager::add_postelimination_transfer(Teuchos::RCP<Data_Transfer_Operator> xfer ) 
{ 
  postelimination_transfers_.push_back(xfer); 

  ID_TYPE sid = (*xfer).source_id();
  ID_TYPE tid = (*xfer).target_id();

  // track which problem ids are involved in transfers
  transfer_setup(post_sids, sid);
  transfer_setup(post_tids, tid);

  // determine if post-elimination dependencies exist
  transfer_details(postelimination_transfers_, post_sids, post_tids);
}

//-----------------------------------------------------------------------------

void Problem_Manager::perform_transfers(vector<Teuchos::RCP<Data_Transfer_Operator> > const & xfers) const
{
  if( residual_active_ )
    copy_composite_solution_into_problem_data();

  vector<Teuchos::RCP<Data_Transfer_Operator> >::const_iterator it = xfers.begin(),
    it_end = xfers.end();
  for (; it != it_end; ++it) 
  {
    bool success = (*it)->perform_data_transfer();
    if( !success )
      throw runtime_error("Transfer operation failed.");
  }
}

//-----------------------------------------------------------------------------

void Problem_Manager::perform_all_transfers() const
{
  if( residual_active_ )
    copy_composite_solution_into_problem_data();

  // First perform all "standard" transfers (which do not involve nonlinear elimination)
  // This is the old transfer system
  perform_transfers(standard_transfers_);
  // ... and this is the new - do both for now.
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    perform_transfers( problem.get_incoming_transfers() ); // do standard transfers needed by this ME
  }

  // Next perform specialized transfers involved in nonlinear elimination
  vector<Teuchos::RCP<Data_Transfer_Operator> >::const_iterator it, end;

  // regular pre elimination transfer step

  it = preelimination_transfers_.begin();
  end = preelimination_transfers_.end();

  if (transfer_type == successive_xfer) {
    // if no complex (cyclic) dependencies
    if (perform_pres.empty() && post_performs.empty()) {

      Teuchos::RCP<Elimination_Module> m; 
      Teuchos::RCP<Data_Transfer_Operator> dto;

      // save off copy of post-elimination transfers
      vector< Teuchos::RCP<Data_Transfer_Operator> > pets(postelimination_transfers_);

      // iterate pre-elimination data transfer targets
      vector<ID_TYPE>::const_iterator id_it = pre_tids.begin(), id_end = pre_tids.end();
      for (; id_it != id_end; ++id_it) {
        // do the pre-elimination transfer
        dto = find_pre_target_xfer_operator(*id_it);
        dto->perform_data_transfer();

        // get the elimination module and invoke it
        m = Teuchos::rcp_dynamic_cast<Elimination_Module>(dto);
        if (!is_null(m)) m->perform_elimination();

        // do the post-elimination transfer, if any
        dto = find_post_source_xfer_operator(*id_it, pets);
        if (!is_null(dto)) dto->perform_data_transfer();
      }

      // iterate post-elimination transfers that didn't have matching pre-transfer
      id_it = post_sids.begin(), id_end = post_sids.end();
      for (; id_it != id_end; ++id_it) {
        // do the post-elimination transfer
        dto = find_post_source_xfer_operator(*id_it, pets);
        if (!is_null(dto)) {
          // do the elimination
          m = Teuchos::rcp_dynamic_cast<Elimination_Module>(get_problem(*id_it));
          m->perform_elimination();
          // post elimination transfer
          dto->perform_data_transfer();
        }
      }
    }
    else {
      assert(0 /* haven't coded up this case yet */);
    }
  }
  else if (transfer_type == fixed_point_xfer) {
    for (; it != end; ++it) {
      (*it)->perform_data_transfer();
    }

    std::vector<pair_om>::const_iterator pp_it = perform_pres.begin(), pp_end = perform_pres.end();
    for (; pp_it != pp_end; ++pp_it) {
      (*pp_it).second->perform_elimination();
      (*pp_it).first->perform_data_transfer();
    }

    // iterate elimination modules performing elimination as we proceed

    vector< Teuchos::RCP<Elimination_Module> >::const_iterator 
      e_it = elimination_modules.begin()
    , e_end = elimination_modules.end();

    for (; e_it != e_end; ++e_it) {
      bool success = (*e_it)->perform_elimination();

    if( !success )
      throw ("Elimination Module failed to solve.");
#if 0
    else {
      // not sure that we want to print every time we successfully perform elimination
      cout << "\n" << "Elimination Module #" << elim_module->first << " succeeded. Final output :" 
        << *(elim_module->second->get_output_vector()) << endl ;
    }
#endif
    }

    pp_it = post_performs.begin();
    pp_end = post_performs.end();

    for (; pp_it != pp_end; ++pp_it) {
      (*pp_it).first->perform_data_transfer();
      (*pp_it).second->perform_elimination();
    }

    // regular post elimination transfer step

    it = postelimination_transfers_.begin();
    end = postelimination_transfers_.end();
    for (; it != end; ++it) {
      (*it)->perform_data_transfer();
    }
  }
}

//-----------------------------------------------------------------------------

void Problem_Manager::setup_fixed_point()
{
  if(problems_.empty())
    throw std::runtime_error("No problems registered with Problem_Manager !!");

  // Configure fixed-point solver using xml file, defauts otherwise
  string convFile = "Problem_Manager_setup_fp.xml";
  Teuchos::RCP<Teuchos::ParameterList> fp_config_list;
  fp_config_list = readParameterList(convFile);
  if( !Teuchos::is_null(fp_config_list) )
  {
    if( fp_config_list->isParameter("Solve Mode") )
    {
      string mode = fp_config_list->get<string>("Solve Mode");
      if( "Jacobi" == mode )
        fp_solve_mode_ = jacobi;
      else if( "Seidel" == mode )
        fp_solve_mode_ = seidel;
      else
        throw runtime_error("Error: Invalid \"Solve Mode\", \"" + mode + "\"; must be \"Jacobi\" or \"Seidel\".");
    }
    if( fp_config_list->isParameter("Maximum Iterations") )
      max_fp_iters_ = fp_config_list->get<int>("Maximum Iterations");
    if( fp_config_list->isParameter("Absolute Tolerance") )
      fp_abs_tolerance_ = fp_config_list->get<double>("Absolute Tolerance");

    if( fp_config_list->isParameter("PreIteration Convergence Check") )
    {
      bool pichk = fp_config_list->get<bool>("PreIteration Convergence Check");
      if(pichk)
        stateful_transfers = false;
      else if(!pichk)
        stateful_transfers = true;
      else
        throw runtime_error("Error in Problem_Manager_setup_fp.xml. \"PreIteration Convergence Check\" must be true of false");
    }
  }

  // Create standalone solvers for each problem that supplies a residual
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    if( problem.supports_standalone_solve() )
    {
      standalone_solvers_[problem_ids_[i]] = problems_[problem_ids_[i]];
    }
#ifdef HAVE_LIME_JFNK
    else if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
    {
      // Use an existing xml file if present, default otherwise
      string nlFile = inArgs.modelEvalDescription() + "_nox.xml";
      Teuchos::RCP<Teuchos::ParameterList> nl_params = readParameterList(nlFile);
      if( Teuchos::is_null(nl_params) )
        nl_params = create_default_nox_params(comm_);

      // Use an existing xml file if present, default otherwise
      string convFile = inArgs.modelEvalDescription() + "_conv.xml";
      Teuchos::RCP<Teuchos::ParameterList> conv_params = readParameterList(convFile);
      if( Teuchos::is_null(conv_params) )
        conv_params = create_default_status_test_list();

      NOX::Utils utils;
      Teuchos::RCP<NOX::StatusTest::Generic> status_test = NOX::StatusTest::buildStatusTests(*conv_params, utils);

      if( Teuchos::is_null(status_test) )
        throw std::runtime_error("Could not create a valid status test for problem \"" + inArgs.modelEvalDescription() + "\"");

      // Create the solver

      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> pinterface = Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(get_problem(problem_ids_[i])));
      Teuchos::RCP<Epetra_Operator> pOp = Teuchos::rcp(new problem_as_operator(get_problem(problem_ids_[i]).getRawPtr()));
      // Configure using an existing xml file if present, default to old mode otherwise
      string setupFile = inArgs.modelEvalDescription() + "_setup.xml";
      Teuchos::RCP<Teuchos::ParameterList> setup_list = readParameterList(setupFile);
      if( !Teuchos::is_null(setup_list) )
      {
        // Allow this problem to supply both preconditioner callback interface and operator
        //   This supports JFNK with Physics-based preconditioning
        setup_list->sublist("Preconditioner").set<NOX::Epetra::Interface::Preconditioner*>("Interface", pinterface.getRawPtr());
        setup_list->sublist("Preconditioner").set<Epetra_Operator*>("Operator", pOp.getRawPtr());
      }
      else // create a default setup list that supports Physics-Based preconditioning
      {
        setup_list = Teuchos::rcp(new Teuchos::ParameterList);
        if( problem.supports_preconditioning() )
        {
          setup_list->sublist("Preconditioner").set<NOX::Epetra::Interface::Preconditioner*>("Interface", pinterface.getRawPtr());
          setup_list->sublist("Preconditioner").set<Epetra_Operator*>("Operator", pOp.getRawPtr());
        }
        else
          setup_list->set("Preconditioner Operator", "None");
      }
      Teuchos::RCP<Epetra_Vector> nonconst_soln = Teuchos::rcp_const_cast<Epetra_Vector>(problem.get_x_init());
      Teuchos::RCP<const NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(nonconst_soln, NOX::Epetra::Vector::CreateView));
      Teuchos::RCP<NOX::Epetra::Group> grpPtr = Teuchos::rcp(create_nox_group(setup_list, nl_params, pinterface, *noxSoln));
      standalone_nox_solvers_[problem_ids_[i]] = NOX::Solver::buildSolver(grpPtr, status_test, nl_params);
      standalone_prob_ids.push_back(problem_ids_[i]);

      // Archive the parameter lists used
      nlFile = pl_dir_name_ + "/" + inArgs.modelEvalDescription() + "_nox_" + create_date_string() + ".xml";
      writeParameterList(nlFile, *nl_params);
      convFile = pl_dir_name_ + "/" + inArgs.modelEvalDescription() + "_conv_" + create_date_string() + ".xml";
      writeParameterList(convFile, *conv_params);
      setupFile = pl_dir_name_ + "/" + inArgs.modelEvalDescription() + "_setup_" + create_date_string() + ".xml";
      if( !Teuchos::is_null(setup_list) )
        writeParameterList(setupFile, *setup_list);
    }
#endif
    else
    {
      // Treat this problem as an elimination module.
      Teuchos::RCP<Elimination_Module> elim_module = Teuchos::rcp_dynamic_cast<Elimination_Module>(get_problem(problem_ids_[i]));;
      if( !Teuchos::is_null(elim_module) )
        elimination_modules_[problem_ids_[i]] = elim_module;
      else
      {
#if 0
        cout << "ERROR: Problems are not all Model_Evaluators, Model_Evaluators or Elimination_Modules." << endl;
        throw "Problem_Manager ERROR";
#else
        // Noel removed fatal error as this module is time independent, not a solver, elimination or evaluator model
#endif
      }
    }
  }

  // Sanity check for the algorithm
  if( (numeric_limits<double>::min() >= fp_abs_tolerance_) &&
      (numeric_limits<int>::max() <= max_fp_iters_ ) )
    throw std::runtime_error(
               "\n*******************************************************\n"
     + (string)"ERROR: Fixed-point has no valid stopping criteria.\n"
     + (string)"       Make sure you have Problem_Manager_setup_fp.xml\n"
     + (string)"       in the run directory.\n"
     + (string)"*******************************************************");
}

//-----------------------------------------------------------------------------
// This LIME routine solves the coupled problem set using a fixed point algorythm
//-----------------------------------------------------------------------------
//
bool Problem_Manager::solve_fixed_point(double & dt)
{
  bool converged = false;

  double original_dt = dt;

  // The pre-iteration convergence check that comes next can be turned off by adding the 
  // following line to the LIME input file "Problem_Manager_setup_fp.xml"
  //   <Parameter name="PreIteration Convergence Check" type="bool" value="false"/>

  if (!stateful_transfers) {
    // Start all problems with synchronized states from transfers
    perform_all_transfers();
    // ... and establish initial convergence state
    converged = compute_fixed_point_convergence();
  }

  int fixed_pt_iters = 0; 
  // Top of the Fixed Point Iteration Loop   
  while( !converged && (fixed_pt_iters < max_fp_iters_) )
  {

    // NOTE: Parts A and B below are mutually exclusive for a given ME
    // TODO : Also, this hard-coded algorithm structure will break a seidel mode
    // whenever a NOX-based ME sove is intended to precede a standalone solve.
    // 
    
    // ------------------------------------------------------------------------
    // A. Set up c++ stl iterator to loop over each of the registered codes 
    //    supporting the stand-alone solve interface
    for( size_t i = 0; i < problem_ids_.size(); ++i )
      if( standalone_solvers_.count(problem_ids_[i]) )
      {
        map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::iterator 
          sas_iter = standalone_solvers_.find(problem_ids_[i]);

        // The following allows us to reference the specific problem we are solving (instead 
        // of a pointer to it. 
        Model_Evaluator & problem  = *(sas_iter->second);

        // TODO: Transfers should have handled this. Consider removing    
        copy_composite_solution_into_problem_data(); // this may need to be reexamined

        // If the solver mode is seidel, do required transfers for next stand-alone solve
        if( seidel == fp_solve_mode_ )
          perform_transfers( problem.get_incoming_transfers() );

        // Perform the solve_standalone for the current physics problem
        problem.solve_standalone(dt);

        // The first time a problem changes the time step during its solve, we 
        // return with converged = false.
        if( has_time_step_changed(original_dt, dt) )
          return false;

        // TODO: This should be reconsidered and perhaps handled with appropriate transfers    
        copy_problem_into_composite( sas_iter->first, SOLUTION_X );
      }


    // ------------------------------------------------------------------------
    // B. Set up c++ stl iterator to loop over each of the registered codes 
    //    supporting the Trilinos NOX stand-alone solver interface
#ifdef HAVE_LIME_JFNK
    map<ID_TYPE, Teuchos::RCP<NOX::Solver::Generic> >::iterator 
      nox_iter     = standalone_nox_solvers_.begin(),
      nox_iter_end = standalone_nox_solvers_.end();

    // Now perform the loop
    for( ; nox_iter_end != nox_iter; ++nox_iter )
    {

      // The following allows us to reference the specific problem we are solving (instead 
      // of a pointer to it. 
      Model_Evaluator & problem  = *get_problem(nox_iter->first);

      // get the companion NOX solver for the ME
      Teuchos::RCP<NOX::Solver::Generic> nox_solver = nox_iter->second;

      // TODO: Replace this with an appropriate transfer
      copy_composite_solution_into_problem_data();
      
      // The following steps are equivalent to a solve_standalone but now using a NOX
      // solver with the ME
      {
         // This contortion is needed to reset the NOX solver with updated solution values
         Teuchos::RCP<NOX::Epetra::Vector> noxVec = Teuchos::rcp(new NOX::Epetra::Vector(problem.get_x_state(), NOX::Epetra::Vector::CreateView));
         nox_solver->reset(*noxVec);

         NOX::StatusTest::StatusType status = nox_solver->solve();

         // Get the Epetra_Vector with the final solution from the solver
         const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(nox_solver->getSolutionGroup());
         const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

         string msg = "";
         if( NOX::StatusTest::Converged == status )
           msg = "\n-- NOX solver converged --";
         else
           msg = "\n-- NOX solver did not converge --";
         echo_tee(msg);
       
         // TODO: This should be reconsidered and perhaps handled with appropriate transfers    
         copy_vector_into_composite_solution( nox_iter->first, finalSolution );
      }
    }
#endif

    // ------------------------------------------------------------------------
    // C. Do a convergence check and update iteration counter
    
    //cout << "\tFor fixed-point iteration # " << fixed_pt_iters << endl;
    converged = compute_fixed_point_convergence();
    //cout << "\t\tConverged --> " << converged << endl;

    ++fixed_pt_iters;
    
  } // Bottom of the Fixed Point Iteration Loop

// Output result of fixed-point iteration  
#if 1
  string msg = "\tFixed-Point iteration " + (string)(converged?"CONVERGED":"DID NOT CONVERGE") 
    + " after " + convert_to_string(fixed_pt_iters) + " iterations.";
  echo_tee(msg);
#endif

  // Note: This routine returns a boolean that is true only if things converged
  return converged;
}

//-----------------------------------------------------------------------------

bool Problem_Manager::compute_fixed_point_convergence()
{
  bool converged = false;

  // Here we refresh data used by all modules in the fixed-point iteration and then have
  // those that can compute its residual.  We then construct a norm from contributions from each.

  if (residual_active_) 
  {
#ifdef HAVE_LIME_JFNK
    bool success = compute_composite_residual();

    // set flag computed_residual to true
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      problem.computed_residual(true);
    }

    if( !success )
      throw std::runtime_error("Could not compute composite residual for use in compute_composite_residual.");

    double normL2 = numeric_limits<double>::max();
    int ierr = composite_residual_->Norm2(&normL2);
    if( ierr )
      throw std::runtime_error("Could not compute L2-norm of composite residual for use in compute_composite_residual.");

    if( fp_abs_tolerance_ > normL2 )
      converged = true;

    if (residual_active_)
      echo_tee("||F|| = "+convert_to_string(normL2)+" - Fixed-point");
    if( converged )
      echo_tee(" - Converged");
//    else
//      echo_tee("");    // adds an extra line after printing current value of composite residual
#endif
  }
  else
  {
    converged = true;

    perform_all_transfers();

    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      if( problem.supports_standalone_solve() ) {
        bool result = problem.is_converged();
        converged = converged && result;
      }
    }
  }

  return converged;
}

//-----------------------------------------------------------------------------

#ifdef HAVE_LIME_JFNK

void Problem_Manager::setup_jfnk()
{
  if(problems_.empty())
    throw runtime_error("No problems registered with Problem_Manager !!");

  // Create a coupled problem NOX solver with us as the callack interface

  // Use default parameters and allow overwrite if xml file is present
  Teuchos::RCP<Teuchos::ParameterList> nl_params = create_default_nox_params(comm_);
  string nlFile = "Problem_Manager_setup_jfnk.xml";
  Teuchos::RCP<Teuchos::ParameterList> mod_nl_params = readParameterList(nlFile);
  if( ! Teuchos::is_null(mod_nl_params) )
    nl_params->setParameters(*mod_nl_params);

  // Redirect detailed solver output to a log file.  The filename could be made user-input.
  nl_params->sublist("Printing").set("Output Stream", out_stream_);
  nl_params->sublist("Printing").set("Error Stream", out_stream_);

  // Use an existing xml file if present, default otherwise
  Teuchos::RCP<Teuchos::ParameterList> status_test_list_ = create_default_status_test_list();
  string convFile = "Problem_Manager_conv.xml";
  Teuchos::RCP<Teuchos::ParameterList> mod_conv_params = readParameterList(convFile);
  if( ! Teuchos::is_null(mod_conv_params) )
    status_test_list_->setParameters(*mod_conv_params);

  // Support a few simple options that are specified per the User Manual andconsistent with Fixed-Point
  // These will take precedence over corresponding options set in the status tests xml file
  if( nl_params->isParameter("Absolute Tolerance") )
  {
    Teuchos::ParameterList& normF_abs = status_test_list_->sublist("Test 0").sublist("Test 0");
    normF_abs.set("Tolerance", nl_params->get<double>("Absolute Tolerance"));
  }
  if( nl_params->isParameter("Maximum Iterations") )
  {
    Teuchos::ParameterList& maxiters = status_test_list_->sublist("Test 2");
    maxiters.set("Maximum Iterations", nl_params->get<int>("Maximum Iterations"));
  }

  NOX::Utils utils;
  status_tests = NOX::StatusTest::buildStatusTests(*status_test_list_, utils);


  // Add a user defined pre/post operator object if verbosity is on
  if( m_verbose )
  {
    Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = Teuchos::rcp(new PrePostOperator());
    nl_params->sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);
  }

  Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(*composite_solution_));

  // Configure our coupling solver using an existing xml file if present, default to old mode otherwise
  string setupFile = m_setupXMLfilename;
  if( "" == setupFile )
    setupFile = "Problem_Manager_setup.xml";
  Teuchos::RCP<Teuchos::ParameterList> setup_list = readParameterList(setupFile);
  if( Teuchos::is_null(setup_list) )  // create a default setup list that supports Physics-Based preconditioning
    setup_list = Teuchos::rcp(new Teuchos::ParameterList);

  // Make ourself (Problem_Manager) available as the preconditioner callback interface and operator
  //   This supports JFNK with Physics-based preconditioning
  Teuchos::RCP<EpetraExt::ModelEvaluator> prob_manager_me = Teuchos::rcp(this, false);
  my_nox_interface_ = Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(prob_manager_me));
  //setup_list->sublist("Preconditioner").set<NOX::Epetra::Interface::Preconditioner*>("Interface", my_nox_interface_.getRawPtr());
  //setup_list->sublist("Preconditioner").set<Epetra_Operator*>("Operator", this);
  //setup_list->set("Preconditioner Operator", "None");

  grp_ = Teuchos::rcp(create_nox_group(setup_list, nl_params, my_nox_interface_, *noxSoln));
  if( Teuchos::is_null(grp_) )
    throw runtime_error("Failed to create a valid NOX::Epetra::Group in Problem_Manager::setup_jfnk.");

  // Create the solver
  composite_nox_solver_ = NOX::Solver::buildSolver(grp_, status_tests, nl_params);

  // Archive the parameter lists used
  nlFile = pl_dir_name_ + "/Problem_Manager_nox_" + create_date_string() + ".xml";
  writeParameterList(nlFile, *nl_params);
  convFile = pl_dir_name_ + "/Problem_Manager_conv_" + create_date_string() + ".xml";
  writeParameterList(convFile, *status_test_list_);
  setupFile = pl_dir_name_ + "/Problem_Manager_setup_" + create_date_string() + ".xml";
  if( !Teuchos::is_null(setup_list) )
    writeParameterList(setupFile, *setup_list);

}

//-----------------------------------------------------------------------------

bool Problem_Manager::solve_jfnk(double & dt)
{
  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;
  status = composite_nox_solver_->solve();

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(composite_nox_solver_->getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  *composite_solution_ = finalSolution; 

  bool converged = true;
  if( NOX::StatusTest::Converged == status )
  {
    int num_iters = composite_nox_solver_->getList().sublist("Output").get<int>("Nonlinear Iterations");
    double finalL2norm = composite_nox_solver_->getList().sublist("Output").get<double>("2-Norm of Residual");
    echo_tee("-- NOX solver converged in "+convert_to_string(num_iters)+" to final L2-norm = "+convert_to_string(finalL2norm)+"  --\n");
  }
  else 
  {
    converged = false;
    echo_tee("-- NOX solver did not converge --\n");
  }

  //cout << "\n" << "Final JFNK-coupled solution :\n" << std::setprecision(16) << finalSolution << endl;

  // And output nonlinear elimination solutions
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    if( !inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
    {
      Teuchos::RCP<Elimination_Module> elim_module = Teuchos::rcp_dynamic_cast<Elimination_Module>(get_problem(problem_ids_[i]));;
#if 0
      if( !Teuchos::is_null(elim_module) )
        cout << "\n" << "Elimination Module #" << i << " final output :\n" 
          << *(elim_module->get_output_vector()) << endl ;
#endif
    }
  }

  // This is a one-off for SIAM CSE studies.
  // I want to determine if/when fixed-point will fail to converge by computing the
  // matrix norm of the fixed-point operator.  For the con1d problem,
  // R = G(T) - T, so that J = dG/dT - I.  I want to compute dG/dT numerically 
  // using a FD operator and then find its matrix norm.  I hope to find some parameters
  // and values that will make or break fixed-point convergence.

  if( 0 )
  {
    Teuchos::RCP<Teuchos::ParameterList> nl_params = create_default_nox_params(comm_);
    Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(composite_solution_, NOX::Epetra::Vector::CreateView));
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = my_nox_interface_;
    Teuchos::RCP<NOX::Epetra::FiniteDifference> myFD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(nl_params->sublist("Printing"), iReq, *noxSoln));

    bool success = myFD->computeJacobian(finalSolution);
    if( !success )
      throw std::runtime_error("Problem_Manager :  FiniteDifference failed to computeJacobian.");

    echo_tee("FD Jacobian using :\n"+convert_to_string(finalSolution)+"\n has inf-norm = "+convert_to_string(myFD->NormInf()));

    cout << "and here's the matrix:\n" << myFD->getUnderlyingMatrix() << endl;

    grp_->setX(finalSolution);
    //composite_nox_solver_->reset(grp_->getX(), status_tests);
    grp_->computeJacobian();
    Teuchos::RCP<const Epetra_Operator> myJac = grp_->getLinearSystem()->getJacobianOperator();
    dump_operator_to_file(*myJac, "PM_Jacobian");
  }

  return converged;
}

//-----------------------------------------------------------------------------

void Problem_Manager::setup_switching()
{
  setup_fixed_point();
  setup_jfnk();

  // Configure switching solver
  string setupFile = m_setupXMLfilename;
  if( "" == setupFile )
    setupFile = "Problem_Manager_setup.xml";
  Teuchos::RCP<Teuchos::ParameterList> setup_list = readParameterList(setupFile);
  if( Teuchos::is_null(setup_list) )
    throw runtime_error("Error: No valid Problem_Manager_setup.xml file needed to configure switching solver.");

  if( !setup_list->isSublist("Switching Solver") )
    throw runtime_error("Error: No valid \"Switching Solver\" sublist in Problem_Manager_setup.xml file.");

  Teuchos::ParameterList & ssolverParams = setup_list->sublist("Switching Solver");
  if( !ssolverParams.isParameter("Number of Stages") )
    throw runtime_error("Error: No specified \"Number of Stages\" in \"Switching Solver\" sublist.");
  int nstages = ssolverParams.get<int>("Number of Stages");

  // This populates the parameter lists in order of stage ids
  for( int i = 0; i < nstages; ++i )
  {
    std::ostringstream int_name;
    int_name << i << std::flush;

    string sublistName = "Stage " + int_name.str();
    if( !ssolverParams.isSublist(sublistName) )
      throw runtime_error("Error: No valid \"" + sublistName + "\" sublist in \"Switching Solver\" sublist in Problem_Manager_setup.xml file.");
    switching_solve_params_.push_back(Teuchos::rcp( new Teuchos::ParameterList(ssolverParams.sublist(sublistName)) ) ); 
  }

  return;
}

//-----------------------------------------------------------------------------

bool Problem_Manager::solve_switching(double & dt)
{
  bool converged = true;

  for( unsigned i = 0; i < switching_solve_params_.size(); ++i )
  {
    Teuchos::ParameterList & solverParams = *switching_solve_params_[i];
    if( !solverParams.isParameter("Solver Type") )
      throw runtime_error("Error: No specified \"Solver Type\" in \"Switching Solver\" sublist.");
    string sType = solverParams.get<string>("Solver Type");

    solve_using su;
    if( "Fixed-Point" == sType )
      su = fixed_point;
    else if( "JFNK" == sType )
      su = jfnk;
    else
      throw runtime_error("Error: Invalid \"Solver Type\", \"" + sType + "\"; must be \"Fixed-Point\" or \"JFNK\".");

    if( fixed_point == su )
    {
      compute_fixed_point_convergence(); // Synchronizes data
      converged = solve_fixed_point(dt) && converged;
      grp_->setX(*composite_solution_);
      copy_composite_solution_into_problem_data();
      grp_->computeF();
    }
    else
    {
      converged = solve_jfnk(dt) && converged;
    }
  }

  return converged;
}

//-----------------------------------------------------------------------------

void Problem_Manager::apply_predictors(solve_using su)
{
  if( su == jfnk )
  {
    if( Teuchos::is_null(grp_) )
      throw runtime_error("Error: No valid NOX::Epetra::Group; one is needed for doing a predictor for JFNK.");

    // Compute overall residual norm before any improvements via predictors.
    // This will be used to compute a NOX residual below.

    grp_->computeF();
    double initial_norm = grp_->getNormF();
    (*out_stream_) << "Predictor L2-norm before any predictors (composite residual) = " << convert_to_string(initial_norm);

    // pre step initialize
    bool predictor = false;
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      // if model evaluator has an initial guess for state vector
      if( problem.supports_predictor() )
      {
        predictor = true;
        problem.solve();
        copy_vector_into_composite_solution(problem_ids_[i], *problem.get_x_init());
      }
    }

    // if someone had an initial guess, compute global residual
    double updated_norm = initial_norm;
    if (predictor) 
    {
      grp_->setX(*composite_solution_);
      copy_composite_solution_into_problem_data();
      composite_nox_solver_->reset(grp_->getX(), status_tests);
      grp_->computeF();
      updated_norm = grp_->getNormF();
      echo_tee("Predictor L2-norm (composite residual) = "+convert_to_string(updated_norm));
    }
  }
  else
  {
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      // if model evaluator has an initial guess for state vector
      if( !problem.supports_standalone_solve() && problem.supports_predictor() )
      {
        problem.solve();
        // Now reset the appropriate solver with the new solution state - for now involves a vector copy
        Teuchos::RCP<NOX::Epetra::Vector> noxVec = Teuchos::rcp(new NOX::Epetra::Vector(*problem.get_x_init()) );
        standalone_nox_solvers_[problem_ids_[i]]->reset(*noxVec);
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------

void Problem_Manager::update_right_scaling_vectors()
{
  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);

    if( problem.supports_right_scaling() )
    {
      problem.compute_right_scaling_vector();
      copy_vector_into_composite_solution(problem_ids_[i], *problem.get_x_init());
    }
    if( problem.supports_left_scaling() )
      problem.compute_left_scaling_vector();
  }

  return;
}
#endif

//-----------------------------------------------------------------------------

void Problem_Manager::const_integrate() const
{
  return const_cast<Problem_Manager*>(this)->integrate();
} 


//---------------------------------------------------------------------------------------------
// This is the highest-level LIME routine that controls the multi-physics integration
// Both transient and steady-state situations are treated. 
//---------------------------------------------------------------------------------------------
//
void Problem_Manager::integrate()
{

  // ==========================================================================================
  // Part 1: Do all set-up tasks required before doing actual integration

  // ---------------------------------------
  // 1.1 Msc. Initialization
  
  // Initialize some variables to be used for time and time step related values:
  unsigned n = 0;
  unsigned max_steps = std::min(max_time_steps_, numeric_limits<unsigned int>::max());
  double t_max = numeric_limits<double>::max();
  double current_time = -1.0;
  double time_step = numeric_limits<double>::max();

  // ---------------------------------------
  // 1.2  Loop over each physics code to determine time and time step related values:
  //      t_max:          Time when time-integration will end if reached
  //      time_step       Time step size
  //      max_time_steps  Maximum number of time steps that will be taken
  //      current_time    The current value of time (must be the same for all codes)
  //

  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);

     // perform "reinitialze"
     problem.reinitialize();  // added to support uq sampling

     // Skip to next physics code if this one does not provide temporal support
     if( !problem.is_transient() ) continue;

     // Get this physics code's time and step size related values
     current_time = problem.get_current_time();
     double t = problem.get_max_time();
     n = problem.get_max_steps();

     // Warning if n is zero. This may indicate the routine get_max_steps was not implemented properly
     if (n == 0) {
       cout << "Problem_Manager WARNING: Value returned from problem.get_max_steps is 0 " << endl;
       cout << "                         Check implementation of get_time_step " << "\n" << endl;
     }

     // Reconciliation of time-related quantities
     // Overall t_max, time_step size, and max_time_steps always taken as mimimum 
     t_max = std::min(t_max, t);
     max_steps = std::min(max_steps, n);
  }

  // ---------------------------------------
  // 1.3 Now make adjustments to time and time step related values as needed to fit 
  //     particular cases and be sure no inconsistencies exist. 
  

  // TODO: Place holder, ADD CODE HERE
  // assert that dt and t_max are positive
  // dt is less than t_max

  // If overall problem is steady-state, setup time loop for single pass
  if( !problem_is_transient )
  {
     t_max = 0.0;
     max_steps = 1;
     current_time = 0.0;
  }

 
  // ==========================================================================================
  // Part 2: Perform time integration (or steady-state solve if problem isn't a transient)
  //         Loop will continue till max_time_steps or t_max are exceeded

  n = 1;
  // Top of Time integration loop 
  // while ( n <= max_steps && current_time < t_max  )  // RCS would like this change, but it breaks the tests
  while ( n <= max_steps && current_time <= t_max  )
  {
     // ---------------------------------------
     // 2.1  A set of checks and pre-solve preparation before integrating one step in time
  
     // Ensure consistent current time among models
     if( !verify_consistent_current_time(current_time) )
       echo_tee("Warning: Current solver times are not consistent among all applications.");

     // Update all physics packages with agreed upon timestep
     time_step = determine_time_step();

     for( size_t i = 0; i < problem_ids_.size(); ++i )
       if( get_problem(problem_ids_[i])->is_transient() )
         get_problem(problem_ids_[i])->set_time_step(time_step);

#ifdef HAVE_LIME_JFNK
     // if the problem supports a predictor step, call the "apply_predictors" routine
     if( use_predictor_ ) apply_predictors(su_);

     // This routine will perform RHS scaling if it is supported
     update_right_scaling_vectors();

     // This is needed for when you are doing a matrix-Free solution.  In this case
     // the NOX solver is reset to reflect any changes to the composite solution state.
     if( !Teuchos::is_null(grp_) )
     {
       grp_->setX(*composite_solution_);
       composite_nox_solver_->reset(grp_->getX(), status_tests);
     }
#endif

     // One-line output for transient problems before doing non-linear integration
     if (problem_is_transient) {
       string msg = "\nDoing nonlinear solve for step " + convert_to_string(n) +
         ", max_time_step = " + convert_to_string(max_steps) + 
         ", start time = " + convert_to_string(current_time) + ", max_time = " + convert_to_string(t_max) + "\n";
       echo_tee(msg);
     }

     // ---------------------------------------
     // 2.2  Perform integration in time for 1 time step, advance time, and 
     // register convergence status

     // Take one step by calling solve_nonlinear
     double test_time_step = time_step;
     bool converged = solve_nonlinear(test_time_step);
     bool dt_has_changed = has_time_step_changed(test_time_step, time_step);
     if( dt_has_changed )
     {
       if( time_step < test_time_step)
         throw std::runtime_error("LIME does not currently allow dt to increase during a solve.");
       else
         continue; // return to top of time stepping and allow renegitiation of time step size
     }

     // Ascertain overall convergence status of coupled problem for use in informing each problem below. 
     // We could throw here.
     Model_Evaluator::convergence status = Model_Evaluator::converged;
     if (!converged) status = Model_Evaluator::not_converged;

     // Advance the current time using current time step
     current_time += time_step;

     // ---------------------------------------
     // 2.3  Collection of msc. post single time-step integration tasks

     //// Loop over each physics code to allow for any post time step processing:

     for( size_t i = 0; i < problem_ids_.size(); ++i )
     {
       Model_Evaluator & problem  = *get_problem(problem_ids_[i]);

       // notify model evaluators of overall convergence status
       problem.step_converged(status);

       // pass solution vector into model evaluator if problem has active variables
       EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
       if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
         problem.finalize_solution( *(create_view_into_composite_vector(problem_ids_[i], SOLUTION_X)) );

       if( problem.is_transient() )
       {
         // update (advance) current time
         problem.update_time();
       }
     }

#ifdef HAVE_LIME_JFNK
     // Reset composite JFNK solver if warranted
     if (su_ == jfnk)
     {
       // reset the convergence criteria so we can run another step
       const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(composite_nox_solver_->getSolutionGroup());
       const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

       // Problem->reset(finalSolution);
       grp_->setX(finalSolution);
       composite_nox_solver_->reset(grp_->getX(), status_tests);
       grp_->computeF();
     }
#endif

     // update step count
     ++n;

  }
  // Bottom of Time interation loop

  // ==========================================================================================
  // Part 3: Final Post integration details
  
  // One-line output for transient problems at end of time loop
  if (problem_is_transient) {
     string msg = "\nTime step loop completed: time = " + convert_to_string(current_time) + "\n";
     echo_tee(msg);
  }

#ifdef HAVE_LIME_JFNK
  // post solve information dump if you are doing jfnk
  if (su_ == jfnk)
  {
     string nlFile = pl_dir_name_ + "/Problem_Manager_nox_final_" + create_date_string() + ".xml";
     writeParameterList(nlFile, composite_nox_solver_->getList());
  }
#endif
}

//---------------------------------------------------------------------------------------------

void Problem_Manager::setup_nonlinear_solve()
{
  switch( su_ )
  {
    case jfnk : {
#ifdef HAVE_LIME_JFNK
      if (!residual_active_)
        throw runtime_error("Error: Solver Strategy JFNK not valid if no physics can compute a residual");
      setup_jfnk();
#else
      throw std::runtime_error("Error: Problem_Manager not configured to support JFNK solves.");
#endif
      break;
    }

    case fixed_point :  
      setup_fixed_point();
      break;

    case switching :
#ifdef HAVE_LIME_JFNK
      setup_switching();
#else
      throw std::runtime_error("Error: Problem_Manager not configured to support JFNK solves.");
#endif
      break;

    default:
      throw std::runtime_error("Error: Problem_Manager configured with unknown coupling solver type.");
  }
}

//-----------------------------------------------------------------------------

bool Problem_Manager::solve_nonlinear(double & dt)
{
  if( !register_complete_called_ )
    throw std::runtime_error("Error: Problem_Manager solve_nonlinear called before calling register_complete called.");

  bool converged = false;

  switch( su_ )
  {
    case jfnk :  
#ifdef HAVE_LIME_JFNK
      converged = solve_jfnk(dt);
#else
      throw std::runtime_error("Error: Problem_Manager not configured to support JFNK solves.");
#endif
      break;

    case fixed_point :  
      converged = solve_fixed_point(dt);
      break;

    case switching :
#ifdef HAVE_LIME_JFNK
      converged = solve_switching(dt);
#else
      throw std::runtime_error("Error: Problem_Manager not configured to support JFNK solves.");
#endif
      break;

    default:
      throw std::runtime_error("Error: Problem_Manager configured with unknown coupling solver type.");
  }

  return converged;
}

//-----------------------------------------------------------------------------

#ifdef HAVE_LIME_JFNK

EpetraExt::ModelEvaluator::InArgs
Problem_Manager::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription("Problem Manager");

  inArgs.setSupports(IN_ARG_x, true);

  if( !Teuchos::is_null(m_model_p_names) )
    inArgs.set_Np(1); // We currently support only a single parameter vector

  return inArgs;
}

//-----------------------------------------------------------------------------

EpetraExt::ModelEvaluator::OutArgs 
Problem_Manager::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription("Problem Manager");

  // For now, just support residual fills
  outArgs.setSupports(OUT_ARG_f, true);

  if( !m_responses.empty() )
  {
    int np = 0;
    if( !Teuchos::is_null(m_model_p_names) )
      np = 1;  // single params vector for now

    outArgs.set_Np_Ng(np,m_responses.size());
  }

  return outArgs;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Map> 
Problem_Manager::get_g_map(int l) const 
{
  if( 0 != l )
    throw "Problem_Manager::get_g_map : Only single response function (index = 0) exists (for now).";

  const Teuchos::RCP<Model_Evaluator> & me = (m_responses.begin()->second).first;
  //Teuchos::RCP<const Epetra_Map> map = Teuchos::rcp((m_responses.begin()->second)->first->get_g_map(l));
  Teuchos::RCP<const Epetra_Map> map = me->get_g_map(l);

  return map;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Vector> 
Problem_Manager::get_p_init(int l) const 
{
  if( 0 != l )
    throw "Problem_Manager::get_p_init : Only single parameter vector (index = 0) exists.";

  return m_model_p_vec;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Teuchos::Array<std::string> >
Problem_Manager::get_p_names(int l) const
{
  if( 0 != l )
    throw "Problem_Manager::get_p_names : Only single parameter vector (index = 0) exists.";

  return m_model_p_names;
}

//-----------------------------------------------------------------------------

void 
Problem_Manager::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  if( outArgs.get_f().get() ) // Signals a residual fill request, eg computeF from NOX
  {
    const Epetra_Vector x(*(inArgs.get_x().get()));
    Epetra_Vector & f = *(outArgs.get_f().get());
    // compute nonlinear residual
    computeF(x, f, NOX::Epetra::Interface::Required::Residual);
  }
  //if( !Teuchos::is_null(inArgs.get_p(0)) ) // Signals the need to set and propagate problem parameters
  if( inArgs.Np() && !Teuchos::is_null(inArgs.get_p(0)) ) // Signals the need to set and propagate problem parameters
  {
    publish_problem_parameters(inArgs.get_p(0));
    for( size_t i = 0; i < problem_ids_.size(); ++i )
      get_problem(problem_ids_[i])->initializeSolution();

    refresh_solution_state();
  }

  if( !m_responses.empty() )
  {
    if( outArgs.get_g(0).get() ) // Signals a response function fill request
    {
      // Fire off the transient solve
      const_integrate();
      compute_problem_responses();
      const map<string, std::pair<Teuchos::RCP<Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > > & g_vals = get_responses();
      if( 1 < g_vals.size() )
        throw std::logic_error("Problem_Manager currently supports only a single problem response function.");

      const Teuchos::RCP<Epetra_Vector> & gvec = (g_vals.begin()->second).second;
      (*outArgs.get_g(0))[0] = (*gvec)[0];
      echo_tee(" Response Function Value = "+convert_to_string((*gvec)[0]));
    }
  }
  else if( outArgs.get_WPrec().get() ) // Signals a computePreconditioner
  {
    const Epetra_Vector x(*(inArgs.get_x().get()));
    Epetra_Operator & prec = *(outArgs.get_WPrec().get());
    Teuchos::ParameterList * nullList = NULL;
    computePreconditioner(x, prec, nullList);
  }
  else if( outArgs.get_W().get() ) // Signals a computeJacobian
  {
    throw "Problem_Manager::evalModel : Jacobian matrix support not available.";
  }
}

//-----------------------------------------------------------------------------

bool 
Problem_Manager::computeF(const Epetra_Vector& x, Epetra_Vector& F, const NOX::Epetra::Interface::Required::FillType fillFlag) const
{
  string context = "";
  switch( fillFlag )
  {
    case NOX::Epetra::Interface::Required::Residual :   context = "Residual";           break;
    case NOX::Epetra::Interface::Required::Jac :        context = "Jacobian";           break;
    case NOX::Epetra::Interface::Required::Prec:        context = "Preconditioner";     break; 
    case NOX::Epetra::Interface::Required::FD_Res:      context = "Finite-Diff-Res";    break;
    case NOX::Epetra::Interface::Required::MF_Res:      context = "Matrix-Free-Res";    break;
    case NOX::Epetra::Interface::Required::MF_Jac:      context = "Matrix-Free-Jac";    break;
    default:                                            context = "Unknown";            break;
  }

#if 0
  //if( m_verbose )
    cout << std::setprecision(16) << "\nProblem_Manager::computeF (" << context << "): solution : " << x << endl;
#endif

  *composite_solution_ = x;

  bool success = compute_composite_residual();

  F = *composite_residual_;

#if 0
  //if( m_verbose )
    cout << std::setprecision(16) << "\nProblem_Manager::computeF (" << context << "): residual : " << F << endl;
#endif

  if(fillFlag == NOX::Epetra::Interface::Required::Residual) 
  {
    // Set flag computed_residual to true 
    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      const Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      problem.computed_residual(true);
    }
  }

  return success;
}

//-----------------------------------------------------------------------------

bool 
Problem_Manager::compute_composite_residual() const
{
  // don't perform transfer or compute residual if not supported
  if (!residual_active_) return true;

  perform_all_transfers();

  for( size_t i = 0; i < problem_ids_.size(); ++i )
  {
    Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
    EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
    if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
    {
      // Create index mapping for this problem into the composite problem
      inArgs.set_x(problem.get_x_state());
      const ModelEvaluator::Evaluation<Epetra_Vector> f_vec(Teuchos::rcp(&(*problem_residual_data_[problem_ids_[i]]), false));
      outArgs.set_f(f_vec);
      problem.evalModel(inArgs, outArgs);

      copy_problem_into_composite(problem_ids_[i], RESIDUAL);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------

bool 
Problem_Manager::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams) const
{
  if( 0 ) // this is a simple test to exercise brute-force FD of the overall problem
  {
    // Make sure we are still dealing with ourself
    if( this != &M )
      throw std::runtime_error("Address of Epetra_Operator for computePreconditioner is not same as Problem_Manager.");

    if( Teuchos::is_null(FD) )
    {
      Teuchos::RCP<Teuchos::ParameterList> nl_params = create_default_nox_params(comm_);
      Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(composite_solution_, NOX::Epetra::Vector::CreateView));
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = my_nox_interface_;
      FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(nl_params->sublist("Printing"), iReq, *noxSoln));

      // Optional and appears to affect the quality of a direct solve via Amesos
      FD->setDifferenceMethod(NOX::Epetra::FiniteDifference::Centered);

    }

    bool ok = FD->computeJacobian(x, *FD);

    // This is a major hack to see if Rio_Rio is producing a ill-conditioned matrix.
    // This test simply adds to the diagonal to attempt to produce a valid Jacobian... RWH
    if( 0 )
    {
      cout << "Problem_Manager::computePreconditioner : doing the diagonal shift ...." << endl;
      Epetra_Vector diagonal_vector(*composite_solution_);
      Epetra_Vector shift_vector(*composite_solution_);
      shift_vector.PutScalar(1.0);
      Epetra_CrsMatrix & matrix = FD->getUnderlyingMatrix();
      FD->ExtractDiagonalCopy(diagonal_vector);
      diagonal_vector.Update(1.0, shift_vector, 1.0);
      matrix.ReplaceDiagonalValues(diagonal_vector);
    }

    if( 1 )
    {
      dump_operator_to_file(*FD, "FD-precon");
    }

    return ok;
  }
  else // this is the real interface we want to support
  {
    *composite_solution_ = x;

    copy_composite_solution_into_problem_data();

    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
      EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
      if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
      {
        if( problem.supports_preconditioning() )
          problem.compute_preconditioner(*problem.get_x_state());
      }
    }

    if( 0 )
    {
      dump_operator_to_file(*this, "Rio-precon");
    }

    return true;
  }
}

//-----------------------------------------------------------------------------

int 
Problem_Manager::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if( Teuchos::is_null(FD) )
  {
    cout << "Problem_Manager::Apply. Should't happen." << endl; 
    return -1; /* We don't support this. */
  }

  return FD->Apply(X, Y);
}

//-----------------------------------------------------------------------------

int 
Problem_Manager::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
#if 0
  cout << std::setprecision(16) << "\nProblem_Manager::ApplyInverse (apply_preconditioner): incoming X : " << *X(0) << endl;
  cout << std::setprecision(16) << "\nProblem_Manager::ApplyInverse (apply_preconditioner): pre-apply Y : " << *Y(0) << endl;
#endif

  if( 0 ) // this is a simple test to exercise brute-force FD of the overall problem
  {
    if( Teuchos::is_null(FD) )
      throw std::logic_error("Our FiniteDifference Operator is NULL.");

    // We choose to do physics-based preconditioning as a direct solve of a numerically computed Jacobian
    Epetra_RowMatrix * rowMatrix = dynamic_cast<Epetra_RowMatrix *>(&(*FD));
    Epetra_MultiVector * nonconst_X = const_cast<Epetra_MultiVector *>(&X);
    Epetra_LinearProblem linear_problem(rowMatrix, &Y, nonconst_X);

    Amesos Factory;
    string default_solver_type("Amesos_Klu");

    Teuchos::RCP<Amesos_BaseSolver> amesos_solver = Teuchos::rcp( Factory.Create(default_solver_type, linear_problem) );

    Teuchos::ParameterList AmesosList;
    AmesosList.set("PrintTiming",false);
    AmesosList.set("PrintStatus",false);
    AmesosList.set("MaxProcs", 1);

    amesos_solver->SetParameters(AmesosList);
    amesos_solver->SetUseTranspose(false); // needed ? - RWH

    // These could be wrapped in timers - RWH

    amesos_solver->SymbolicFactorization();
    amesos_solver->NumericFactorization();
    amesos_solver->Solve();
  }
  else // this is the real interface we want to support
  {
    *composite_solution_ = *X(0);
    copy_composite_solution_into_problem_data();

    for( size_t i = 0; i < problem_ids_.size(); ++i )
    {
      Model_Evaluator & problem  = *get_problem(problem_ids_[i]);
      EpetraExt::ModelEvaluator::InArgs  inArgs  = problem.createInArgs();
      EpetraExt::ModelEvaluator::OutArgs outArgs = problem.createOutArgs();
      if( inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x) )
      {
        if( problem.supports_preconditioning() )
        {
          problem.apply_preconditioner(*problem.get_x_state(), *(problem_residual_data_[problem_ids_[i]]));
          copy_problem_into_composite(problem_ids_[i], RESIDUAL, Y(0));
        }
        else
        {
          // Simply return the incoming vector which amounts to preconditioning with identity matrix
          *(problem_residual_data_[problem_ids_[i]]) = *problem.get_x_state();
          copy_problem_into_composite(problem_ids_[i], RESIDUAL, Y(0));
        }
      }
    }
  }

#if 0
  cout << std::setprecision(16) << "\nProblem_Manager::ApplyInverse (apply_preconditioner): post-apply Y : " << *Y(0) << endl;
#endif

  return 0;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Teuchos::ParameterList>
Problem_Manager::create_default_status_test_list()
{
  // Create the convergence tests parameter list
  Teuchos::RCP<Teuchos::ParameterList> convergence_paramsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList &convergence_params = *convergence_paramsPtr;

  convergence_params.set("Test Type", "Combo");
  convergence_params.set("Combo Type", "OR");
  convergence_params.set("Number of Tests", 3);
  Teuchos::ParameterList& converged    = convergence_params.sublist("Test 0");
  Teuchos::ParameterList& finite_value = convergence_params.sublist("Test 1");
  Teuchos::ParameterList& maxiters     = convergence_params.sublist("Test 2");
  
  converged.set("Test Type", "Combo");
  converged.set("Combo Type", "AND");
  converged.set("Number of Tests", 1);
  Teuchos::ParameterList& normF_abs = converged.sublist("Test 0");
  normF_abs.set("Test Type", "NormF");
  normF_abs.set("Tolerance", 1.0e-6);
  normF_abs.set("Norm Type", "Two Norm");
  normF_abs.set("Scale Type", "Unscaled");

  finite_value.set("Test Type", "FiniteValue");
  finite_value.set("Vector Type", "F Vector");
  finite_value.set("Norm Type", "Two Norm");
  
  maxiters.set("Test Type", "MaxIters");
  maxiters.set("Maximum Iterations", 5);
  
  return convergence_paramsPtr;
}

//-----------------------------------------------------------------------------

Teuchos::RCP<Teuchos::ParameterList>
Problem_Manager::create_default_nox_params(const Epetra_Comm & comm)
{
  // We'll create a single solver nonlinear parameters list for use with each problem.  For loose
  // coupling, each problem could be configured differently.
  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nl_params.get());
  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", comm.MyPID()); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
                        NOX::Utils::LinearSolverDetails +
			//NOX::Utils::Parameters + 
			//NOX::Utils::Details + 
			NOX::Utils::Warning);
  // start definition of nonlinear solver parameters
  nlParams.sublist("Line Search").set("Method", "Full Step");
  nlParams.sublist("Direction").set("Method", "Newton");
  nlParams.sublist("Direction").sublist("Newton").set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = nlParams.sublist("Direction").sublist("Newton").sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 4);  
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 1);    
  lsParams.set("Aztec Preconditioner", "ilu"); 
  lsParams.set("Preconditioner", "None");

  return nl_params;
}

//-----------------------------------------------------------------------------

NOX::Epetra::Group *
Problem_Manager::create_nox_group( Teuchos::RCP<Teuchos::ParameterList> setupPL,
                                   Teuchos::RCP<Teuchos::ParameterList> nlPL,
                                   Teuchos::RCP<NOX::Epetra::Interface::Required> iReq, 
                                   const NOX::Epetra::Vector& noxSoln)
{
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian>       iJac   = Teuchos::null;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec  = Teuchos::null;
  Teuchos::RCP<Epetra_Operator>                        jacOp  = Teuchos::null;
  Teuchos::RCP<Epetra_Operator>                        precOp = Teuchos::null;

  // A convenient refernece to the linear solver parameter list
  Teuchos::ParameterList & lsParams = nlPL->sublist("Direction").sublist("Newton").sublist("Linear Solver");

  // Create the Jacobian Operator
  string jacobianOp = "";
  if( !Teuchos::is_null(setupPL) )
    jacobianOp = setupPL->get("Jacobian Operator", "Matrix-Free");
  if( ("" == jacobianOp) || ("Matrix-Free" == jacobianOp) )
  {
    jacOp = Teuchos::rcp(new NOX::Epetra::MatrixFree(nlPL->sublist("Printing"), iReq, noxSoln));
    iJac = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(jacOp);
    if( setupPL->isSublist("Matrix-Free") )
    {
      Teuchos::ParameterList & mfParams = setupPL->sublist("Matrix-Free");
      if( mfParams.isParameter("Constant Perturbation Value") )
      {
        double eps = mfParams.get<double>("Constant Perturbation Value");
        Teuchos::rcp_dynamic_cast<NOX::Epetra::MatrixFree>(jacOp)->setPerturbation(eps);
      }
      if( mfParams.isParameter("Lambda Value") )
      {
        double lambda = mfParams.get<double>("Lambda Value");
        Teuchos::rcp_dynamic_cast<NOX::Epetra::MatrixFree>(jacOp)->setLambda(lambda);
      }
    }
    //Teuchos::rcp_dynamic_cast<NOX::Epetra::MatrixFree>(jacOp)->setDifferenceMethod(NOX::Epetra::MatrixFree::Centered);
    //Teuchos::rcp_dynamic_cast<NOX::Epetra::MatrixFree>(jacOp)->setPerturbation(1.e-4);
  }
  else if( "Finite-Difference" == jacobianOp )
  {
    if( setupPL->isSublist("Finite-Difference") )
    {
      Teuchos::ParameterList & fdParams = setupPL->sublist("Finite-Difference");
      if( fdParams.isParameter("Constant Perturbation Value") )
      {
        double eps = fdParams.get<double>("Constant Perturbation Value");
        jacOp = Teuchos::rcp(new NOX::Epetra::FiniteDifference(nlPL->sublist("Printing"), iReq, noxSoln, 0.0, eps));
      }
    }
    else
      jacOp = Teuchos::rcp(new NOX::Epetra::FiniteDifference(nlPL->sublist("Printing"), iReq, noxSoln));

    iJac = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(jacOp);
  }
  else
    throw std::runtime_error("ERROR: There is currently no support for \"Jacobian Operator\" = " + jacobianOp);

  // Create the Preconditioner Operator
  bool construct_with_prec_op = false;
  string precondOp = "";
  if( !Teuchos::is_null(setupPL) )
    precondOp = setupPL->get("Preconditioner Operator", "None"); // default if option not specified in xml
  if( ("" == precondOp) || ("None" == precondOp) )
  {
    lsParams.set("Preconditioner", "None");
    construct_with_prec_op = false;
  }
  else if( "Physics-Based" == precondOp )
  {
    // The user must supply thhe needed objects in the setup parameter list, setupPL
    Teuchos::ParameterList & precOpParams = setupPL->sublist("Preconditioner");
    NOX::Epetra::Interface::Preconditioner * iPrec_ptr  = precOpParams.get<NOX::Epetra::Interface::Preconditioner*>("Interface");
    if( NULL == iPrec_ptr )
      throw std::runtime_error("Could not obtain a valid preconditioner callback interface in setup ParameterList.");

    Epetra_Operator * precOp_ptr  = precOpParams.get<Epetra_Operator *>("Operator");
    if( NULL == precOp_ptr )
      throw std::runtime_error("Could not obtain a valid preconditioner operator in setup ParameterList.");

    iPrec = Teuchos::rcp(iPrec_ptr, false);
    precOp = Teuchos::rcp(precOp_ptr, false);

    lsParams.set("Preconditioner", "User Defined");
    construct_with_prec_op = true;
  }
  else if( "Finite-Difference" == precondOp )
  {
    lsParams.set("Preconditioner", "AztecOO");

    // If Jacobian is Finite-Difference, don't create another one
    if( "Finite-Difference" == jacobianOp )
    {
      lsParams.set("Preconditioner Operator", "Use Jacobian");
      construct_with_prec_op = false;
    }
    else
    {
      precOp = Teuchos::rcp(new NOX::Epetra::FiniteDifference(nlPL->sublist("Printing"), iReq, noxSoln));
      iPrec = Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Preconditioner>(precOp);
      construct_with_prec_op = true;
    }
  }
  else
    throw std::runtime_error("There is currently no support for \"Preconditioner Operator\" = " + precondOp);

  const Epetra_Vector & soln = noxSoln.getEpetraVector();

  // Would like to add support for finite-difference coloring
  // This will likely require adding support to read in an existing matrix-graph.
  {
    //Epetra_CrsMatrix * crsMatrix = NULL;
    //EpetraExt::MatrixMarketFileToCrsMatrix("fd_crs_matrix", comm_, crsMatrix);
    ////cout << "Matrix read from file : " << endl << *crsMatrix << endl;
    ////exit(1);

    //EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = EpetraExt::CrsGraph_MapColoring::GREEDY;
    //EpetraExt::CrsGraph_MapColoring tmpMapColoring( algType, false );
    //Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(const_cast<Epetra_CrsGraph*>(&crsMatrix->Graph()), false);
    //Teuchos::RCP<Epetra_MapColoring> colorMap = Teuchos::rcp(&tmpMapColoring(*graph));
    //cout << "Using " << colorMap->NumColors() << " colors for " << crsMatrix->Map().NumMyElements() << " unknowns" << endl;
    //EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
    //Teuchos::RCP< vector<Epetra_IntVector> > columns = Teuchos::rcp(&colorMapIndex(*graph));

    //Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoring> FDC = 
    //  Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(nlPL->sublist("Printing"),
    //                                                         iReq, soln,
    //                                                         graph,
    //                                                         colorMap, 
    //                                                         columns,
    //                                                         false));
    //Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = FDC;
    ////Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = FDC;
    //linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPL->sublist("Printing"), 
    //                                                  lsParams,
    //                                                  iJac, MF, 
    //                                                  iPrec, FDC, 
    //                                                  *composite_solution_));
  }

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys;

  if( construct_with_prec_op )
    linSys =  Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPL->sublist("Printing"), 
                                                        lsParams,
                                                        iJac, jacOp, 
                                                        iPrec, precOp, 
                                                        soln));
  else
    linSys =  Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPL->sublist("Printing"), 
                                                        lsParams,
                                                        iReq,
                                                        iJac, jacOp, 
                                                        soln));
    
  // Create the Group
  NOX::Epetra::Group * group = new NOX::Epetra::Group(nlPL->sublist("Printing"), iReq, noxSoln, linSys);

  if( !Teuchos::is_null(setupPL) && setupPL->isParameter("Dump Jacobian Operator") )
  {
    string dumpName = setupPL->get<string>("Dump Jacobian Operator");
    dump_operator_to_file(*jacOp, dumpName);
  }

  if( !Teuchos::is_null(setupPL) && setupPL->isParameter("Dump Preconditioner Operator") )
  {
    string dumpName = setupPL->get<string>("Dump Preconditioner Operator");
    dump_operator_to_file(*jacOp, dumpName);
  }

  if( !Teuchos::is_null(setupPL) && setupPL->isParameter("Dump Residual Vector") )
  {
    group->computeF();
    const Epetra_Vector& resid = (dynamic_cast<const NOX::Epetra::Vector&>(group->getF())).getEpetraVector();
    string dumpName = setupPL->get<string>("Dump Residual Vector");
    dump_vector_to_file(resid, dumpName);
  }


  return group;
}
#endif

//-----------------------------------------------------------------------------

Teuchos::RCP<Teuchos::ParameterList>
Problem_Manager::readParameterList( const string & name )
{
  FILE * fp = fopen(name.c_str(), "r");
  if( fp )
  {
    fclose(fp);
    return Teuchos::getParametersFromXmlFile(name);
  }

  return Teuchos::null;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::writeParameterList( const string & name, const Teuchos::ParameterList & plist )
{
  // Output the status test list and nonlinear parameter list to XML
  Teuchos::XMLParameterListWriter xml2pl;
  Teuchos::XMLObject xmllist = xml2pl.toXML(plist);
  std::ofstream file(name.c_str());
  file << xmllist << std::endl;
  file.close();

  return;
}

//-----------------------------------------------------------------------------

string
Problem_Manager::create_date_string()
{
  // get a string representation of the current date/time
  time_t current_time = time(NULL);
  char time_str[26];
  std::strcpy(time_str, ctime(&current_time));
  time_str[24]='\0';
  // ... and then make a legal filename extension out of it
  string time_out(time_str);
  string::size_type i;
  while((i = time_out.find(" ")) != string::npos)
    time_out.replace(i, 1, "_");

  return time_out;
}

//-----------------------------------------------------------------------------

void
Problem_Manager::echo_tee(const string & msg) const
{
  if (m_use_stdout) cout << msg << endl;
  (*out_stream_) << msg << endl;
}

//-----------------------------------------------------------------------------

void 
dump_vector_to_file(const Epetra_Vector & vector, string filename )
{
  EpetraExt::VectorToMatrixMarketFile(filename.c_str(), vector);
  filename += "_map";
  EpetraExt::BlockMapToMatrixMarketFile(filename.c_str(), vector.Map());

  return;
}

//-----------------------------------------------------------------------------

#ifdef HAVE_LIME_JFNK
void
dump_operator_to_file( const Epetra_Operator & op, const string & filename)
{
  static int times_called = 0;
  std::ostringstream suffix;
  suffix << times_called << std::flush;
  string suffixed_name = filename + "_" + suffix.str();

  Teuchos::RCP<const Epetra_RowMatrix> matrix = get_matrix_from_operator(op);

  EpetraExt::RowMatrixToMatrixMarketFile(suffixed_name.c_str(), *matrix);

  times_called++;

  return;
}
#endif

//-----------------------------------------------------------------------------

#ifdef HAVE_LIME_JFNK
Teuchos::RCP<const Epetra_RowMatrix> 
get_matrix_from_operator(const Epetra_Operator & op)
{
  // See if we have an actual matrix
  Teuchos::RCP<const Epetra_RowMatrix> matrix = Teuchos::rcp(dynamic_cast<const Epetra_RowMatrix *>(&op), false);

  if( Teuchos::is_null(matrix) )
  {
    // Or an operator with an underlying matrix
    const NOX::Epetra::FiniteDifference * fd_op = dynamic_cast<const NOX::Epetra::FiniteDifference*>(&op);
    if( fd_op )
      matrix = Teuchos::rcp(&(fd_op->getUnderlyingMatrix()), false);
  }

  if( Teuchos::is_null(matrix) )
  {
    // Right now this works only in serial - RWH
    const Epetra_Map & map = op.OperatorDomainMap();
    Epetra_CrsMatrix * crs_matrix = new Epetra_CrsMatrix(Copy, map, 0);
    Epetra_MultiVector workVec1(map, 1);
    Epetra_MultiVector workVec2(map, 1);

    for( int col = 0; col < map.NumMyElements(); ++col )
    {
      if( col && (col%100 == 0) )
        cout << "\t\tMatrix-Free numerical Jacobian: " << col << " of " << map.NumMyElements() << " dofs completed." << endl;
      workVec1.PutScalar(0.0);
      (*workVec1(0))[col] = 1.0;
      op.Apply(workVec1, workVec2);
      for( int row = 0; row < map.NumGlobalElements(); ++row )
        if( 1.e-10 < fabs((*workVec2(0))[row]) )
          crs_matrix->InsertGlobalValues(row, 1, &(*workVec2(0))[row], &col);
    }
    crs_matrix->FillComplete();
    matrix = Teuchos::rcp(crs_matrix);
  }

  return matrix;
}
#endif

} // namespace LIME

