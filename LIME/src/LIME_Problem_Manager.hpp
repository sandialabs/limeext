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
// LIME_Problem_Manager.hpp
//
// Description: 
// 
// ********************************************************************************************

#ifndef NOX_EPETRA_PROBLEM_MANAGER_H
#define NOX_EPETRA_PROBLEM_MANAGER_H

#include "LIME_config.hpp"

#include <ostream>
#include <vector>
#include <utility>

#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "LIME_Elimination_Module.hpp"
#include "LIME_Model_Evaluator.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifdef HAVE_LIME_JFNK
#include "NOX_Epetra.H"
#endif

using std::map;
using std::multimap;
using std::pair;
using std::set;
using std::vector;
using std::runtime_error;

//-----------------------------------------------------------------------------

class Amesos_BaseSolver;

// These would be better as template parameters, but for now ....
// Actually, these would be better heap allocated, not as template
// parameters.  Template parameters are compile time values and a
// realistic physics simulation may be very dynamic meaning the 
// actual number of hierarchical solves may be problem (input)
// dependent and not known at compile time.  For example, I'd like
// to use LIME for Melcor but the numbers of nested solves and transfers
// won't be known until we read the NPP input deck so fixing nesting
// depth as template parameters is not a good idea.

#define LIME_MAX_HIERARCHIES 100
#define LIME_MAX_PROBLEMS_PER_HIERARCHY 100

//-----------------------------------------------------------------------------

namespace LIME {

  // A utility that might be better in a LIME_Utils file.
  template<typename T>
    string
    convert_to_string( T id )
    { 
      std::ostringstream name;
      name << id << std::flush;
      return name.str();
    }

  // More candidates for a utility
  Teuchos::RCP<const Epetra_RowMatrix> get_matrix_from_operator(const Epetra_Operator &);
  void dump_vector_to_file( const Epetra_Vector &, string filename );
#ifdef HAVE_LIME_JFNK
  void dump_operator_to_file( const Epetra_Operator &, const string & filename );
#endif
  inline bool has_time_step_changed(double dt1, double dt2)
    { return (abs(dt1 - dt2) > 1.e-10); }


  class Data_Transfer_Operator;
  class FixedPoint_Accelerator;

  typedef unsigned ID_TYPE;

  /*! A manager class for a collection of registered problems.
      
      This class drives various coupling algorithm strategies.  Current
      support for various time integrations and/or operator splitting
      is limited.
   */

  class Problem_Manager 
#ifdef HAVE_LIME_JFNK
                        : public EpetraExt::ModelEvaluator,
                          public Epetra_Operator
#endif
  {
    public:

      friend class LIME::Model_Evaluator;

      /*! \brief Used to specify which composite vector to copy from or to 
       */
      enum VECTOR_TYPE { SOLUTION_X, RESIDUAL };

      /*! \brief Used to specify which coupling strategy to employ
       */
      typedef enum {
            jfnk = 0
          , fixed_point
          , switching
      } solve_using ;

      /*! \brief Used to specify the context of a data tranfer between applications being coupled
       */
      typedef enum {
        fixed_point_xfer = 10
          , successive_xfer
      } transfer_types;

      /*! \brief Used to specify the context of a data tranfer between applications being coupled
       */
      typedef enum {
        jacobi = 100 ,
        seidel
      } fixed_point_mode;

      Problem_Manager(const Epetra_Comm& Comm, bool verbose=false, bool use_stdout = true);

      virtual ~Problem_Manager();

      void problem_has_stateful_transfers();

      /*! \brief Add a problem to the Manager

        Returns an ID_TYPE identifier used as a key in various inter-problem interactions.
       */
      virtual ID_TYPE add_problem(Teuchos::RCP<LIME::Model_Evaluator>);

      /*! \brief Get the command line options for this Problem_Manager
       */
      Teuchos::CommandLineProcessor & get_command_line_options()
      { return m_clp; }

      /*! \brief Get the verbosity setting
       */
      unsigned get_verbosity() const
      { return m_verbose; }

      /*! \brief Get the hierarchy level for this manager
       */
      unsigned get_hierarchy() const
      { return my_hierarchy_; }

      /*! \brief Get a problem from the Manager given its integer id
       */
      Teuchos::RCP<LIME::Model_Evaluator> & get_problem(ID_TYPE);
      const Teuchos::RCP<LIME::Model_Evaluator> & get_problem(ID_TYPE) const;

      /*! \brief Get all leaf problems for this manager
       */
      vector<ID_TYPE> get_all_descendants() const;

      /*! \brief Get a problem from the Manager given its integer id
       */
       const multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> > & get_all_registered_transfers() const
       {
         if( 0 != my_hierarchy_ )
           throw runtime_error("Can only call get_all_registered_transfers on top-level Problem_Manager.");

         return all_registered_transfers_;
       }

      /*! \brief Allow user to specify a setup XML file other than Problem_Manager_setup.xml
       */
      void set_xml_file(const string & name)
      { m_setupXMLfilename = name; }

      /*! \brief Allow user to specify a parameters configuration XML file for exposing parameters
       */
      void set_params_xml_file(const string & name)
      { m_paramsXMLfilename = name; }

      /*! \brief Allow user to directly specify coupling algorithm
       */
      void set_coupling_algorithm(solve_using su)
      { 
        verify_not_register_complete(); 
        su_ = su;
      }

      /*! \brief Allow user to directly specify max fixed-point iterations
       */
      void set_max_fixed_point_iters(int iters)
      { max_fp_iters_ = iters; }

      /*! \brief Allow user to directly specify fixed-point convergence tolerance
       */
      void set_fixed_point_tolerance(double tol)
      { fp_abs_tolerance_ = tol; }

      /*! \brief Allow user to directly specify fixed-point mode
       */
      void set_fixed_point_mode(fixed_point_mode mode)
      { fp_solve_mode_ = mode; }

      /*! \brief Associate a name with a problem parameter vector
       */

      /*! \brief Associate a name with a problem response vector
       */
      void post_problem_response(const string &, Teuchos::RCP<LIME::Model_Evaluator>, Teuchos::RCP<Epetra_Vector>);

      // Add normal data transfer
      void add_transfer(Teuchos::RCP<Data_Transfer_Operator>&);

      // Add specialized data transfer
      void add_preelimination_transfer(Teuchos::RCP<Data_Transfer_Operator>);

      // Add specialized data transfer
      void add_postelimination_transfer(Teuchos::RCP<Data_Transfer_Operator>);

      /*! \brief Signal completion of setup

        This call triggers construction of nonlinear solvers
        and other data and operators as appropriate.
       */
      void register_complete();

      /*! \brief Configure parameters for problems involving them
       */
      void configure_model_parameters(const string &);

      /*! \brief Refresh composite solution from registered problem solutions
       */
      void refresh_solution_state() const;

      /*! \brief Copy vectors from each problem into a composite problem vector
       */
      void copy_problems_into_composite(VECTOR_TYPE) const;

      /*! \brief Copy a vector from a single problem into the appropriate portion of a composite vector
       */
      void copy_problem_into_composite(ID_TYPE prob_id, VECTOR_TYPE) const;

      /*! \brief Copy a vector from a single problem into the appropriate portion of an arbitrary vector sized conformal to the composite problem vector
       */
      void copy_problem_into_composite(ID_TYPE prob_id, VECTOR_TYPE, Epetra_Vector *) const;

      /*! \brief Copy an arbitrary vector sized conformal to prob_id into the appropriate portion the composite solution vector
       */
      void copy_vector_into_composite_solution(ID_TYPE prob_id, const Epetra_Vector &);

      /*! \brief Copy composite solution vector into all associated problem data vectors
       */
      void copy_composite_solution_into_problem_data() const;

      /*! \brief Copy composite solution vector into a specified problem data vector
       */
      void copy_composite_solution_into_problem_data(ID_TYPE) const;

      /*! \brief Create a view into the appropriate portion of the specified composite vector

        Because this allows direct modification of the composite vector, it should be used with care. 
       */
      const Epetra_Vector * create_view_into_composite_vector(ID_TYPE prob_id, VECTOR_TYPE);

      /*! \brief Configure the nonlinear solver used to advance each time step
       */
      void setup_nonlinear_solve();

      /*! \brief Solve the nonlinear problem for a given time step
       */
      bool solve_nonlinear(double & dt);

      /*! \brief Configure a fixed-point coupling strategy to solve the coupled system
       */
      void setup_fixed_point();

      /*! \brief Solve the coupled system using a fixed-point strategy

        Returns true if converged, false if not.
       */
      bool solve_fixed_point(double & dt);

      /*! \brief Determine overall convergence status of the coupled system solved using fixed-point
       */
      bool compute_fixed_point_convergence();

      /*! \brief Configure a Jacobian-Free Newton-Krylov coupling strategy to solve the coupled system
       */
      void setup_jfnk();
      /*! \brief Solve the coupled system using a Jacobian-Free Newton-Krylov strategy

        Returns true if converged, false if not.
       */
      bool solve_jfnk(double & dt);

      /*! \brief Configure a hybrid "Switching" coupling strategy to solve the coupled system

        Details are configured via XML input
       */
      void setup_switching();
      /*! \brief Solve the coupled system using a hybrid "Switching" coupling strategy

        Returns true if converged, false if not.
       */
      bool solve_switching(double & dt);

      /*! \brief Compute predicted solutions for each problem at the beginning of each time step
       */
      void apply_predictors(solve_using);

      /*! \brief Update scaling vectors used to reduce differneces in magnitude of problems variables 
       */
      void update_right_scaling_vectors();

      /*! \brief Begin solution of the coupled problem 
       */
      void integrate();
      void const_integrate() const;

      /*! \brief Provide a summary of the coupled system
       */
      void output_status( ostream & os );

      /*! \brief Verify consistency of current time state of all problems
       */
      bool verify_consistent_current_time( double );

      static Teuchos::RCP<Teuchos::ParameterList> readParameterList( const string & filename );
      static void writeParameterList( const string & filename, const Teuchos::ParameterList &);
      static string create_date_string();

#ifdef HAVE_LIME_JFNK
      // Inherited from NOX::Epetra::Interface::Required and allows us
      // to be a callback interface to a coupled problem NOX solver
      virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const NOX::Epetra::Interface::Required::FillType fillFlag) const;

      // Allows composite residual fill from contexts outside NOX solver
      bool compute_composite_residual() const;

      // Inherited from NOX::Epetra::Interface::Preconditioner and allows us
      // to be a callback interface for physics-based preconditioning
      virtual bool computePreconditioner(const Epetra_Vector& x, 
          Epetra_Operator& M,
          Teuchos::ParameterList* precParams = 0) const;

      // Inherited from Epetra_Operator and allows us to handle the 
      // ApplyInverse method in which we perform physics-based preocnditrioning

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      // Pure virtuals inherited from EpetraExt::ModelEvaluator.
      virtual InArgs createInArgs() const;
      virtual OutArgs createOutArgs() const;
      virtual Teuchos::RCP<const Epetra_Map> get_x_map() const
      { return composite_map_; }
      virtual Teuchos::RCP<const Epetra_Map> get_f_map() const
      { return composite_map_; }
      virtual Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
      virtual Teuchos::RCP<const Epetra_Vector> get_x_init() const
      { return composite_solution_; }
      //virtual Teuchos::RCP<Epetra_Operator> create_W() const
      //{ return Teuchos::null; }
      virtual Teuchos::RCP<const Epetra_Vector> get_p_init(int) const;
      virtual Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
      virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;
#endif

      void publish_problem_parameters(Teuchos::RCP<const Epetra_Vector>) const;

      
#ifdef HAVE_LIME_JFNK
      Teuchos::RCP<NOX::Epetra::Group> get_nox_group()
      { return grp_; }

      /*! \brief Get the solver parameter list (valid only for JFNK solves)
       */
      const Teuchos::ParameterList & get_jfnk_solver_params() const;
#endif

      void compute_problem_responses() const;

      const map<string, pair<Teuchos::RCP<LIME::Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > > & get_responses() const
      { return m_responses; }
      
      // Also some boilerplate to fulfill pure virtual methods of Epetra_Operator
#ifdef HAVE_LIME_JFNK
      virtual int SetUseTranspose(bool UseTranspose)
      { return FD->SetUseTranspose(UseTranspose); }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      virtual double NormInf() const
      { return FD->NormInf(); }
#endif

      virtual const char * Label() const
      { static const char * me = "Physics-Based via NOX::Epetra::Problem_Manager"; return me; }

#ifdef HAVE_LIME_JFNK
      virtual bool UseTranspose() const
      { return FD->UseTranspose(); }

      virtual bool HasNormInf() const
      { return FD->HasNormInf(); }
#endif

      virtual const Epetra_Comm & Comm() const
      { return comm_; }

      virtual const Epetra_Map & OperatorDomainMap() const
      { return *composite_map_; }

      virtual const Epetra_Map & OperatorRangeMap() const
      { return *composite_map_; }

    protected:

      const Epetra_Comm & comm_;
      Teuchos::CommandLineProcessor m_clp;
      bool m_verbose;
      string m_setupXMLfilename;
      string m_paramsXMLfilename;
      bool use_predictor_;
      map<unsigned, unsigned> num_hierarchy_problems_; 
      unsigned my_hierarchy_; 
      int num_total_dofs_; 
      solve_using su_;
      transfer_types transfer_type;
      bool register_complete_called_;

      string pl_dir_name_;
      Teuchos::RCP<std::ostream> out_stream_;

      Teuchos::RCP<Epetra_Map> composite_map_;
      mutable Teuchos::RCP<Epetra_Vector> composite_solution_; 
      mutable Teuchos::RCP<Epetra_Vector> composite_residual_; 
      set<ID_TYPE> provisional_problem_ids_;
      vector<ID_TYPE> problem_ids_;
      mutable map<ID_TYPE, Teuchos::RCP<LIME::Model_Evaluator> > problems_;
      map<ID_TYPE, Teuchos::RCP<Teuchos::ParameterList> > nl_params_;
#ifdef HAVE_LIME_JFNK
      map<ID_TYPE, Teuchos::RCP<NOX::StatusTest::Generic> > status_test_;
#endif
      map<ID_TYPE, string> id_to_name_;
      map<string, ID_TYPE> name_to_id_;

#ifdef HAVE_LIME_JFNK
      Teuchos::RCP<NOX::Epetra::Group> grp_;
#endif
      // A wrapper to present myself to NOX as an interface originating from a ModelEvaluator
#ifdef HAVE_LIME_JFNK
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> my_nox_interface_;
#endif
      multimap<ID_TYPE, ID_TYPE> problem_dependencies_;
      mutable map<ID_TYPE, Teuchos::RCP<Epetra_Vector> > problem_residual_data_; 
      mutable map<ID_TYPE, Teuchos::RCP<Epetra_IntVector> > problem_LID_to_composite_GID_map_;

#ifdef HAVE_LIME_JFNK
      // Needed for use with physics-based preconditioning work - RWH
      mutable Teuchos::RCP<NOX::Epetra::FiniteDifference> FD;
      Teuchos::RCP<NOX::Solver::Generic> composite_nox_solver_;
      Teuchos::RCP<NOX::StatusTest::Generic> status_tests;
#endif

      // data for fixed point solve
      fixed_point_mode fp_solve_mode_;
      int max_fp_iters_;
      double fp_abs_tolerance_;
      vector<ID_TYPE> standalone_prob_ids;
      map<ID_TYPE, Teuchos::RCP<LIME::Model_Evaluator> > standalone_solvers_;
#ifdef HAVE_LIME_JFNK
      map<ID_TYPE, Teuchos::RCP<NOX::Solver::Generic> > standalone_nox_solvers_;
#endif
      vector<Teuchos::RCP<Data_Transfer_Operator> > standard_transfers_;
      Teuchos::RCP<FixedPoint_Accelerator> fp_accelerator_;

      // These only get populated by the top-level PM
      map<ID_TYPE, Teuchos::RCP<LIME::Model_Evaluator> > all_registered_problems_;
      multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> > all_registered_transfers_;

      // Data for supporting nonlinear elimination with specific pre- and post-elimination
      // transfer operation requirements
      map<ID_TYPE, Teuchos::RCP<Elimination_Module> > elimination_modules_;
      vector<Teuchos::RCP<Data_Transfer_Operator> > preelimination_transfers_;
      vector<Teuchos::RCP<Data_Transfer_Operator> > postelimination_transfers_;

      // Stopping criterion for time stepping
      unsigned int max_time_steps_;

      // data for switching solver
      vector<Teuchos::RCP<Teuchos::ParameterList> > switching_solve_params_;

      // Helper creation functions
#ifdef HAVE_LIME_JFNK
      static NOX::Epetra::Group * create_nox_group(
          Teuchos::RCP<Teuchos::ParameterList> setupPL,
          Teuchos::RCP<Teuchos::ParameterList> nlPL,
          Teuchos::RCP<NOX::Epetra::Interface::Required>,
          const NOX::Epetra::Vector& soln);
#endif
      static Teuchos::RCP<Teuchos::ParameterList> create_default_status_test_list();

    public:
#ifdef HAVE_LIME_JFNK
      static Teuchos::RCP<Teuchos::ParameterList> create_default_nox_params(const Epetra_Comm &);
#endif

      // Get num problems for a hierarchy level
      unsigned num_problems(unsigned) const;

      // Get num problems for a hierarchy level
      virtual void increment_num_problems(unsigned hier);

      // Get unique problem id
      virtual bool get_next_problem_id(unsigned hier, ID_TYPE & new_id) const;

      // An output utility to send output to multiple streams
      void echo_tee(const string &) const;

    protected:

      // Add command line arguments with supported options
      void add_command_line_options();

      // Add a provisional problem id to the Manager
      virtual void add_provisional_problem_id(ID_TYPE id)
      { provisional_problem_ids_.insert(id); }

      void verify_not_register_complete() const
      { if(register_complete_called_) throw runtime_error("ERROR: register_complete has already been called."); }

      double determine_time_step() const;

      // these functions handle all dependency checking, ordering
      // and the actual running of the pre / perform / post operations
      // within the context of nonlinear elimination

      void transfer_setup(vector<ID_TYPE>& ids, ID_TYPE id);
      void transfer_details(vector<Teuchos::RCP<Data_Transfer_Operator> >& ops
          , vector<ID_TYPE>& sids, vector<ID_TYPE>& tids);
      void perform_transfers(vector<Teuchos::RCP<Data_Transfer_Operator> > const &) const;
      void perform_all_transfers() const;
      Teuchos::RCP<Data_Transfer_Operator> find_pre_target_xfer_operator(ID_TYPE id) const;
      Teuchos::RCP<Data_Transfer_Operator> find_post_source_xfer_operator(ID_TYPE id
          , vector< Teuchos::RCP<Data_Transfer_Operator> >& pets) const;

      vector< Teuchos::RCP<Elimination_Module> > elimination_modules;
      vector<ID_TYPE> pre_sids, pre_tids;
      vector<ID_TYPE> post_sids, post_tids;
      typedef pair<Teuchos::RCP<Data_Transfer_Operator>, Teuchos::RCP<Elimination_Module> > pair_om;
      vector<pair_om> perform_pres, post_performs;
      bool residual_active_;
      bool problem_is_transient;

      Teuchos::RCP<Epetra_LocalMap> m_model_p_map;
      Teuchos::RCP<Epetra_Vector> m_model_p_vec;
      Teuchos::RCP<Teuchos::Array<std::string> > m_model_p_names;
      mutable multimap<int, pair<int,int> > m_model_p_mappings;
      map<string, pair<Teuchos::RCP<LIME::Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > > m_responses;

      bool stateful_transfers;
      bool m_use_stdout;
  };

}

#endif
