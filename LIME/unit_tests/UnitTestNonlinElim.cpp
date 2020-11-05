
#ifdef HAVE_LIME_JFNK

// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
// We will incluide the serial comm header regardless; see comments below.
#include "Epetra_SerialComm.h"

namespace {
  const int NUM_ACTIVE_VARS = 2;
}

// -----------------------------------------------------------------------------------
//
// A simple three-equation nonlinear system to test exposed residuals along with 
// nonlinear elimination
//
//               2*x(0)   -   x(1) + x(2)^3   - 31 = 0 
//                 x(0)   + 2*x(1) - 2*x(2)^2 + 21 = 0 
//                 x(0)^2 +   x(1) + 2*x(2)   - 5  = 0 
//
// The system has an exact result of x(0) = 1, x(1) = -2, x(2) = 3
// 
// -----------------------------------------------------------------------------------


// ----------------------------------------------------------------
// A Simple Model Evaluator exposing a system of residual equations
// and which also supports standalone solves.  This version is used
// in the unit test that does not employ Trilinos solvers.  We will
// derive off this one later to disable standalone solves and instead
// use NOX for the second unit test.
// ----------------------------------------------------------------
class ResEQ_ModelEval : public LIME::Model_Evaluator
{
  public:

    ResEQ_ModelEval(LIME::Problem_Manager & pm) :
      LIME::Model_Evaluator(pm, "ResEQ_ModelEval")
    {
      epetra_map_ = Teuchos::rcp( new Epetra_Map(NUM_ACTIVE_VARS, 0, pm.Comm()) );
      soln_vec_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
    }

    virtual ~ResEQ_ModelEval() {}

    // Signal use of standalone solve
    bool supports_standalone_solve() const { return true; }

    // And provide the ability to solve the 2-equation system standalone
    void solve_standalone()
    {
      // The 2-equation system is linear wrt its active variables.
      // The inverse matrix is :
      //    _    _          _      _     _           _  
      //   |  x0  |    1   |  2   1 |   | 31 - x2^3   |
      //   |      | =  - * |        | * |             |
      //   |  x1  |    5   | -1  1  |   | 2*x2^2 - 21 |
      //    -    -          -      -     -           - 

      const double x2 = x2_from_xfer_;
      double rhs0 = 31.0 - pow(x2,3);
      double rhs1 = 2.0*x2*x2 - 21.0;
      //double rhs0 = 7.0 - x2;
      //double rhs1 = 2.0*x2 - 9.0;
      (*soln_vec_)[0] = (2.0*rhs0 + rhs1)/5.0;
      (*soln_vec_)[1] = (2.0*rhs1 - rhs0)/5.0;
      //cout << "Standalone solve for x2 = " << x2 << endl;
      //cout << "x0 = " << (*soln_vec_)[0] << endl;
      //cout << "x1 = " << (*soln_vec_)[1] << endl;
    }

    // Indicate to PM that we expose residuals 
    virtual bool supports_residual() const { return true; }

    // Methods needed to support residual fills
    EpetraExt::ModelEvaluator::InArgs createInArgs() const
    {
      EpetraExt::ModelEvaluator::InArgsSetup inArgs;

      // Register our identify
      inArgs.setModelEvalDescription(my_name_);

      // Signal that we can do calculations with incoming state x
      inArgs.setSupports(IN_ARG_x, true);

      return inArgs;
    }

    EpetraExt::ModelEvaluator::OutArgs createOutArgs() const
    {
      EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

      // Register our identify - consistent with createInArgs
      outArgs.setModelEvalDescription(my_name_);

      // Signal that we can compute a residual vector
      outArgs.setSupports(OUT_ARG_f, true);

      return outArgs;
    }

    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
    {
      // Copy incoming solution state into our own solution vector
      (*soln_vec_) = *(inArgs.get_x());
      const Epetra_Vector & x = *soln_vec_;

      if( outArgs.get_f().get() ) // A non-NULL f-vector signals a residual fill request
      {
        // Get a reference to the vector we will populate
        Epetra_Vector & f = *(outArgs.get_f().get());

        const double x2 = x2_from_xfer_;
        // compute our application's nonlinear residual
        f[0] = 2.0*x[0] - x[1] + pow(x2,3) - 31.0;
        f[1] = x[0] + 2.0*x[1] - 2.0*pow(x2,2) + 21.0;
        //f[0] = 2.0*x[0] - x[1] + x2 - 7.0;
        //f[1] = x[0] + 2.0*x[1] - 2.0*x2 + 9.0;

        //cout << "Doing evalModel with\n\t\t"
        //     << "\tx[0] = " << x[0]
        //     << "\tx[1] = " << x[1]
        //     << "\tx[2] = " << x2_from_xfer_ << endl;
        //cout << "Residual is \n\t\t"
        //     << "\tf[0] = " << f[0]
        //     << "\tf[1] = " << f[1] << endl;
      }
      else
        // We don't support any other callbacks and so should never get
        // called for anything other than a residual fill.
        throw std::runtime_error("ResEQ_ModelEval::evalModel : Only residual fills are supported.");
    }

    Teuchos::RCP<const Epetra_Map> get_x_map() const
      { return epetra_map_; }
    
    Teuchos::RCP<const Epetra_Map> get_f_map() const
      { return epetra_map_; }

    Teuchos::RCP<const Epetra_Vector> get_x_init() const
      { return soln_vec_; }

    Teuchos::RCP<Epetra_Vector> get_x_state()
      { return soln_vec_; }

    Teuchos::RCP<Epetra_Vector> get_x_state() const
      { return soln_vec_; }

    virtual void initializeSolution()
    {
      (*soln_vec_)[0] = 0.0;
      (*soln_vec_)[1] = 0.0;
    }

    // Needed for transfer of data to Elimination module
    virtual double get_x0() const
    { return (*soln_vec_)[0]; }
    virtual double get_x1() const
    { return (*soln_vec_)[1]; }

    // Needed for transfer of data from Elimination module
    virtual void set_x2(double x)
    { x2_from_xfer_ = x; }

  protected:

    string my_name_;
    Teuchos::RCP<Epetra_Map> epetra_map_;
    Teuchos::RCP<Epetra_Vector> soln_vec_;
    double x2_from_xfer_;
};

// ----------------------------------------------------------------
// A Simple Model Evaluator for performing nonlinear elimination
// ----------------------------------------------------------------
class ElimEQ_ModelEval : public LIME::Model_Evaluator,
                         public Elimination_Module
{
  public:

    ElimEQ_ModelEval(LIME::Problem_Manager & pm) :
      LIME::Model_Evaluator(pm, "ElimEQ_ModelEval"),
      Elimination_Module()
  { }

    ~ElimEQ_ModelEval() {}

    // Disable use of standalone solve
    //bool supports_standalone_solve() const { return true; }

    bool perform_elimination()
    {
      const double x0 = x0_from_xfer_;
      const double x1 = x1_from_xfer_;

      my_x_ = (5.0 - x0*x0 - x1) / 2.0;

      //cout << "ElimEQ_ModelEval::perform_elimination with x0 = " << x0 << ", x1 = " << x1 << " and x2 = " << my_x_ << endl;

      return true;
    }

    void initializeSolution()
      { my_x_ = 0.0; }

    // Needed for transfer of data to Residual EQ module
    double get_my_x() const
      { return my_x_; }

    // Needed for transfer of data from Residual EQ module
    void set_x0(double x)
      { x0_from_xfer_ = x; }

    void set_x1(double x)
      { x1_from_xfer_ = x; }

  private:

    const string my_name_;
    double my_x_;
    double x0_from_xfer_;
    double x1_from_xfer_;
};

// -------------------------------------------------------------------------
// A Simple Transfer operator for moving data 
// -------------------------------------------------------------------------
class Elim_To_Res_Transfer : public LIME::Data_Transfer_Operator 
{
  public:

    Elim_To_Res_Transfer(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to) :
      LIME::Data_Transfer_Operator(from, to)
  {}

    bool perform_data_transfer() const
    { 
      const Teuchos::RCP<ElimEQ_ModelEval> fromME = Teuchos::rcp_dynamic_cast<ElimEQ_ModelEval>(source());
      const Teuchos::RCP<ResEQ_ModelEval> toME = Teuchos::rcp_dynamic_cast<ResEQ_ModelEval>(target());
      toME->set_x2( fromME->get_my_x() );
      //cout << "Doing Elim_To_Res_Transfer: sending x2 = " << fromME->get_my_x() << endl;
      return true;
    }
};

// -------------------------------------------------------------------------
// A Simple Transfer operator for moving data
// -------------------------------------------------------------------------
class Res_To_Elim_Transfer : public LIME::Data_Transfer_Operator 
{
  public:

    Res_To_Elim_Transfer(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to) :
      LIME::Data_Transfer_Operator(from, to)
  {}

    bool perform_data_transfer() const
    { 
      const Teuchos::RCP<ResEQ_ModelEval> fromME = Teuchos::rcp_dynamic_cast<ResEQ_ModelEval>(source());
      Teuchos::RCP<ElimEQ_ModelEval> toME = Teuchos::rcp_dynamic_cast<ElimEQ_ModelEval>(target());
      toME->set_x0( fromME->get_x0() );
      toME->set_x1( fromME->get_x1() );
      //cout << "Doing Res_To_Elim_Transfer: sending x0 = " << fromME->get_x0() << " and x1 = " << fromME->get_x1() << endl;
      return true;
    }
};

/* **************************************************************************
   This version of Nonlinear Elimination involves a standalone solve of the
   2-equation system with elimination of the 1-equation problem that occurs
   at the end of each coupling fixed-point iteration when overall
   convergence is assessed.  Fixed-point (any flavor) will not converge and
   this is tested.  No Trilinos solvers are used.
   ************************************************************************** */
TEUCHOS_UNIT_TEST(LIMEUnit, UnitTestNonlinElimSSolve)
{
  // Because we are not using any solvers in Trilinos, we can configure
  // the Problem_Manager using a serial Communicator even when we have
  // built with parallel (mpi) enabled.
  Epetra_SerialComm comm;

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<ResEQ_ModelEval> resME   = Teuchos::rcp(new ResEQ_ModelEval(pm));
  Teuchos::RCP<ElimEQ_ModelEval> elimME = Teuchos::rcp(new ElimEQ_ModelEval(pm));

  pm.add_problem(resME);
  pm.add_problem(elimME);

  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp( new Res_To_Elim_Transfer(resME, elimME) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp( new Elim_To_Res_Transfer(elimME, resME) );

  pm.add_preelimination_transfer(xfer1);
  pm.add_postelimination_transfer(xfer2);

  // This is the default, but we want to make explicit our use of fixed-point coupling
  pm.set_coupling_algorithm(LIME::Problem_Manager::fixed_point);
  pm.set_max_fixed_point_iters(20);

  // Do setup of coupled problem solver
  pm.register_complete(); 

  // Fire off the solve
  pm.integrate();

  double x0 = resME->get_x0();
  double x1 = resME->get_x1();
  double x2 = elimME->get_my_x();

  //cout << "Final result :\n" 
  //     << "\tx[0] = " << x0 << endl
  //     << "\tx[1] = " << x1 << endl
  //     << "\tx[2] = " << x2 << endl;


  // Fixed-point should not converge, and error will actually grow to the
  // point of generating NaNs; we want to test for this.
  bool testNaN0 = Teuchos::ScalarTraits<double>::isnaninf(x0);
  bool testNaN1 = Teuchos::ScalarTraits<double>::isnaninf(x1);
  bool testNaN2 = Teuchos::ScalarTraits<double>::isnaninf(x2);
  TEST_ASSERT(testNaN0);
  TEST_ASSERT(testNaN1);
  TEST_ASSERT(testNaN2);

  // Test that we make it to the end
  TEST_ASSERT(true);
}



/* **************************************************************************
   This version of Nonlinear Elimination has the 2-equation system expose
   its residual equations but NOT support a standalone solve.  This scenario
   demonstrates use of JFNK using NOX with default solver settings.  This is
   the first segue to using solver algorithms (both linear and nonlinear)
   within Trilinos.
   ************************************************************************** */

// First we derive from the ResEQ_ModelEval class and override the
// standalone solve query.  We will also explicitly configure 
// Problem_Manager to perform a jfnk coupling solve instead of fixed-point.
class ResEQ_ModelEval_no_SSolve : public ResEQ_ModelEval
{
  public:

    ResEQ_ModelEval_no_SSolve(LIME::Problem_Manager & pm) :
      ResEQ_ModelEval(pm)
  {
    my_name_ = "ResEQ_ModelEval_no_SSolve";
  }

    ~ResEQ_ModelEval_no_SSolve() {}

    // Disable use of standalone solves
    bool supports_standalone_solve() const { return false; }
};

TEUCHOS_UNIT_TEST(LIMEUnit, UnitTestNonlinElimNOX)
{
 // Now we are using solvers in Trilinos and this requires that we use the 
 // proper mpi communicator depending on hether or not we have done
 // a serial or parallel build.
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<ResEQ_ModelEval_no_SSolve> resME   = Teuchos::rcp(new ResEQ_ModelEval_no_SSolve(pm));
  Teuchos::RCP<ElimEQ_ModelEval> elimME = Teuchos::rcp(new ElimEQ_ModelEval(pm));

  pm.add_problem(resME);
  pm.add_problem(elimME);

  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp( new Res_To_Elim_Transfer(resME, elimME) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp( new Elim_To_Res_Transfer(elimME, resME) );

  pm.add_preelimination_transfer(xfer1);
  pm.add_postelimination_transfer(xfer2);

  // Set JFNK coupling algorithm
  pm.set_coupling_algorithm(LIME::Problem_Manager::jfnk);

  // Do setup of coupled problem solver
  pm.register_complete(); 

  // Fire off the solve
  pm.integrate();

  double x0 = resME->get_x0();
  double x1 = resME->get_x1();
  double x2 = elimME->get_my_x();

  //cout << "Final result :\n" 
  //     << "\tx[0] = " << x0 << endl
  //     << "\tx[1] = " << x1 << endl
  //     << "\tx[2] = " << x2 << endl;


  TEST_FLOATING_EQUALITY(  1.0, x0, 1.e-4);
  TEST_FLOATING_EQUALITY( -2.0, x1, 1.e-4);
  TEST_FLOATING_EQUALITY(  3.0, x2, 1.e-4);

  // Test that we make it to the end
  TEST_ASSERT(true);
}
#endif
