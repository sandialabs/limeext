
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

using LIME::Data_Transfer_Operator;
using LIME::Hierarchical_Problem_Manager;
using LIME::Model_Evaluator;
using LIME::Problem_Manager;

const double FP_UNIT_TEST_TOL = 1.e-4;

// -------------------------------------------------------------
// A Simple Model Evaluator for each equation in a linear system
// -------------------------------------------------------------
class EQ_ModelEval : public Model_Evaluator
{
  public:

    EQ_ModelEval(Problem_Manager & pm, double (*f_ptr)(double), bool verbose = false) :
      Model_Evaluator(pm, "EQ_ModelEval"),
      verbose_(verbose),
      my_func_(f_ptr)
  {
    // Assigning these values will test that initializeSolution gets called by PM
    my_x_ = -1.e20;
    other_x_ = -1.e20;
  }

    ~EQ_ModelEval() {}
    bool supports_standalone_solve() const { return true; }

    void initializeSolution() {
      my_x_ = 0.0;
      other_x_ = 0.0;
    }

    void reinitialize() {
      iterate_vals_.clear();
      initializeSolution();
      iterate_vals_.push_back(my_x_);
    }

    void solve_standalone() {
      if(verbose_)
        cout << "Doing solve with other_x_ = " << other_x_ << " ... ";
      my_x_ = (*my_func_)(other_x_);
      if(verbose_)
        cout << "x = " << my_x_ << endl;
      // Store our solution states for later retrieval and comparison
      iterate_vals_.push_back(my_x_);
    }

    bool is_converged() { 
      return ( fabs(1.0 - my_x_) < FP_UNIT_TEST_TOL );
    }

    const double & get_my_x() { return my_x_; }
    void set_other_x(double x) { other_x_ = x; }
    const vector<double> & get_my_iterates() { return iterate_vals_; }

  private:

    bool verbose_;
    double my_x_;
    double other_x_;
    double (*my_func_)(double);
    vector<double> iterate_vals_;
};

// -------------------------------------------------------------------------
// A Simple Transfer operator for moving data between EQ_ModelEval instances
// -------------------------------------------------------------------------
class Simple_EQ_Transfer : public Data_Transfer_Operator 
{
  public:

    Simple_EQ_Transfer(Teuchos::RCP<Model_Evaluator> from, Teuchos::RCP<Model_Evaluator> to, bool verbose = false) :
      Data_Transfer_Operator(from, to),
      verbose_(verbose)
  {}

    bool perform_data_transfer() const
    { 
      const Teuchos::RCP<EQ_ModelEval> & fromME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(source());
      const Teuchos::RCP<EQ_ModelEval> & toME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(target());
      toME->set_other_x( fromME->get_my_x() );
      if(verbose_)
        cout << "Transferring " << fromME->get_my_x() << " from " << fromME->get_id() << " to " << toME->get_id() << endl;
      return true;
    }

  private:

    bool verbose_;
};


// -----------------------------------------------------------------------------------
// A Two-equation linear system to test Jacobi and Seidel mode fixed-point
// 
//    This is a contrived linear system broken into two standalone equations (ME s):
// 
//                     4.0*x1 - 3.0*x2 =  1.0
//                         x1 - 2.0*x2 = -1.0
// 
//    which should produce the solution (1.0, 1.0) to within tolerance.
// -----------------------------------------------------------------------------------
namespace
{
  double solve_x1(double x2)
  {
    double x1 = (1.0 + 3.0*x2)/4.0;
    return x1;
  }

  double solve_x2(double x1)
  {
    double x2 = (1.0 + x1)/2.0;
    return x2;
  }

}


TEUCHOS_UNIT_TEST(UnitTestFixedPoint, testUnit1)
{
  Epetra_SerialComm comm;

  bool verbose = false;
  Problem_Manager pm(comm, verbose);

  Teuchos::RCP<EQ_ModelEval> eq1ME = Teuchos::rcp(new EQ_ModelEval(pm, solve_x1));
  Teuchos::RCP<EQ_ModelEval> eq2ME = Teuchos::rcp(new EQ_ModelEval(pm, solve_x2));

  pm.add_problem(eq1ME);
  pm.add_problem(eq2ME);

  Teuchos::RCP<Data_Transfer_Operator> xfer1 = Teuchos::rcp( new Simple_EQ_Transfer(eq1ME, eq2ME) );
  Teuchos::RCP<Data_Transfer_Operator> xfer2 = Teuchos::rcp( new Simple_EQ_Transfer(eq2ME, eq1ME) );

  pm.add_transfer(xfer1);
  pm.add_transfer(xfer2);

  pm.register_complete(); // Trigger setup of coupled problem solver

  // Create gold results for later comparison
  size_t jacobi_iters = 0, seidel_iters = 0;
  vector<double> jacobi_x1, jacobi_x2, seidel_x1, seidel_x2;
  {
    double x1 = 0.0, x2 = 0.0;
    jacobi_x1.push_back(x1);
    jacobi_x2.push_back(x2);
    seidel_x1.push_back(x1);
    seidel_x2.push_back(x2);
    for( size_t i = 1; i < 25; ++i )
    {
      x1 = solve_x1(jacobi_x2[i-1]);
      x2 = solve_x2(jacobi_x1[i-1]);
      jacobi_x1.push_back(x1);
      jacobi_x2.push_back(x2);
      if( 0 == jacobi_iters )
      {
        double test = fabs(1.0-x1);
        test = std::max(test, fabs(1.0-x2));
        if( test < FP_UNIT_TEST_TOL )
          jacobi_iters = i;
      }

      x1 = solve_x1(seidel_x2[i-1]);
      x2 = solve_x2(x1);
      seidel_x1.push_back(x1);
      seidel_x2.push_back(x2);
      if( 0 == seidel_iters )
      {
        double test = fabs(1.0-x1);
        test = std::max(test, fabs(1.0-x2));
        if( test < FP_UNIT_TEST_TOL )
          seidel_iters = i;
      }
    }
  }
  // Need to include initial guess in count of ME iterates
  jacobi_iters++;
  seidel_iters++;

  // Allow sufficient iterations for either Jacobi or Seidel modes
  pm.set_max_fixed_point_iters(25);

  // The default fixed-point mode should be Jacobi; we will test and see
  pm.integrate();

  const vector<double> & eq1_iterates = eq1ME->get_my_iterates();
  const vector<double> & eq2_iterates = eq2ME->get_my_iterates();

  TEST_EQUALITY(eq1_iterates.size(), eq2_iterates.size());
  TEST_EQUALITY(jacobi_iters, eq1_iterates.size());
  double error = 0.0;
  for( size_t i = 0; i < eq1_iterates.size(); ++i )
    error += std::max(fabs(eq1_iterates[i]-jacobi_x1[i]), fabs(eq2_iterates[i]-jacobi_x2[i]));
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12 );

  // Now solve using Seidel fixed-point 
  pm.set_fixed_point_mode(Problem_Manager::seidel);
  pm.integrate(); // also tests that problems get reinitialized

  TEST_EQUALITY(eq1_iterates.size(), eq2_iterates.size());
  TEST_EQUALITY(seidel_iters, eq1_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq1_iterates.size(); ++i )
    error += std::max(fabs(eq1_iterates[i]-seidel_x1[i]), fabs(eq2_iterates[i]-seidel_x2[i]));
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12 );


  // Now test that Jacobi convergenece behavior is invariant wrt problem ordering
  // and that Seidel IS affected by order of problem registration
  Problem_Manager pm2(comm, verbose);

  Teuchos::RCP<EQ_ModelEval> eq3ME = Teuchos::rcp(new EQ_ModelEval(pm2, solve_x1));
  Teuchos::RCP<EQ_ModelEval> eq4ME = Teuchos::rcp(new EQ_ModelEval(pm2, solve_x2));

  Teuchos::RCP<Data_Transfer_Operator> xfer3 = Teuchos::rcp( new Simple_EQ_Transfer(eq3ME, eq4ME) );
  Teuchos::RCP<Data_Transfer_Operator> xfer4 = Teuchos::rcp( new Simple_EQ_Transfer(eq4ME, eq3ME) );

  // Flip the previous order of problem registration
  pm2.add_problem(eq4ME);
  pm2.add_problem(eq3ME);

  // Flip the order of transfer registration; shouldn't matter
  pm2.add_transfer(xfer4);
  pm2.add_transfer(xfer3);

  pm2.register_complete(); // Trigger setup of coupled problem solver

  pm2.set_max_fixed_point_iters(25); // Allow sufficient iterations

  // The default fixed-point mode should be Jacobi; we will test and see
  pm2.integrate();

  const vector<double> & eq3_iterates = eq3ME->get_my_iterates();
  const vector<double> & eq4_iterates = eq4ME->get_my_iterates();

  TEST_EQUALITY(eq3_iterates.size(), eq4_iterates.size());
  TEST_EQUALITY(jacobi_iters, eq3_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq3_iterates.size(); ++i )
    error += std::max(fabs(eq3_iterates[i]-jacobi_x1[i]), fabs(eq4_iterates[i]-jacobi_x2[i]));
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12);

  // Now solve using Seidel fixed-point 
  pm2.set_fixed_point_mode(Problem_Manager::seidel);
  pm2.integrate(); // also tests that problems get reinitialized

  TEST_EQUALITY(eq3_iterates.size(), eq4_iterates.size());
  TEST_INEQUALITY(seidel_iters, eq3_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq3_iterates.size(); ++i )
    error += std::max(fabs(eq3_iterates[i]-seidel_x1[i]), fabs(eq4_iterates[i]-seidel_x2[i]));
  error += 1.0; // shift needed for equality test
  TEST_INEQUALITY( 1.0, error);

  // Test that we make it to the end
  TEST_ASSERT(true);
}


TEUCHOS_UNIT_TEST(UnitTestFixedPoint, testUnit2)
{
  Epetra_SerialComm comm;

  bool verbose = false;
  Problem_Manager pm(comm, verbose);

  Teuchos::RCP<Hierarchical_Problem_Manager> pm2 = Teuchos::rcp(new Hierarchical_Problem_Manager(pm, verbose));
  Teuchos::RCP<Hierarchical_Problem_Manager> pm3 = Teuchos::rcp(new Hierarchical_Problem_Manager(pm, verbose));

  Teuchos::RCP<EQ_ModelEval> eq1ME = Teuchos::rcp(new EQ_ModelEval(*pm2, solve_x1));
  Teuchos::RCP<EQ_ModelEval> eq2ME = Teuchos::rcp(new EQ_ModelEval(*pm3, solve_x2));

  pm2->add_problem(eq1ME);
  pm3->add_problem(eq2ME);
  pm2->register_complete();
  pm3->register_complete();

  pm.add_problem(pm2);
  pm.add_problem(pm3);

  Teuchos::RCP<Data_Transfer_Operator> xfer1 = Teuchos::rcp( new Simple_EQ_Transfer(eq1ME, eq2ME) );
  Teuchos::RCP<Data_Transfer_Operator> xfer2 = Teuchos::rcp( new Simple_EQ_Transfer(eq2ME, eq1ME) );

  pm.add_transfer(xfer1);
  pm.add_transfer(xfer2);

  pm.register_complete(); // Trigger setup of coupled problem solver

  // Create gold results for later comparison
  size_t jacobi_iters = 0, seidel_iters = 0;
  vector<double> jacobi_x1, jacobi_x2, seidel_x1, seidel_x2;
  {
    double x1 = 0.0, x2 = 0.0;
    jacobi_x1.push_back(x1);
    jacobi_x2.push_back(x2);
    seidel_x1.push_back(x1);
    seidel_x2.push_back(x2);
    for( size_t i = 1; i < 25; ++i )
    {
      x1 = solve_x1(jacobi_x2[i-1]);
      x2 = solve_x2(jacobi_x1[i-1]);
      jacobi_x1.push_back(x1);
      jacobi_x2.push_back(x2);
      if( 0 == jacobi_iters )
      {
        double test = fabs(1.0-x1);
        test = std::max(test, fabs(1.0-x2));
        if( test < FP_UNIT_TEST_TOL )
          jacobi_iters = i;
      }

      x1 = solve_x1(seidel_x2[i-1]);
      x2 = solve_x2(x1);
      seidel_x1.push_back(x1);
      seidel_x2.push_back(x2);
      if( 0 == seidel_iters )
      {
        double test = fabs(1.0-x1);
        test = std::max(test, fabs(1.0-x2));
        if( test < FP_UNIT_TEST_TOL )
          seidel_iters = i;
      }
    }
  }
  // Need to include initial guess in count of ME iterates
  jacobi_iters++;
  seidel_iters++;

  // Allow sufficient iterations for either Jacobi or Seidel modes
  pm.set_max_fixed_point_iters(25);

  // The default fixed-point mode should be Jacobi; we will test and see
  pm.integrate();

  const vector<double> & eq1_iterates = eq1ME->get_my_iterates();
  const vector<double> & eq2_iterates = eq2ME->get_my_iterates();

  TEST_EQUALITY(eq1_iterates.size(), eq2_iterates.size()+1);
  TEST_EQUALITY(jacobi_iters, eq1_iterates.size());
  double error = 0.0;
  for( size_t i = 0; i < eq1_iterates.size(); ++i )
    error += fabs(eq1_iterates[i]-jacobi_x1[i]);
  for( size_t i = 0; i < eq2_iterates.size(); ++i )
    error += fabs(eq2_iterates[i]-jacobi_x2[i]);
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12);

  // Now solve using Seidel fixed-point 
  pm.set_fixed_point_mode(Problem_Manager::seidel);
  pm.integrate(); // also tests that problems get reinitialized

  TEST_EQUALITY(eq1_iterates.size(), eq2_iterates.size()+1);
  TEST_EQUALITY(seidel_iters, eq1_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq1_iterates.size(); ++i )
    error += fabs(eq1_iterates[i]-seidel_x1[i]);
  for( size_t i = 0; i < eq2_iterates.size(); ++i )
    error += fabs(eq2_iterates[i]-seidel_x2[i]);
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12);


  // Now test that Jacobi convergenece behavior is invariant wrt problem ordering
  // and that Seidel IS affected by order of problem registration
  Problem_Manager pm5(comm, verbose);

  Teuchos::RCP<Hierarchical_Problem_Manager> pm6 = Teuchos::rcp(new Hierarchical_Problem_Manager(pm5, verbose));
  Teuchos::RCP<Hierarchical_Problem_Manager> pm7 = Teuchos::rcp(new Hierarchical_Problem_Manager(pm5, verbose));

  Teuchos::RCP<EQ_ModelEval> eq3ME = Teuchos::rcp(new EQ_ModelEval(*pm6, solve_x1));
  Teuchos::RCP<EQ_ModelEval> eq4ME = Teuchos::rcp(new EQ_ModelEval(*pm7, solve_x2));

  pm6->add_problem(eq3ME);
  pm7->add_problem(eq4ME);
  pm6->register_complete();
  pm7->register_complete();

  // Flip the previous order of problem registration
  pm5.add_problem(pm7);
  pm5.add_problem(pm6);

  Teuchos::RCP<Data_Transfer_Operator> xfer3 = Teuchos::rcp( new Simple_EQ_Transfer(eq3ME, eq4ME) );
  Teuchos::RCP<Data_Transfer_Operator> xfer4 = Teuchos::rcp( new Simple_EQ_Transfer(eq4ME, eq3ME) );

  // Flip the order of transfer registration; shouldn't matter
  pm5.add_transfer(xfer4);
  pm5.add_transfer(xfer3);

  pm5.register_complete(); // Trigger setup of coupled problem solver

  pm5.set_max_fixed_point_iters(25); // Allow sufficient iterations

  // The default fixed-point mode should be Jacobi; we will test and see
  pm5.integrate();

  const vector<double> & eq3_iterates = eq3ME->get_my_iterates();
  const vector<double> & eq4_iterates = eq4ME->get_my_iterates();

  TEST_EQUALITY(eq3_iterates.size(), eq4_iterates.size()+1);
  TEST_EQUALITY(jacobi_iters, eq3_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq3_iterates.size(); ++i )
    error += fabs(eq3_iterates[i]-jacobi_x1[i]);
  for( size_t i = 0; i < eq4_iterates.size(); ++i )
    error += fabs(eq4_iterates[i]-jacobi_x2[i]);
  error += 1.0; // shift needed for equality test
  TEST_FLOATING_EQUALITY( 1.0, error, 1.e-12);

  // Now solve using Seidel fixed-point 
  pm5.set_fixed_point_mode(Problem_Manager::seidel);
  pm5.integrate(); // also tests that problems get reinitialized

  TEST_EQUALITY(eq3_iterates.size(), eq4_iterates.size());
  TEST_INEQUALITY(seidel_iters, eq3_iterates.size());
  error = 0.0;
  for( size_t i = 0; i < eq3_iterates.size(); ++i )
    error += fabs(eq3_iterates[i]-seidel_x1[i]);
  for( size_t i = 0; i < eq4_iterates.size(); ++i )
    error += fabs(eq4_iterates[i]-seidel_x2[i]);
  error += 1.0; // shift needed for equality test
  TEST_INEQUALITY( 1.0, error);

  // Test that we make it to the end
  TEST_ASSERT(true);
}
