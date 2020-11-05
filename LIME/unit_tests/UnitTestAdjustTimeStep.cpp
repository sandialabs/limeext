
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

namespace {
  class Mock_ModelEval : public LIME::Model_Evaluator
  {

    public:

      Mock_ModelEval(LIME::Problem_Manager & pm, const string & name) :
        LIME::Model_Evaluator(pm, name),
        times_solvestandalone_called_(0),
        max_time_(LIME::Model_Evaluator::get_max_time()),
        dt_(LIME::Model_Evaluator::get_time_step()),
        time_(LIME::Model_Evaluator::get_current_time()),
        max_steps_(LIME::Model_Evaluator::get_max_steps()),
        adjust_dt_enabled_(false)
    {}

      ~Mock_ModelEval() {}

      // Inherited from Model_Evaluator
      bool supports_standalone_solve() const { return true; }

      bool solve_standalone(double & dt)
      { 
        times_solvestandalone_called_++;
        if( adjust_dt_enabled_ && (time_ >= adjust_at_time_) )
        {
          dt = adjust_to_value_;
          dt_ = dt;
          disable_dt_adjust();
        }

        return true; 
      }
      // Need to return "false" in order for solve_standalone to get called
      // though this will trigger a warning from Problem_Manager
      bool is_converged() { return false; }

      bool is_transient() const { return true; }
      double get_max_time() const { return max_time_; }

      double get_time_step() const 
      { 
        return dt_; 
      }

      double get_current_time() const { return time_; }
      unsigned get_max_steps() const { return max_steps_; }

      void set_time_step(double dt) { dt_ = dt; }
      void update_time() 
      { 
        solved_times_.push_back(time_); 
        time_ += dt_;
      }

      // Methods to facilitate testing
      int get_times_called()
      { return times_solvestandalone_called_; }

      const vector<double> & get_solved_times()
      { return solved_times_; }

      void reset( double max_time, double dt, double time, int max_steps )
      {
        max_time_ = max_time;
        dt_ = dt;
        time_ = time;
        max_steps_ = max_steps;
        times_solvestandalone_called_ = 0;
        solved_times_.clear();
      }

      void enable_dt_adjust(double time, double new_dt)
      {
        adjust_dt_enabled_ = true;
        adjust_at_time_ = time;
        adjust_to_value_ = new_dt;
      }

      void disable_dt_adjust()
      { adjust_dt_enabled_ = false; }

    private:

      int times_solvestandalone_called_;
      vector<double> solved_times_;
      double max_time_;
      double dt_;
      double time_;
      int max_steps_;

      bool adjust_dt_enabled_;
      double adjust_at_time_;
      double adjust_to_value_;
  };
}


TEUCHOS_UNIT_TEST(UnitTestAdjustTimeStep, testUnit)
{
  Epetra_SerialComm comm;

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<Mock_ModelEval> mockME1 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 1"));
  Teuchos::RCP<Mock_ModelEval> mockME2 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 2"));
  Teuchos::RCP<Mock_ModelEval> mockME3 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 3"));

  pm.add_problem(mockME1);
  pm.add_problem(mockME2);
  pm.add_problem(mockME3);

  // Trigger setup of coupled problem solver
  pm.register_complete(); 

  // Test default settings, ie no solves attempted because max_steps = 0, etc.
  pm.integrate();

  TEST_EQUALITY_CONST( mockME1->get_times_called(), 0 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 0 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 0 );

  // Test 0 steps but with non-zero final time and time step
  mockME1->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 0 /* max_steps */ );
  mockME2->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );

  pm.integrate();

  TEST_EQUALITY_CONST( mockME1->get_times_called(), 0 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 0 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 0 );

  // Test stepping constrained to 5 steps
  mockME1->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 5 /* max_steps */ );
  mockME2->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );

  pm.integrate();

  TEST_EQUALITY_CONST( mockME1->get_times_called(), 5 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 5 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 5 );

  // Test synchronous stepping constrained to "10" steps (see comment below)
  mockME1->reset( 1.1 /* max_time */, 0.01 /* dt */, 1.0 /* time */, 100 /* max_steps */ );
  mockME2->reset( 2.0 /* max_time */, 0.01 /* dt */, 1.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 2.0 /* max_time */, 0.01 /* dt */, 1.0 /* time */, 100 /* max_steps */ );

  pm.integrate();

  /* Note that at present, time stepping happens
               while time <= t_max
     which causes a solve to be called for time = t_max
     and produces one more than expected.  This might need to 
     be revisited. - RWH 10/25/2011
   */
  TEST_EQUALITY_CONST( mockME1->get_times_called(), 11 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 11 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 11 );
  TEST_EQUALITY( mockME1->get_solved_times().size(), 11 );
  TEST_EQUALITY( mockME2->get_solved_times().size(), 11 );
  TEST_EQUALITY( mockME3->get_solved_times().size(), 11 );
  TEST_FLOATING_EQUALITY( mockME1->get_current_time(), 1.11, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME2->get_current_time(), 1.11, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME3->get_current_time(), 1.11, 1.e-12 );

  // Test adjusted stepping - verify that increasing dt during solve_standalone is disallowed
  mockME1->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME2->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 0.1 /* max_time */, 0.01 /* dt */, 0.0 /* time */, 100 /* max_steps */ );

  mockME1->enable_dt_adjust(0.03, 0.02);

  TEST_THROW(pm.integrate(), std::runtime_error);
  TEST_EQUALITY_CONST( mockME1->get_times_called(), 4 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 3 ); // called one less time because of immediate exit upon chage in dt
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 3 ); // called one less time because of immediate exit upon chage in dt

  mockME1->disable_dt_adjust();


  // Test adjusted stepping - test decrease in dt during solve_standalone
  mockME1->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME2->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );

  mockME2->enable_dt_adjust(0.3, 0.05);

  pm.integrate();
  // Times : 0, 0.1, 0.2, 0.3, 0.35, 0.40, 0.45, 0.50, 0.55
  TEST_EQUALITY_CONST( mockME1->get_times_called(),  9 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(),  9 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(),  8 ); // called one less time because of immediate exit upon chage in dt
  TEST_EQUALITY( mockME1->get_solved_times().size(), 8 );
  TEST_EQUALITY( mockME2->get_solved_times().size(), 8 );
  TEST_EQUALITY( mockME3->get_solved_times().size(), 8 );
  TEST_FLOATING_EQUALITY( mockME1->get_current_time(), 0.55, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME2->get_current_time(), 0.55, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME3->get_current_time(), 0.55, 1.e-12 );

  mockME2->disable_dt_adjust();


  // Test adjusted stepping - test decrease in dt during solve_standalone from multiple MEs
  mockME1->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME2->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );
  mockME3->reset( 0.5 /* max_time */, 0.1 /* dt */, 0.0 /* time */, 100 /* max_steps */ );

  mockME2->enable_dt_adjust(0.3, 0.05);
  mockME3->enable_dt_adjust(0.45, 0.01); // mockME3->solve_standalone gets called twice for time = 0.2
                                      // so the second adjustment to dt happens at time 0.45

  pm.integrate();
  // Times : 0, 0.1, 0.2, 0.30, 0.35, 0.40, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51
  TEST_EQUALITY_CONST( mockME1->get_times_called(), 14 );
  TEST_EQUALITY_CONST( mockME2->get_times_called(), 14 );
  TEST_EQUALITY_CONST( mockME3->get_times_called(), 13 );
  TEST_EQUALITY( mockME1->get_solved_times().size(), 12 );
  TEST_EQUALITY( mockME2->get_solved_times().size(), 12 );
  TEST_EQUALITY( mockME3->get_solved_times().size(), 12 );
  TEST_FLOATING_EQUALITY( mockME1->get_current_time(), 0.51, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME2->get_current_time(), 0.51, 1.e-12 );
  TEST_FLOATING_EQUALITY( mockME3->get_current_time(), 0.51, 1.e-12 );
  double gold_values_array[] = { 0.0, 0.1, 0.2, 0.30, 0.35, 0.40, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50 };
  vector<double> gold_values( gold_values_array, gold_values_array + sizeof(gold_values_array) / sizeof(double) );
  const vector<double> & test_times = mockME1->get_solved_times();
  TEST_COMPARE_FLOATING_ARRAYS( test_times, gold_values, 1.e-12 );


  // Test that we make it to the end
  TEST_ASSERT(true);
}
