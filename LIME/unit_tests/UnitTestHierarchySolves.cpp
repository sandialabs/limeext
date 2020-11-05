
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

#include <fixtures/FiveLinearEQs.hpp>

using namespace LIME::Five_Linear_Equations;


TEUCHOS_UNIT_TEST(UnitTestHierarchySolves, testUnit1)
{
  Epetra_SerialComm comm;

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<EQ_ModelEval> eq1ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_ONE, verbose);
  Teuchos::RCP<EQ_ModelEval> eq2ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_TWO, verbose);
  Teuchos::RCP<EQ_ModelEval> eq3ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_THREE, verbose);
  Teuchos::RCP<EQ_ModelEval> eq4ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_FOUR, verbose);
  Teuchos::RCP<EQ_ModelEval> eq5ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_FIVE, verbose);

  pm.add_problem(eq1ME);
  pm.add_problem(eq2ME);
  pm.add_problem(eq3ME);
  pm.add_problem(eq4ME);
  pm.add_problem(eq5ME);

  // Get all the needed transfers for this coupled system
  vector<Teuchos::RCP<EQ_ModelEval> > MEs;
  MEs.push_back(eq1ME);
  MEs.push_back(eq2ME);
  MEs.push_back(eq3ME);
  MEs.push_back(eq4ME);
  MEs.push_back(eq5ME);

  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > xfers = get_all_transfers(MEs);
  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> >::iterator iter = xfers.begin(),
                                                            iter_end = xfers.end();
  for( ; iter_end != iter; ++iter )
    pm.add_transfer(*iter);

  pm.register_complete(); // Trigger setup of coupled problem solver

  // Allow sufficient iterations for either Jacobi or Seidel modes
  pm.set_max_fixed_point_iters(50);

  // The default fixed-point mode should be Jacobi; we will test and see
  pm.integrate();

  cout << "Soln EQ 1 = " << eq1ME->get_my_x() << endl;
  cout << "Soln EQ 2 = " << eq2ME->get_my_x() << endl;
  cout << "Soln EQ 3 = " << eq3ME->get_my_x() << endl;
  cout << "Soln EQ 4 = " << eq4ME->get_my_x() << endl;
  cout << "Soln EQ 5 = " << eq5ME->get_my_x() << endl;

  // Allow sufficient iterations for either Jacobi or Seidel modes
  pm.set_fixed_point_mode(LIME::Problem_Manager::seidel);

  // The default fixed-point mode should be Jacobi; we will test and see
  pm.integrate();

  cout << "Soln EQ 1 = " << eq1ME->get_my_x() << endl;
  cout << "Soln EQ 2 = " << eq2ME->get_my_x() << endl;
  cout << "Soln EQ 3 = " << eq3ME->get_my_x() << endl;
  cout << "Soln EQ 4 = " << eq4ME->get_my_x() << endl;
  cout << "Soln EQ 5 = " << eq5ME->get_my_x() << endl;

  // Test that we make it to the end
  TEST_ASSERT(true);
}


TEUCHOS_UNIT_TEST(UnitTestHierarchySolves, testUnit2)
{
  Epetra_SerialComm comm;

  bool verbose = false;

  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm2 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(pm, verbose));
  Teuchos::RCP<EQ_ModelEval> eq4ME = EQ_ModelEval::create_instance(*pm2, EQ_ModelEval::EQ_FOUR, verbose);

  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm3 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(*pm2, verbose));
  Teuchos::RCP<EQ_ModelEval> eq1ME = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_ONE, verbose);
  Teuchos::RCP<EQ_ModelEval> eq2ME = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_TWO, verbose);
  Teuchos::RCP<EQ_ModelEval> eq3ME = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_THREE, verbose);

  //Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm4 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(pm, verbose));
  Teuchos::RCP<EQ_ModelEval> eq5ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_FIVE, verbose);


  // setup the lower hierarchy solver
  pm3->add_problem(eq1ME);
  pm3->add_problem(eq2ME);
  pm3->add_problem(eq3ME);
  pm3->register_complete();
  pm3->set_max_fixed_point_iters(10);

  pm2->add_problem(eq4ME);
  pm2->add_problem(pm3);
  pm2->register_complete();
  pm2->set_max_fixed_point_iters(20);

  //pm4->add_problem(pm3);
  //pm4->register_complete();
  //pm4->set_max_fixed_point_iters(20);

  pm.add_problem(pm2);
  pm.add_problem(eq5ME);

  // Get all the needed transfers for this coupled system
  vector<Teuchos::RCP<EQ_ModelEval> > MEs;
  MEs.push_back(eq1ME);
  MEs.push_back(eq2ME);
  MEs.push_back(eq3ME);
  MEs.push_back(eq4ME);
  MEs.push_back(eq5ME);
  cout << "EQ 1 has id = " << eq1ME->get_id() << endl;
  cout << "EQ 2 has id = " << eq2ME->get_id() << endl;
  cout << "EQ 3 has id = " << eq3ME->get_id() << endl;
  cout << "EQ 4 has id = " << eq4ME->get_id() << endl;
  cout << "EQ 5 has id = " << eq5ME->get_id() << endl;

  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > xfers = get_all_transfers(MEs, verbose);
  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> >::iterator iter = xfers.begin(),
                                                            iter_end = xfers.end();
  for( ; iter_end != iter; ++iter )
    pm.add_transfer(*iter);

  pm.register_complete();
  pm.set_max_fixed_point_iters(75);
  pm.set_fixed_point_mode(LIME::Problem_Manager::jacobi);

  pm.integrate();

  cout << "Soln EQ 1 = " << eq1ME->get_my_x() << endl;
  cout << "Soln EQ 2 = " << eq2ME->get_my_x() << endl;
  cout << "Soln EQ 3 = " << eq3ME->get_my_x() << endl;
  cout << "Soln EQ 4 = " << eq4ME->get_my_x() << endl;
  cout << "Soln EQ 5 = " << eq5ME->get_my_x() << endl;

  // Allow sufficient iterations for either Jacobi or Seidel modes
  pm.set_fixed_point_mode(LIME::Problem_Manager::seidel);
  pm2->set_fixed_point_mode(LIME::Problem_Manager::seidel);
  pm3->set_fixed_point_mode(LIME::Problem_Manager::seidel);
  //pm4->set_fixed_point_mode(LIME::Problem_Manager::seidel);

  // The default fixed-point mode should be Jacobi; we will test and see
  pm.integrate();

  cout << "Soln EQ 1 = " << eq1ME->get_my_x() << endl;
  cout << "Soln EQ 2 = " << eq2ME->get_my_x() << endl;
  cout << "Soln EQ 3 = " << eq3ME->get_my_x() << endl;
  cout << "Soln EQ 4 = " << eq4ME->get_my_x() << endl;
  cout << "Soln EQ 5 = " << eq5ME->get_my_x() << endl;

  // Test that we make it to the end
  TEST_ASSERT(true);
}
