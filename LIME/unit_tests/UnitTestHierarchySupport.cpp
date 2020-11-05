
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

#include "fixtures/FiveLinearEQs.hpp"

// Trilinos Objects
#include <Epetra_SerialComm.h>

using namespace LIME::Five_Linear_Equations;

TEUCHOS_UNIT_TEST(UnitTestHierarchySupport, testUnit1)
{
  Epetra_SerialComm comm;

  bool verbose = false;

  // Create the top-level Manager
  LIME::Problem_Manager pm(comm, verbose);

  // Create some lower-level hierarchy managers
  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm2 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(pm, verbose));
  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm3 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(*pm2, verbose));
  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm4 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(*pm3, verbose));
  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> pm5 = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(*pm4, verbose));

  // Test each hierarchy label
  TEST_EQUALITY(0, pm.get_hierarchy());
  TEST_EQUALITY(1, pm2->get_hierarchy());
  TEST_EQUALITY(2, pm3->get_hierarchy());
  TEST_EQUALITY(3, pm4->get_hierarchy());
  TEST_EQUALITY(4, pm5->get_hierarchy());

  // Test ability to acquire top-level manager
  TEST_EQUALITY(0, pm2->get_top_manager().get_hierarchy());
  TEST_EQUALITY(0, pm3->get_top_manager().get_hierarchy());
  TEST_EQUALITY(0, pm4->get_top_manager().get_hierarchy());
  TEST_EQUALITY(0, pm5->get_top_manager().get_hierarchy());

  // Test assignment of unique problem ids at construction
  Teuchos::RCP<EQ_ModelEval> eq1_1 = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_ONE);
  Teuchos::RCP<EQ_ModelEval> eq1_2 = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_ONE);
  Teuchos::RCP<EQ_ModelEval> eq1_3 = EQ_ModelEval::create_instance(*pm3, EQ_ModelEval::EQ_ONE);
  TEST_INEQUALITY(eq1_1->get_id(), eq1_2->get_id());
  TEST_INEQUALITY(eq1_1->get_id(), eq1_3->get_id());
  TEST_INEQUALITY(eq1_2->get_id(), eq1_3->get_id());


  // Test that we make it to the end
  TEST_ASSERT(true);
}
