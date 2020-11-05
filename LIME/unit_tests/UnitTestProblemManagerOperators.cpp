
#ifdef HAVE_LIME_JFNK

// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <fixtures/ThreeNonlinEQsRes.hpp>

using namespace LIME::Three_Nonlinear_Equations_w_Resid;


TEUCHOS_UNIT_TEST(UnitTestProblemManagerOperators, testOperators1)
{
#ifdef HAVE_MPI
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<EQ_ModelEval> eq1ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_ONE, verbose);
  Teuchos::RCP<EQ_ModelEval> eq2ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_TWO, verbose);
  Teuchos::RCP<EQ_ModelEval> eq3ME = EQ_ModelEval::create_instance(pm, EQ_ModelEval::EQ_THREE, verbose);

  pm.add_problem(eq1ME);
  pm.add_problem(eq2ME);
  pm.add_problem(eq3ME);

  // Get all the needed transfers for this coupled system
  vector<Teuchos::RCP<EQ_ModelEval> > MEs;
  MEs.push_back(eq1ME);
  MEs.push_back(eq2ME);
  MEs.push_back(eq3ME);

  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > xfers = get_all_transfers(MEs, verbose);
  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> >::iterator iter = xfers.begin(),
                                                            iter_end = xfers.end();
  for( ; iter_end != iter; ++iter )
    pm.add_transfer(*iter);

  pm.set_coupling_algorithm(LIME::Problem_Manager::jfnk);
  // Trigger setup of objects
  pm.register_complete();

  // Set solution to known values
  Epetra_Vector x_state(*pm.get_x_init());
  TEST_ASSERT( MEs.size() == (size_t)x_state.MyLength() );
  for( size_t i = 0; i < MEs.size(); ++i )
    x_state[i] = MEs[i]->get_known_x();

  // Compute a numerical jacobian
  Teuchos::RCP<NOX::Epetra::Group> nox_group = pm.get_nox_group();
  Teuchos::RCP<const NOX::Epetra::LinearSystem> linsys = nox_group->getLinearSystem();
  Teuchos::RCP<const Epetra_Operator> jacOp = linsys->getJacobianOperator();
  TEST_ASSERT( x_state.Map().SameAs(jacOp->OperatorDomainMap()) );
  nox_group->setX(x_state);
  nox_group->computeJacobian();

  // Test correctness of numerical jacobian wrt gold coeffs at known solution
  const vector<vector<double> > & gold_coeffs = get_gold_jacobian_coeffs();
  Teuchos::RCP<const Epetra_RowMatrix> jacMatrix = LIME::get_matrix_from_operator(*jacOp);
  int numRows = jacMatrix->NumMyRows();
  TEST_EQUALITY(3, numRows);
  int maxEntries = jacMatrix->MaxNumEntries();
  TEST_EQUALITY(3, maxEntries);
  vector<double> values(maxEntries);
  vector<int> indices(maxEntries);
  int numEntries = 0;
  for( int row = 0; row < numRows; ++row )
  {
    jacMatrix->ExtractMyRowCopy(row, maxEntries, numEntries, &values[0], &indices[0]);
    TEST_EQUALITY(3, numEntries);
    for( int col = 0; col < numEntries; ++col )
      TEST_FLOATING_EQUALITY( values[col], gold_coeffs[row][col], 1.e-5);
  }

  // Test applyInverse of operator
  NOX::Epetra::Vector randVec(x_state);
  randVec.random();
  NOX::Epetra::Vector RHSVec(x_state);
  jacMatrix->Apply(randVec.getEpetraVector(), RHSVec.getEpetraVector());
  RHSVec.scale(1.234);
  Teuchos::ParameterList emptyList;
  emptyList.set("Tolerance", 1.e-8);
  NOX::Epetra::Vector testVec(x_state);
  nox_group->applyJacobianInverse(emptyList, RHSVec, testVec);
  testVec.update(1.234, randVec, -1.0);
  double norm = 1.0 + testVec.norm();
  TEST_FLOATING_EQUALITY(1.0, norm, 1.e-4);
  //emptyList.print(std::cout);

  // Test that we make it to the end
  TEST_ASSERT(true);
}
#endif
