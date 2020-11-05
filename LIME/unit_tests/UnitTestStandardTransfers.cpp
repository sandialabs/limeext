
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

class Mock_ModelEval : public LIME::Model_Evaluator
{

  public:

    Mock_ModelEval(LIME::Problem_Manager & pm, const string & name) :
      LIME::Model_Evaluator(pm, name)
  {}

    ~Mock_ModelEval() {}

    //virtual bool supports_standalone_solve() const { return true; }
};

class Mock_Transfer : public LIME::Data_Transfer_Operator 
{
  public:

    Mock_Transfer(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to) :
      LIME::Data_Transfer_Operator(from, to)
  {}

    bool perform_data_transfer() const
    { return false; }
};

TEUCHOS_UNIT_TEST(UnitTestStandardTransfers, testUnit)
{
  Epetra_SerialComm comm;

  bool verbose = false;
  LIME::Problem_Manager pm(comm, verbose);

  Teuchos::RCP<Mock_ModelEval> mockME1 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 1"));
  Teuchos::RCP<Mock_ModelEval> mockME2 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 2"));
  Teuchos::RCP<Mock_ModelEval> mockME3 = Teuchos::rcp(new Mock_ModelEval(pm, "Mock ME 3"));

  LIME::ID_TYPE mock_id1 = pm.add_problem(mockME1);
  LIME::ID_TYPE mock_id2 = pm.add_problem(mockME2);
  pm.add_problem(mockME3);

  // Test that we disallow self transfers
  TEST_THROW( Mock_Transfer(mockME1, mockME1), std::runtime_error);

  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp( new Mock_Transfer(mockME1, mockME2) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp( new Mock_Transfer(mockME2, mockME1) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer3 = Teuchos::rcp( new Mock_Transfer(mockME1, mockME3) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer4 = Teuchos::rcp( new Mock_Transfer(mockME2, mockME3) );

  pm.add_transfer(xfer1);
  pm.add_transfer(xfer2);
  pm.add_transfer(xfer3);
  pm.add_transfer(xfer4);

  pm.register_complete(); // Trigger setup of coupled problem solver

  // Test that we get the correct upstream transfers for each ME
  vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > upstream_xfers = mockME1->get_incoming_transfers();
  TEST_ASSERT((1 == upstream_xfers.size()));
  TEST_ASSERT((mock_id2 == upstream_xfers[0]->source_id()));
  //
  upstream_xfers = mockME2->get_incoming_transfers();
  TEST_ASSERT((1 == upstream_xfers.size()));
  TEST_ASSERT((mock_id1 == upstream_xfers[0]->source_id()));
  //
  upstream_xfers = mockME3->get_incoming_transfers();
  TEST_ASSERT((2 == upstream_xfers.size()));
  TEST_ASSERT((mock_id1 == upstream_xfers[0]->source_id()));
  TEST_ASSERT((mock_id2 == upstream_xfers[1]->source_id()));

  // Test that we make it to the end
  TEST_ASSERT(true);
}
