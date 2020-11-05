
// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

// The Mock DeCART/Star problem interface
#include <MockModelEvals.hpp>

using namespace LIME;

int
main(int argc, char *argv[])
{
  // Boilerplate for setting up the LIME::Problem_Manager (currently using Epetra objects)

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  bool verbose = true;

  // Top-level problem manager
  LIME::Problem_Manager top_pm(comm, verbose);

  // Create instances of the native Mock Apps we will wrap as Model_Evaluators
  Mock_DeCART_APP decart_app(   /* H_rxn = */ -0.8, 
                                /* alpha = */  0.95, 
                                /* beta = */   0.73, 
                                /* a = */      0.82, 
                                /* b = */      0.56);

  Mock_Star_APP star_app(       /* rho0 = */   1.23, 
                                /* gamma = */  0.33);

  // Level-1 Hierarchical problem manager for steady-state coupling of DeCART(Energy) with Star
  Teuchos::RCP<LIME::Hierarchical_Problem_Manager> hier1_pm = Teuchos::rcp(new LIME::Hierarchical_Problem_Manager(top_pm, verbose));

  // Instantiate the DeStar Model Evaluators
  Teuchos::RCP<Mock_DeCART_Energy_ME>    decart_energy_me    = Teuchos::rcp( new Mock_DeCART_Energy_ME   (*hier1_pm, decart_app, 0.0, verbose) );
  Teuchos::RCP<Mock_DeCART_Depletion_ME> decart_depletion_me = Teuchos::rcp( new Mock_DeCART_Depletion_ME( top_pm  , decart_app, 5.0, verbose) );
  Teuchos::RCP<Mock_Star_ME            > star_me             = Teuchos::rcp( new Mock_Star_ME            (*hier1_pm, star_app,   verbose) );

  // setup the hierarchy solver

  // for level-1 hierarchy, order of registration matters for Seidel-based Fixed-point
  // but does not matter for Jacobi-based
  hier1_pm->add_problem(decart_energy_me);
  hier1_pm->add_problem(star_me);
  hier1_pm->register_complete();
  hier1_pm->set_max_fixed_point_iters(20);
  //hier1_pm->set_fixed_point_mode(Problem_Manager::jacobi);
  hier1_pm->set_fixed_point_mode(Problem_Manager::seidel);

  // Again for top-level, order of registration matters for Seidel but not for Jacobi
  top_pm.add_problem(hier1_pm);
  top_pm.add_problem(decart_depletion_me);

  // Add our data transfers
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp( new LIME::Transfer_Mock_DeCART_To_Star(decart_energy_me, star_me, verbose) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp( new LIME::Transfer_Mock_Star_To_DeCART(star_me, decart_energy_me, verbose) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer3 = Teuchos::rcp( new LIME::Transfer_Mock_DeCART_Energy_To_Depletion(decart_energy_me, decart_depletion_me, verbose) );
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer4 = Teuchos::rcp( new LIME::Transfer_Mock_DeCART_Depletion_To_Energy(decart_depletion_me, decart_energy_me, verbose) );
  top_pm.add_transfer(xfer1);
  top_pm.add_transfer(xfer2);
  top_pm.add_transfer(xfer3);
  top_pm.add_transfer(xfer4);

  // Signal completion of registration (all problems and data transfers
  top_pm.register_complete();

  // Manually set any additional solver options (can also get from command-line and/or xml setup files -- sould settle on one!)
  top_pm.set_max_fixed_point_iters(50);
  top_pm.set_fixed_point_mode(Problem_Manager::seidel);

  // The default fixed-point mode should be Jacobi; we will test and see
  top_pm.integrate();

  // Need to add some checks for correctness...

  cout << "Test Passed!" << endl;
}
