// Yeckel, Pandy & Derby IJNME 2006.

// LIME headers
#include "LIME_Problem_Manager.hpp"

// For parsing command line
#include "Teuchos_CommandLineProcessor.hpp"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "ConvDiff_EQ.hpp"

using namespace std;

int main(int argc, char *argv[])
{
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

  Teuchos::CommandLineProcessor clp( false );

  // Default run-time options that can be changed from the command line
  bool                verbose         = false ;
  int                 numElements     = 20    ;
  int                 probSizeRatio   = 1     ;
  // Coupling parameters
  double              alpha           = 0.50  ;
  double              beta            = 0.40  ;
  // Physical parameters
  //double              radiation       = 5.67  ;
  double              radiation       = 0.00  ;

  clp.setOption( "verbose", "no-verbose", &verbose, "Verbosity on or off." );
  clp.setOption( "n", &numElements, "Number of elements" );
  clp.setOption( "nratio", &probSizeRatio, "Ratio of size of problem 2 to problem 1" );
  clp.setOption( "alpha", &alpha, "Interfacial coupling coefficient, alpha" );
  clp.setOption( "beta", &beta, "Interfacial coupling coefficient, beta" );
  clp.setOption( "radiation", &radiation, "Radiation source term coefficient, R" );

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) 
    return parse_return;

  numElements++; // convert #elements to #nodes

  // Create the Problem Manager
  LIME::Problem_Manager problemManager(comm, verbose);

  // Domain boundary temperaures
  double Tleft          = 0.98          ;
  double Tright         = 1.0           ;

  // Distinguish certain parameters needed for T1_analytic
  double peclet_1     	= 9.0           ;
  double peclet_2     	= 0.0           ;
  double kappa_1      	= 1.0           ;
  double kappa_2	= 0.1           ;

  double T1_analytic = ConvDiff_EQ::computeAnalyticInterfaceTemp( radiation, Tleft, Tright, kappa_2, peclet_1 );

  // Create Region 1 PDE
  string myName         = "Region_1"    ;
  double radiation_reg1 = 0.0           ;
  double xmin  		= 0.0           ;
  double xmax  		= 1.0           ;

  Teuchos::RCP<ConvDiff_EQ> reg1_eq = Teuchos::rcp(
      new ConvDiff_EQ (
        problemManager, 
        peclet_1,
        radiation_reg1,
        kappa_1,
        alpha,
        ConvDiff_EQ::RIGHT,
        xmin,
        xmax, 
        Tleft,
        T1_analytic,
        numElements, 
        myName  ) );

  problemManager.add_problem(reg1_eq);


  // Create Region 2 PDE
  myName 		        = "Region_2"    ;
  xmin  		        = 1.0           ;
  xmax  		        = 2.0           ;

  Teuchos::RCP<ConvDiff_EQ> reg2_eq = Teuchos::rcp(
      new ConvDiff_EQ (
        problemManager, 
        peclet_2,
        radiation,
        kappa_2,
        beta,
        ConvDiff_EQ::LEFT,
        xmin,
        xmax, 
        T1_analytic,
        Tright,
        probSizeRatio*numElements, 
        myName  ));

  problemManager.add_problem(reg2_eq);

  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp(new ConvDiff_Xfer(reg1_eq, reg2_eq));
  Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp(new ConvDiff_Xfer(reg2_eq, reg1_eq));
  problemManager.add_transfer(xfer1);
  problemManager.add_transfer(xfer2);

  problemManager.register_complete();

  // Solve the coupled problem

  // We need to override certain default settings that get applied to the default NOX solvers
  // constructed for each ConvDiff_EQ
  {
    Teuchos::RCP<Teuchos::ParameterList> default_nox_params = LIME::Problem_Manager::create_default_nox_params(comm);
    Teuchos::ParameterList& printParams = default_nox_params->sublist("Printing");
    printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			//NOX::Utils::OuterIterationStatusTest + 
			//NOX::Utils::InnerIteration +
                        //NOX::Utils::LinearSolverDetails +
			//NOX::Utils::Parameters + 
			//NOX::Utils::Details + 
			NOX::Utils::Warning);
    Teuchos::ParameterList& lsParams = default_nox_params->sublist("Direction").sublist("Newton").sublist("Linear Solver");
    lsParams.set("Max Iterations", 22);
    string nlFile1 = reg1_eq->createInArgs().modelEvalDescription() + "_nox.xml";
    string nlFile2 = reg2_eq->createInArgs().modelEvalDescription() + "_nox.xml";
    LIME::Problem_Manager::writeParameterList(nlFile1, *default_nox_params);
    LIME::Problem_Manager::writeParameterList(nlFile2, *default_nox_params);
  }

  problemManager.set_fixed_point_mode(LIME::Problem_Manager::jacobi);
  problemManager.set_max_fixed_point_iters(50);
  //problemManager.set_coupling_algorithm(LIME::Problem_Manager::jfnk);
  problemManager.integrate();

  //const Epetra_Vector & soln1 = *reg1_eq->get_x_state();
  //const Epetra_Vector & soln2 = *reg2_eq->get_x_state();
  //cout << setprecision(10) << "\n\tFinal Solution set:\n" << soln1 << "\n\n" << soln2 << endl;
  //cout << setprecision(10) << "\n\tFinal Solution set:\n" << soln2 << endl;
  //cout << setprecision(10) << "\n\tFinal Solution set:\n" << soln1 << endl;

  double chk_T1 = reg1_eq->get_interface_temp();
  double chk_T2 = reg1_eq->get_interface_temp();

  cout << setprecision(10) << "\nAnalytic T1 = " << T1_analytic << ", T1 from reg1 = " << chk_T1 << ", T1 from reg2 = " << chk_T2 << endl;

  double tol = 1.e-5;
  bool passed = ( (fabs(T1_analytic - chk_T1) < tol) && (fabs(T1_analytic - chk_T2) < tol) );

  if( passed )
    cout << "Test Passed!" << endl;
  else
    cout << "Test Failed!" << endl;

  // Finalize MPI
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0 ;
}

