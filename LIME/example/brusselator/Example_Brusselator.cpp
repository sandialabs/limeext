// 1D Finite Element Brusselator Example Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T    
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2   
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C    
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2   
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 *
 * This problem is examined with a variety of time integration schemes in:
 * "Studies on the Convergence of Various Time-Integration Schemes for the
 * Radiation-Diffusion Problem," Curtis C. Ober & John N. Shadid, in prep.
 *
 * In this example, only a 1st-order fully implicit (backward Euler)
 * time integration scheme is considered currently.
 *
 * Values for time step size and finite spatial extent are specified in
 * the constructor initialization list in Brusselator.C using
 * variables dt  and xmin,xmax, respectively.
 * The number of time steps to be taken is specified by variable
 * maxTimeSteps below.
 */

// LIME headers
#include "LIME_Problem_Manager.hpp" 

// For parsing command line
#include "Teuchos_CommandLineProcessor.hpp"

// For mkdir
#include <sys/stat.h>

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Brusselator_Species.hpp"              
#include "Brusselator_Temperature.hpp"              

using namespace std;

int main(int argc, char *argv[])
{

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm Comm;
#endif

    // Create the Problem Manager
    LIME::Problem_Manager problemManager(Comm);

    Teuchos::CommandLineProcessor & clp = problemManager.get_command_line_options();

    // Default run-time options that can be changed from the command line
    int           numNodes        = 20    ;
    string        outputDir       = "."   ;

    clp.setOption( "n", &numNodes, "Number of elements" );
    clp.setOption( "outputdir", &outputDir, "Directory to output mesh and results into. Default is \"./\"" );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) 
      return parse_return;

    // Go through motions to create an output subdir, even if it already exists
    mkdir(outputDir.c_str(), 0777 );
    outputDir += "/";

    // Get the process ID and the total number of processors
    //int MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    numNodes++; // convert #elements to #nodes

    // The number of unknowns must be at least equal to the number of processors.
    if (numNodes < NumProc) 
      throw ("numNodes cannot be < number of processors");

    // Create the physics components of the coupled Brusselator problem
    Teuchos::RCP<Brusselator_Temperature> temperature_eq = 
      Teuchos::rcp(new Brusselator_Temperature(problemManager, numNodes, outputDir, problemManager.get_verbosity()));
    Teuchos::RCP<Brusselator_Species> species_eq = 
      Teuchos::rcp(new Brusselator_Species(problemManager, numNodes, outputDir, problemManager.get_verbosity()));
    //Burgers  burgers   (Comm, NumGlobalNodes, "Burgers"     );

    // An interesting note: the order of solving each problem is based on the
    // order of adding.  For this decoupled problem, Species is linear
    // with respect to its variables, whereas Temperature is nonlinear wrt to its
    // variables.  The order of solution appears to strongly affect the rate
    // of convergence of the decoupled Brusselator.  Solving Temperature first
    // dramatically reduces the number of total iterations.

    // We could allow problemManager to call an initialize for each problem furing setup
    temperature_eq->initialize();
    species_eq->initialize();
    problemManager.add_problem(temperature_eq);
    problemManager.add_problem(species_eq);

    Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 = Teuchos::rcp(new Brusselator_EQ_Xfer(temperature_eq, species_eq));
    Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 = Teuchos::rcp(new Brusselator_EQ_Xfer(species_eq, temperature_eq));
    problemManager.add_transfer(xfer1);
    problemManager.add_transfer(xfer2);

    problemManager.register_complete(); // Trigger setup of groups, solvers, etc.

    problemManager.output_status(std::cout);

    // Fire off the transient solve
    problemManager.integrate();

    cout << "Test Passed!" << endl;

    // shutdown
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  }
  catch (exception& e) {
    cout << e.what() << endl;
    cout << "Test FAILED!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  }
  catch (...) {
    cout << "Error: caught unknown exception, exiting." << endl;
    cout << "Test FAILED!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  }


  return 0 ;
}

