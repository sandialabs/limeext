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
    Teuchos::RCP<LIME::Problem_Manager> problemManager = Teuchos::rcp( new LIME::Problem_Manager(Comm) );

    Teuchos::CommandLineProcessor & clp = problemManager->get_command_line_options();

    // Default run-time options that can be changed from the command line
    int           numNodes        = 20    ;
    string        outputDir       = "."   ;
    string        paramsXMLfile   = "";

    clp.setOption( "n", &numNodes, "Number of elements" );
    clp.setOption( "outputdir", &outputDir, "Directory to output mesh and results into. Default is \"./\"" );
    clp.setOption( "paramsxml", &paramsXMLfile, "Filename for XML used to setup exposed parameters." );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) 
      return parse_return;

    // Go through motions to create an output subdir, even if it already exists
    mkdir(outputDir.c_str(), 0777 );
    outputDir += "/";
    //goldDir   += "/";

    // Get the process ID and the total number of processors
    //int MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    numNodes++; // convert #elements to #nodes

    // The number of unknowns must be at least equal to the number of processors.
    if (numNodes < NumProc) 
      throw std::runtime_error("numNodes cannot be < number of processors");

    // Set the exposed params file
    if( "" == paramsXMLfile )
      throw std::runtime_error("It appears you forgot to specify a filename for XML used to setup exposed parameters.");
    problemManager->set_params_xml_file(paramsXMLfile);

    // Create the physics components of the coupled Brusselator problem
    Teuchos::RCP<Brusselator_Temperature> temperature_eq = 
      Teuchos::rcp(new Brusselator_Temperature(*problemManager, numNodes, outputDir, problemManager->get_verbosity()));
    Teuchos::RCP<Brusselator_Species> species_eq = 
      Teuchos::rcp(new Brusselator_Species(*problemManager, numNodes, outputDir, problemManager->get_verbosity()));

    temperature_eq->initialize();
    species_eq->initialize();
    problemManager->add_problem(temperature_eq);
    problemManager->add_problem(species_eq);

    Teuchos::RCP<LIME::Data_Transfer_Operator> xfer1 =
      Teuchos::rcp(new Brusselator_EQ_Xfer(temperature_eq, species_eq));
    Teuchos::RCP<LIME::Data_Transfer_Operator> xfer2 =
      Teuchos::rcp(new Brusselator_EQ_Xfer(species_eq, temperature_eq));
    problemManager->add_transfer(xfer1);
    problemManager->add_transfer(xfer2);

    problemManager->register_complete(); // Trigger setup of groups, solvers, etc.

    problemManager->output_status(std::cout);

    // Emulate what TriKota will likely do ....

    double tol = 1.e-8;
    int isim = 0;
    double dgdp = 0.0;
    double p_val = 0.001;
    double g_val = 1.0e12;
    double old_p_val = p_val;
    double old_g_val = 0.0;

    // Get a reference to the parameter entry we wish to manipulate
    Teuchos::RCP<const Epetra_Vector> p_init = problemManager->get_p_init(0);
    Teuchos::RCP<const Teuchos::Array<std::string> > p_names = problemManager->get_p_names(0);
    const string freeParamName("shared_Dcoeff"); // this must be consisteent with XML file
    int i = 0;
    for( ; i < p_init->MyLength(); ++i )
      if( freeParamName == (*p_names)[i] )
        break;
    if( i == p_init->MyLength() )
      throw std::runtime_error("Could not find free parameter \"" + freeParamName + "\".  Check your XML file.");

    // Setup the path for setting input model parameters via the ME interface 
    Teuchos::RCP<Epetra_Vector> working_p_vec = Teuchos::rcp(new Epetra_Vector(*p_init));
    double & adjustableParam = (*working_p_vec)[i];
    EpetraExt::ModelEvaluator::InArgs inArgs = problemManager->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = problemManager->createOutArgs();
    inArgs.set_p(0, working_p_vec);

    // The optimization solve loop that TriKota will drive
    while( fabs(g_val) > tol )
    {
      adjustableParam = p_val;
      problemManager->evalModel(inArgs, outArgs);

      // Update problem solution states to reflect any parameter changes (BCs in this case)
      temperature_eq->initializeSolution();
      species_eq->initializeSolution();
      // Update problemManager composite solution with changes to problem solution states
      problemManager->refresh_solution_state();


      // Fire off the transient solve
      problemManager->integrate();

      // Now deal with responses
      problemManager->compute_problem_responses();
      const map<string, std::pair<Teuchos::RCP<LIME::Model_Evaluator>, Teuchos::RCP<Epetra_Vector> > > & g_vals = problemManager->get_responses();
      if( 1 < g_vals.size() )
        throw std::logic_error("This example should have only a single problem response function.");

      const Teuchos::RCP<Epetra_Vector> & gvec = (g_vals.begin()->second).second;
      g_val = (*gvec)[0];
      cout << " Response Function Value = " << g_val << " for Dcoeff = " << p_val << endl;


      // Use a secant method to update param value
      if( isim > 1 )
      {
        dgdp = (g_val - old_g_val)/(p_val - old_p_val);
        old_p_val = p_val;
        p_val -= g_val/dgdp;
      }
      else
      {
        old_p_val = p_val;
        p_val *= 2.0;
      }

      old_g_val = g_val;

      isim++;
    }

    // Tolerance on p_val needs to be less strict than on g_val
    bool passed = ( (fabs(g_val) < tol) && (fabs(p_val - 1.96911e-3) < 100*tol) );

    if( passed )
      cout << "Test Passed!" << endl;
    else
      cout << "Test Failed!" << endl;

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

