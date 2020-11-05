// ********************************************************************************************
//
// lime_mpd3.cpp
// This is a LIME driver for coupling the "super-simple" codes con1d and neutron
// using the following Model Evaluators: CON1D_ModelEval_w_Resid and NEUTRON_ModelEval_w_Elim
// 
// ********************************************************************************************
//
// "#include" directives for the various header files that are needed
//
// c++ 
#include <exception>
#include <iostream>
#include <string>

// con1d physics
#include "CON1D_ModelEval_w_Resid.hpp"

// neutron physics
#include "NEUTRON_ModelEval_w_Elim.hpp"

// LIME headers
#include <LIME_Problem_Manager.hpp>

// Trilinos Objects
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//
// define symbols to be used without qualifying prefix
//
using std::cout;
using std::endl;
using std::exception;
using std::string;

// --------------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int rc = 0;  // return code will be set to a non zero value on errors

  try {    // the "try" construct is used to allow for catching "exceptions" that might
           // unexpectedly occur.  They can occur at any point in the program's
           // call stack so our exception handling policy is to catch all thrown excptions
           // here in main, output a diagnostic message, and terminate the program cleanly.


    // 1  Create a pointer to an instance of an "Epetra_SerialComm" object, called "comm",  
    //    which is a specfic type of "Epetra_Comm" object.  (on the heap).  Note that Trilinos  
    //    must be constructed with either an MPI communicator or, if we are not using MPI, a  
    //    class that has the same interface as the MPI communicator but does nothing.    

#ifdef HAVE_MPI
          MPI_Init(&argc, &argv);
          Epetra_MpiComm * comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
          Epetra_SerialComm * comm = new Epetra_SerialComm;
#endif

    // 2  Create an instance of a "Problem_Manager" object called "pm"  (on the stack)
    //    Note: set verbosity on/off by defining the boolean "verbose" as true or false
    
          bool verbose = true;
          LIME::Problem_Manager pm(*comm, verbose);


    // 3  Create pointers to an instance of a "CON1D_ModelEval_w_Resid" object called "con1d" 
    //    and an instance of a "NEUTRON_ModelEval_w_Elim" object called "NEUTRON_ModelEval_w_Elim" 
    //    (on the heap)

          Teuchos::RCP<CON1D_ModelEval_w_Resid> con1d =
                   Teuchos::rcp(new CON1D_ModelEval_w_Resid(pm,"CON1DModelEval"));

          Teuchos::RCP<NEUTRON_ModelEval_w_Elim> neutron =
                   Teuchos::rcp(new NEUTRON_ModelEval_w_Elim(pm,"NEUTRONModelEval"));

    // 4  Register or "add" our physics problems "con1d" and "neutron", with
    //    the Problem manager. When doing this, we get back values for the 
    //    identifier "con1d_id" and "neutron_id" from the Problem Manager. 

          pm.add_problem(con1d);
          pm.add_problem(neutron);

    // 5  Create and setup two LIME data transfer operators that we will need
     
          LIME::Data_Transfer_Operator* p_n2c = new neutronics_2_conduction(neutron, con1d);
          LIME::Data_Transfer_Operator* p_c2n = new conduction_2_neutronics(con1d, neutron);

          Teuchos::RCP<LIME::Data_Transfer_Operator> n2c_op = Teuchos::rcp(p_n2c);
          Teuchos::RCP<LIME::Data_Transfer_Operator> c2n_op = Teuchos::rcp(p_c2n);

 
    // 6  Register or "add" our two transfer operators "c2n" and "n2c" with 
    //    the Problem manager

          pm.add_preelimination_transfer(c2n_op);
          pm.add_postelimination_transfer(n2c_op);

          pm.register_complete(); // Trigger setup of groups, solvers, etc.

          // After doing register_complete, we can manually set some options in pm
          //pm.set_coupling_algorithm(LIME::Problem_Manager::jfnk);

          // Optional way to directly set max_fixed_point_iterations and fixed point tolerance
          // This over-rides default values and/or values set in input file Problem_Manager_setup_fp.xml
          pm.set_max_fixed_point_iters(10);
          pm.set_fixed_point_tolerance(0.001);

    // 7  We are now ready to let the problem manager drive the coupled problem

           pm.integrate();

    // 8 Test that we converged to the "correct" answer
           float * temp = con1d->t;
           const float * gold_temps = con1d->gold_temps;
           float * qs = neutron->q;
           const float * gold_qs = neutron->gold_qs;
           float error = 0.0;
           for( int i = 0; i < 3; ++i )
           {
             error = std::max( error, abs((temp[i] - gold_temps[i])/gold_temps[i]) );
             error = std::max( error, abs((qs[i] - gold_qs[i])/gold_qs[i]) );
           }
           cout << "Max relative error = " << error << endl;
           if( 1.e-4 > error )
             cout << "Test Passed!" << endl;
  }
  catch (exception& e)
  {
    cout << e.what() << endl
         << "\n"
         << "Test FAILED!" << endl;
         rc = -1;
  }
  catch (...) {
    cout << "Error: caught unknown exception, exiting." << endl
         << "\n"
         << "Test FAILED!" << endl;
             rc = -2;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return rc;
}
