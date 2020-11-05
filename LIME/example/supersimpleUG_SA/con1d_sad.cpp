// ********************************************************************************************
//
// con1d_sad.cpp
// This is the LIME stand-alone driver for the super-simple con1d code
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
#include "CON1D_ModelEval.hpp"

// LIME headers
#include <LIME_Problem_Manager.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

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

          Epetra_Comm* comm = new Epetra_SerialComm;

    // 2  Create an instance of a "Problem_Manager" object called "pm"  (on the stack)
    //    Note: set verbosity on/off by defining the boolean "verbose" as true or false
    
          bool verbose = false;
          LIME::Problem_Manager pm(*comm, verbose);


    // 3  Create a pointer to an instance of a "CON1D_ModelEval" object called "con1d" 
    //    (on the heap)
    //    Note: There are some Trilinos-specific implementation details here. For example,
    //          Teuchos is a Trilinos library and RCP is a reference counted pointer

          Teuchos::RCP<CON1D_ModelEval> con1d =
                   Teuchos::rcp(new CON1D_ModelEval(pm,"CON1DModelEval"));


    // 4  Register or "add" our physics problem "con1d", with
    //    the Problem manager. When doing this, we get back values for the 
    //    identifier "con1d_id" from the Problem Manager. 

           int con1d_id = pm.add_problem(con1d);
           pm.register_complete(); // Trigger setup of groups, solvers, etc.
           pm.output_status(cout);

    // 5  We are now ready to let the problem manager drive the coupled problem

           pm.integrate();

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

  return rc;
}

