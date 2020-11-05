// ********************************************************************************************
//
// thinwall_sad.cpp
// This is the LIME stand-alone driver for the super-simple thinwall code
// 
// ********************************************************************************************
//
// "#include" directives for the various header files that are needed
//
// c++ 
#include <exception>
#include <iostream>
#include <string>

// thinwall physics
#include "THINWALL_ModelEval.hpp"

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


    // 3  Create a pointer to an instance of a "THINWALL_ModelEval" object called "thinwall" 
    //    (on the heap)
    //    Note: There are some Trilinos-specific implementation details here. For example,
    //          Teuchos is a Trilinos library and RCP is a reference counted pointer

          Teuchos::RCP<THINWALL_ModelEval> thinwall =
                   Teuchos::rcp(new THINWALL_ModelEval(pm,"THINWALLModelEval"));


    // 4  Register or "add" our physics problem "thinwall", with
    //    the Problem manager. When doing this, we get back values for the 
    //    identifier "thinwall_id" from the Problem Manager. 

           int thinwall_id = pm.add_problem(thinwall);
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

