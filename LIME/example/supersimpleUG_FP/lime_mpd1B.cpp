// ********************************************************************************************
//
// lime_mpd1B.cpp
// This is a LIME driver for fixed pt coupling of super-simple codes ss_con1d and ss_thinwall
// where convergence is determined using a global residual based convergence check.
// The following Model Evaluators are used (note the last 2 are used indirectly): 
//            CON1D_ModelEval0_w_Resid, and THINWALL_ModelEval0_w_Resid
//            CON1D_ModelEval0,         and THINWALL_ModelEval0  
// ********************************************************************************************
//
// "#include" directives for the various header files that are needed
//
// c++ 
#include <exception>
#include <iostream>
#include <string>

// con1d and thinwall physics header files
#include "CON1D_ModelEval0_w_Resid.hpp"
#include "THINWALL_ModelEval0_w_Resid.hpp"

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


    // 3  Create pointers to an instance of a "CON1D_ModelEval0_w_Resid" object called "con1d"
    //    and an instance of a "THINWALL_ModelEval0_w_Resid" object called "thinwall" 
    //    Note: There are some Trilinos-specific implementation details here. For example,
    //          Teuchos is a Trilinos library and RCP is a reference counted pointer

          Teuchos::RCP<CON1D_ModelEval0_w_Resid> con1d =
                   Teuchos::rcp(new CON1D_ModelEval0_w_Resid(pm,"CON1DModelEval0_w_Resid"));

          Teuchos::RCP<THINWALL_ModelEval0_w_Resid> thinwall =
                   Teuchos::rcp(new THINWALL_ModelEval0_w_Resid(pm,"THINWALLModelEval0_w_Resid"));

    // 4  Register or "add" our physics problems "con1d" and "thinwall" with
    //    the Problem manager. When doing this, we get back values for the 
    //    identifier "con1d_id" and "thinwall_id" from the Problem Manager. 

          pm.add_problem(con1d);
          pm.add_problem(thinwall);

    // 5  Create and setup the two LIME data transfer operators that we will need
     
          LIME::Data_Transfer_Operator* p_t2c = new thinwall_2_conduction(thinwall, con1d);
          LIME::Data_Transfer_Operator* p_c2t = new conduction_2_thinwall(con1d, thinwall);

          Teuchos::RCP<LIME::Data_Transfer_Operator> t2c_op = Teuchos::rcp(p_t2c);
          Teuchos::RCP<LIME::Data_Transfer_Operator> c2t_op = Teuchos::rcp(p_c2t);

 
    // 6  Register or "add" our two transfer operators "t2c_op" and "c2t_op" with
    //    the Problem manager

          pm.add_preelimination_transfer(c2t_op);
          pm.add_postelimination_transfer(t2c_op);

          pm.register_complete(); // Trigger setup of groups, solvers, etc.


    // 7  We are now ready to let the problem manager drive the coupled problem

           pm.integrate();

    // 8 Test that we converged to the "correct" answer
    //   We will use known "gold" values for comparison.
    //       float * temp = con1d->t;
    //       const float gold_temps[] = {500.837, 473.952, 387.25};
    //       const float & Twall = thinwall->t[0];
    //       const float & gold_w_temp = 332.372;
    //       float error = 0.0;
    //       for( int i = 0; i < 3; ++i )
    //       {
    //         cout << "Temp[" << i << "] = " << temp[i] << ",\tgold_temp = " << gold_temps[i] << endl;
    //         error = max( error, abs((temp[i] - gold_temps[i])/gold_temps[i]) );
    //       }
    //       error = max( error, abs((Twall - gold_w_temp)/gold_w_temp) );
    //       cout << "Twall = " << Twall << ",\tgold_temp_wall = " << gold_w_temp << endl;
    //       cout << "Max relative error = " << error << endl;
    //       if( 1.e-4 > error )
    //         cout << "Test Passed!" << endl;
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
