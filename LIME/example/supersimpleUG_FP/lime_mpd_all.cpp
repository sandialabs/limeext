// ********************************************************************************************
//
// lime_mpd_all.cpp
// This is a LIME driver for coupling the "super-simple" codes ss_con1d, ss_neutron0, and 
// ss_ thinwall  using the following Model Evaluators: 
//             CON1D_ModelEval0, NEUTRON_ModelEval0, and THINWALL_ModelEval0
// 
// ********************************************************************************************
//
// "#include" directives for the various header files that are needed
//
// c++ 
#include <exception>
#include <iostream>
#include <string>

// con1d physics header file
#include "CON1D_ModelEval0.hpp"

// neutron physics
#include "NEUTRON_ModelEval0.hpp"

// thinwall physics header file
#include "THINWALL_ModelEval0.hpp"

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


    // 3  Create pointers to 
    //    an instance of a "CON1D_ModelEval0"    object called "con1d", and 
    //    an instance of a "NEUTRON_ModelEval0"  object called "neutron", and 
    //    an instance of a "THINWALL_ModelEval0" object called "thinwall" 
    //    (on the heap)
    //    Note: There are some Trilinos-specific implementation details here. For example,
    //          Teuchos is a Trilinos library and RCP is a reference counted pointer

          Teuchos::RCP<CON1D_ModelEval0> con1d =
                   Teuchos::rcp(new CON1D_ModelEval0(pm,"CON1DModelEval0"));

          Teuchos::RCP<NEUTRON_ModelEval0> neutron =
                   Teuchos::rcp(new NEUTRON_ModelEval0(pm,"NEUTRONModelEval0"));

          Teuchos::RCP<THINWALL_ModelEval0> thinwall =
                   Teuchos::rcp(new THINWALL_ModelEval0(pm,"THINWALLModelEval0"));

    // 4  Register or "add" our physics problems "con1d", "neutron", and "thinwall" with
    //    the Problem manager. When doing this, we get back values for the 
    //    identifier "con1d_id", "neutron_id" and "thinwall_id" from the Problem Manager. 

          pm.add_problem(con1d);
          pm.add_problem(neutron);
          pm.add_problem(thinwall);

    // 5  Create and setup the four LIME data transfer operators that we will need
     
          LIME::Data_Transfer_Operator* p_n2c = new neutronics_2_conduction(neutron, con1d);
          LIME::Data_Transfer_Operator* p_c2n = new conduction_2_neutronics(con1d, neutron);
          LIME::Data_Transfer_Operator* p_t2c = new thinwall_2_conduction(thinwall, con1d);
          LIME::Data_Transfer_Operator* p_c2t = new conduction_2_thinwall(con1d, thinwall);

          Teuchos::RCP<LIME::Data_Transfer_Operator> n2c_op = Teuchos::rcp(p_n2c);
          Teuchos::RCP<LIME::Data_Transfer_Operator> c2n_op = Teuchos::rcp(p_c2n);
          Teuchos::RCP<LIME::Data_Transfer_Operator> t2c_op = Teuchos::rcp(p_t2c);
          Teuchos::RCP<LIME::Data_Transfer_Operator> c2t_op = Teuchos::rcp(p_c2t);

 
    // 6  Register or "add" our four transfer operators with
    //    the Problem manager

          pm.add_preelimination_transfer(c2n_op);
          pm.add_preelimination_transfer(c2t_op);
          pm.add_postelimination_transfer(n2c_op);
          pm.add_postelimination_transfer(t2c_op);

          pm.register_complete(); // Trigger setup of groups, solvers, etc.


    // 7  We are now ready to let the problem manager drive the coupled problem

           pm.integrate();

    // 8 Test that we converged to the "correct" answer
    //       float * temp = con1d->t;
    //       const float * gold_temps = con1d->gold_temps_w_thinwall;
    //       const float & Twall = thinwall->t[0];
    //       const float & gold_w_temp = thinwall->gold_w_temp;
    //       float error = 0.0;
    //       for( int i = 0; i < 3; ++i )
    //         error = max( error, abs((temp[i] - gold_temps[i])/gold_temps[i]) );
    //       error = max( error, abs((Twall - gold_w_temp)/gold_w_temp) );
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
