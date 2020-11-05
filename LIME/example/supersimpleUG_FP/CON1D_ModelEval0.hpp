// ******************************************************************************************** 
// 
// CON1D_ModelEval0.hpp 
// This is the ss_con1d LIME model evaluator used for a simple coupling with ss_neutron0
// or ss_thinwall
// 
// ********************************************************************************************
#ifndef LIME_EXAMPLE_CON1D_MODELEVAL0_HPP
#define LIME_EXAMPLE_CON1D_MODELEVAL0_HPP

#include "LIME_Model_Evaluator.hpp"

class CON1D_ModelEval0 : public LIME::Model_Evaluator
{

  public:

    CON1D_ModelEval0(LIME::Problem_Manager & pm, const string & name);   // C++ constructor

    virtual ~CON1D_ModelEval0();                                         // C++ destructor

    // Tells the LIME problem manger that this code supports a stand-alone solve
    virtual bool supports_standalone_solve() const { return true; }
    
    // Tells the LIME problem manger that ss_con1d solves a transient problem LIME must direct
    virtual bool is_transient() const { return true; }


    // The implementation of the following functions must be defined here or in the .cpp file
    virtual void solve_standalone();
    virtual bool is_converged();
    virtual double get_time_step() const;
    virtual double get_max_time() const;
    virtual double get_current_time() const;
    virtual unsigned int get_max_steps() const { return 50; }
    virtual void update_time();

   // data required for transfer operators that will be implemented
   float* qs; // will be transferred from ss_neutron0   to ss_con1d
   float* t;  // will be transferred   to ss_neutron0 from ss_con1d
   float& trbc; // thinwall will write this data

};
#endif

