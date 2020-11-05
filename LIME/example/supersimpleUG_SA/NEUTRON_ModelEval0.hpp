// ******************************************************************************************** 
// 
// NEUTRON_ModelEval0.hpp 
// This is a LIME model evaluator header file for a standalone wrap of the ss_neutron0.f90 code 
// 
// ********************************************************************************************
#ifndef LIME_EXAMPLE_NEUTRON_MODELEVAL0_HPP
#define LIME_EXAMPLE_NEUTRON_MODELEVAL0_HPP

#include "LIME_Model_Evaluator.hpp"

class NEUTRON_ModelEval0 : public LIME::Model_Evaluator
{

  public:

    NEUTRON_ModelEval0(LIME::Problem_Manager & pm, const string & name);    // C++ constructor

    virtual ~NEUTRON_ModelEval0();                                          // C++ destructor
 
    // Tells the LIME Problem Manager that this code supports a stand-alone solve
    virtual bool supports_standalone_solve() const { return true; }

    // The actual implementation of "solve_standalone" must be defined in the c++ source file
    virtual void solve_standalone();

};
#endif
