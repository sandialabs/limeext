// *************************************************************************************** 
// CON1D_ModelEval.hpp 
// This is a LIME model evaluator header file for a standalone wrap of the ss_con1d.f90 code 
// *************************************************************************************** 
//
#ifndef LIME_EXAMPLE_CON1D_MODELEVAL_HPP
#define LIME_EXAMPLE_CON1D_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class CON1D_ModelEval : public LIME::Model_Evaluator
{

  public:

    CON1D_ModelEval(LIME::Problem_Manager & pm, const string & name);   // C++ constructor

    virtual ~CON1D_ModelEval();                                         // C++ destructor

    // Tells the LIME Problem Manager that this code supports a stand-alone solve
    virtual bool supports_standalone_solve() const { return true; }

    // The actual implementation of "solve_standalone" must be defined in the c++ source file
    virtual void solve_standalone();

};
#endif

