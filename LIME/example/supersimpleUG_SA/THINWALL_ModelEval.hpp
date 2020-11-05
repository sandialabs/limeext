// **************************************************************************************
//
// THINWALL_ModelEval.hpp
// This is a LIME model evaluator header file for a standalone wrap of the super-simple thinwall code 
// *************************************************************************************** 
//
#ifndef LIME_EXAMPLE_THINWALL_MODELEVAL_HPP
#define LIME_EXAMPLE_THINWALL_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class THINWALL_ModelEval : public LIME::Model_Evaluator
{

  public:

    THINWALL_ModelEval(LIME::Problem_Manager & pm, const string & name);   // C++ constructor

    virtual ~THINWALL_ModelEval();                                         // C++ destructor

    // Tells the LIME problem manger that this code supports a stand-alone solve
    virtual bool supports_standalone_solve() const { return true; }

    // The actual implementation of "solve_standalone" must be defined in the c++ source file
    virtual void solve_standalone();

};
#endif
