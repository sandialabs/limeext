#ifndef LIME_EXAMPLE_THINWALL_MODELEVAL_HPP
#define LIME_EXAMPLE_THINWALL_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class THINWALL_ModelEval : public LIME::Model_Evaluator
{

  public:

    THINWALL_ModelEval(const LIME::Problem_Manager & pm, const string & name);

    virtual ~THINWALL_ModelEval();

    virtual bool supports_standalone_solve() const { return true; }
    virtual void solve_standalone();

};
#endif

