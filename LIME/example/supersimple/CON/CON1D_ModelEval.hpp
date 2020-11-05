#ifndef LIME_EXAMPLE_CON1D_MODELEVAL_HPP
#define LIME_EXAMPLE_CON1D_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class CON1D_ModelEval : public LIME::problem_manager_api
{

  public:

    CON1D_ModelEval(const LIME::Problem_Manager & pm, const string & name);

    virtual ~CON1D_ModelEval();

    virtual bool supports_standalone_solve() const { return true; }
    virtual void solve_standalone();

};
#endif

