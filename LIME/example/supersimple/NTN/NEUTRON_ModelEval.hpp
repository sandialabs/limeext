#ifndef LIME_EXAMPLE_NEUTRON_MODELEVAL_HPP
#define LIME_EXAMPLE_NEUTRON_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class NEUTRON_ModelEval : public LIME::problem_manager_api
{

  public:

    NEUTRON_ModelEval(const LIME::Problem_Manager & pm, const string & name);

    virtual ~NEUTRON_ModelEval();

    virtual bool supports_standalone_solve() const { return true; }
    virtual void solve_standalone();

};
#endif

