#ifndef LIME_EXAMPLE_NEUTRON__W_ELIM_MODELEVAL_HPP
#define LIME_EXAMPLE_NEUTRON__W_ELIM_MODELEVAL_HPP

#include "NEUTRON_ModelEval.hpp"
#include "LIME_Elimination_Module.hpp"

class NEUTRON_ModelEval_w_Elim : public NEUTRON_ModelEval,
      public Elimination_Module
{

  public:

    NEUTRON_ModelEval_w_Elim(LIME::Problem_Manager & pm, const string & name);

    virtual ~NEUTRON_ModelEval_w_Elim() {}

    // Disable standalone solve mode so we can serve as an elimination module
    virtual bool supports_standalone_solve() const { return false; }

    virtual bool perform_elimination();
};

#endif
