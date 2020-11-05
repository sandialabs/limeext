#ifndef LIME_EXAMPLE_NEUTRON_MODELEVAL_HPP
#define LIME_EXAMPLE_NEUTRON_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"
#include "LIME_Data_Transfer_Operator.hpp"

class NEUTRON_ModelEval : public LIME::Model_Evaluator
{

  public:

    NEUTRON_ModelEval(LIME::Problem_Manager & pm, const string & name);

    virtual ~NEUTRON_ModelEval();

    virtual bool supports_standalone_solve() const { return true; }
    virtual void solve_standalone();
    virtual bool is_converged();

    float *q, *t;

    // "correct" solutions to verify answers
    static const float gold_qs[];
};

class conduction_2_neutronics : public LIME::Data_Transfer_Operator {
public:
    conduction_2_neutronics(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

  virtual bool perform_data_transfer() const;
};

class neutronics_2_conduction : public LIME::Data_Transfer_Operator {
public:
    neutronics_2_conduction(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

  virtual bool perform_data_transfer() const;
};

#endif
