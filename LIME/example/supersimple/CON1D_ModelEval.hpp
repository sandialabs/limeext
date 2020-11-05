#ifndef LIME_EXAMPLE_CON1D_MODELEVAL_HPP
#define LIME_EXAMPLE_CON1D_MODELEVAL_HPP

#include "LIME_Model_Evaluator.hpp"

class CON1D_ModelEval : public LIME::Model_Evaluator
{

  public:

    CON1D_ModelEval(LIME::Problem_Manager & pm, const string & name);

    virtual ~CON1D_ModelEval();

    virtual bool supports_standalone_solve() const { return true; }
//    virtual bool supports_residual() const { return true; }
    virtual void solve_standalone();
    virtual bool is_converged();
    virtual void residual_con1d();

    virtual bool is_transient() const { return true; }
    virtual double get_time_step() const;
    virtual double get_max_time() const;
    virtual double get_current_time() const;
//    virtual unsigned int get_max_steps() const { return 100000; }
    virtual unsigned int get_max_steps() const { return 50; }
    virtual void update_time();

   // neutronics needs access to this data
   float* qs; // neutronics will write this data
   float* t; // neutronics will read this data

   // "correct" solutions to verify answers
   static const float gold_temps[];
};
#endif

