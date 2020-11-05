#ifndef LIME_EXAMPLE_CON1D_MODELEVAL1_HPP
#define LIME_EXAMPLE_CON1D_MODELEVAL1_HPP

#include "LIME_Model_Evaluator.hpp"

class CON1D_ModelEval0 : public LIME::Model_Evaluator
{

  public:

    CON1D_ModelEval0(LIME::Problem_Manager & pm, const string & name);

    virtual ~CON1D_ModelEval0();

    virtual bool supports_standalone_solve() const { return true; }
    virtual void solve_standalone();
    virtual bool is_converged();

    virtual bool is_transient() const { return true; }
    virtual double get_time_step() const;
    virtual double get_max_time() const;
    virtual double get_current_time() const;
    virtual unsigned int get_max_steps() const { return 50; }
    virtual void update_time();

   // neutronics needs access to this data
   float* qs; // neutronics will write this data
   float* t;  // neutronics will read this data
   float& trbc; // thinwall will write this data

   // "correct" solutions to verify answers
   static const float gold_temps[];
   static const float gold_temps_w_thinwall[];
};
#endif

