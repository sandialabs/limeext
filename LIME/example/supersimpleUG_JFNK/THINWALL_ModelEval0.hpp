// ********************************************************************************************
//
// THINWALL_ModelEval0.hpp
// This is the ss_thinwall LIME model evaluator used for a simple coupling with ss_con1d 
// 
// ********************************************************************************************
#ifndef LIME_EXAMPLE_THINWALL_MODELEVAL0_HPP
#define LIME_EXAMPLE_THINWALL_MODELEVAL0_HPP

#include "LIME_Model_Evaluator.hpp"
#include "LIME_Data_Transfer_Operator.hpp"

class THINWALL_ModelEval0 : public LIME::Model_Evaluator
{

  public:

    THINWALL_ModelEval0(LIME::Problem_Manager & pm, const string & name);   // C++ constructor

    virtual ~THINWALL_ModelEval0();                                         // C++ destructor

    // Tells the LIME problem manger that this code supports a stand-alone solve
    virtual bool supports_standalone_solve() const { return true; }

    // Tells the LIME problem manger that ss_con1d solves a transient problem LIME must direct
    virtual bool is_transient() const { return true; }

    // The implementation of the following functions must be defined here or in the .cpp file
    virtual void solve_standalone();
    virtual bool is_converged();
    virtual double get_time_step() const;
    virtual double get_max_time() const;
    virtual double get_current_time() const;
    virtual unsigned int get_max_steps() const { return 50; }
    virtual void update_time();

   // data required for transfer operators that will be implemented
   float* t;       // will transfer this array   to con1d from thinwall
   float& trad;    // will transfer this value from con1d   to thinwall

};


// -- ss_con1d to ss_thinwall transfer operator -------
//
class conduction_2_thinwall : public LIME::Data_Transfer_Operator {
public:
    conduction_2_thinwall(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

    // The implimentation of this function is defined in the c++ source file
    virtual bool perform_data_transfer() const;
};


// -- ss_thinwall to ss_con1d transfer operator -------
//
class thinwall_2_conduction : public LIME::Data_Transfer_Operator {
public:
    thinwall_2_conduction(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

    // The implimentation of this function is defined in the c++ source file
    virtual bool perform_data_transfer() const;
};

#endif
