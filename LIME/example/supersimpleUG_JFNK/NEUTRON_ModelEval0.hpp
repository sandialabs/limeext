// ******************************************************************************************** 
// 
// NEUTRON_ModelEval0.hpp 
// This is the ss_neutron LIME model evaluator used for a simple coupling with ss_con1d
// 
// ********************************************************************************************
#ifndef LIME_EXAMPLE_NEUTRON_MODELEVAL1_HPP
#define LIME_EXAMPLE_NEUTRON_MODELEVAL1_HPP

#include "LIME_Model_Evaluator.hpp"
#include "LIME_Elimination_Module.hpp"
#include "LIME_Data_Transfer_Operator.hpp"

class NEUTRON_ModelEval0 : public LIME::Model_Evaluator, public Elimination_Module
{

  public:

    NEUTRON_ModelEval0(LIME::Problem_Manager & pm, const string & name);    // C++ constructor

    virtual ~NEUTRON_ModelEval0();                                          // C++ destructor


    // The implementation of the following functions are defined here in the .cpp file
    virtual bool perform_elimination();
    virtual void solve_standalone();

   // data required for transfer operators that will be implemented
    float *q, *t;

};


// -- ss_con1d to ss_neutron transfer operator -------
//
class conduction_2_neutronics : public LIME::Data_Transfer_Operator {
public:
    conduction_2_neutronics(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

    // The implimentation of this function is defined in the c++ source file
    virtual bool perform_data_transfer() const;
};



// -- ss_neutron to ss_con1d transfer operator -------
//
class neutronics_2_conduction : public LIME::Data_Transfer_Operator {
public:
    neutronics_2_conduction(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to);

    // The implimentation of this function is defined in the c++ source file
    virtual bool perform_data_transfer() const;
};

#endif
