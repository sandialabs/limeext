#include <algorithm>
#include <iostream>

#include "NEUTRON_ModelEval_w_Elim.hpp"

#include "LIME_Problem_Manager.hpp"

extern "C" {

#if defined(__CYGWIN__)
   extern float __neutron_mod_MOD_q[4];
   extern float __neutron_mod_MOD_t[4];
   extern float __neutron_mod_MOD_r_max;
#define VAR_Q     __neutron_mod_MOD_q
#define VAR_T     __neutron_mod_MOD_t
#define VAR_R_MAX __neutron_mod_MOD_r_max
#elif defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__CYGWIN__)
   extern float __neutron_mod__q[4];
   extern float __neutron_mod__t[4];
   extern float __neutron_mod__r_max;
#define VAR_Q __neutron_mod__q
#define VAR_T __neutron_mod__t
#define VAR_R_MAX __neutron_mod__r_max
#elif defined(__INTEL_COMPILER)
   extern float neutron_mod_mp_q_[4];
   extern float neutron_mod_mp_t_[4];
   extern float neutron_mod_mp_r_max;
#define VAR_Q neutron_mod_mp_q_
#define VAR_T neutron_mod_mp_t_
#define VAR_R_MAX neutron_mod_mp_r_max
#endif

   void setup_neutron_(); 
   void residual_neutron_(); 
   void solve_neutron_(); 
   void finish_neutron_(); 
}

//------------NEUTRON w Elim setup -----------------------------------------------------------

NEUTRON_ModelEval_w_Elim::NEUTRON_ModelEval_w_Elim(LIME::Problem_Manager & pm, const string & name) :
  NEUTRON_ModelEval(pm, name)
{
}

//------------NEUTRON w Elim eliminate internal variables ------------------------------------
 
bool NEUTRON_ModelEval_w_Elim::perform_elimination()
{
  solve_neutron_();
  return is_converged();
}

//--------------------------------------------------------------------------------------------
