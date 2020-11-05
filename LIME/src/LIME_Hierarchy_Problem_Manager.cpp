// ********************************************************************************************
// LIME 1.0:  Lightweight Integrating Multiphysics Environment for coupling codes, Version 1.0
// Copyright (c) 2012, Sandia National Laboratories
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted 
// provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright notice, this list of 
//      conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright notice, this list of 
//      conditions and the following disclaimer in the documentation and/or other materials  
//      provided with the distribution.
//    * Neither the name of Sandia National Laboratories nor the names of its contributors may 
//      be used to endorse or promote products derived from this software without specific prior 
//      written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// (Note: This is a BSD 3-Clause License.  For info. see www.opensource.org/licenses/BSD-3-Clause )
// --------------------------------------------------------------------------------------------
//
// LIME_Hierarchy_Problem_Manager.cpp
//
// Description: 
// 
// ********************************************************************************************

//-----------------------------------------------------------------------------
// Include Directives

#include <iostream>
#include <vector>
#include <utility>

#include "LIME_Hierarchy_Problem_Manager.hpp"

using std::cout;
using std::endl;
using std::vector;

using LIME::Model_Evaluator;
using LIME::Problem_Manager;

//-----------------------------------------------------------------------------

namespace LIME {

Hierarchical_Problem_Manager::Hierarchical_Problem_Manager(Problem_Manager & pm, bool verbose) :
  Problem_Manager(pm.Comm(), verbose),
  Model_Evaluator(pm, "Hierarchical_Problem_Manager")
{
  my_hierarchy_ = pm.get_hierarchy() + 1;
}

//-----------------------------------------------------------------------------

const Problem_Manager &
Hierarchical_Problem_Manager::get_top_manager() const
{ 
  if( 0 == owning_manager_.get_hierarchy() )
    return owning_manager_;

  return dynamic_cast<const Hierarchical_Problem_Manager&>(owning_manager_).get_top_manager();
}

//-----------------------------------------------------------------------------

Problem_Manager &
Hierarchical_Problem_Manager::get_top_manager() 
{ 
  if( 0 == owning_manager_.get_hierarchy() )
    return owning_manager_;

  return dynamic_cast<Hierarchical_Problem_Manager&>(owning_manager_).get_top_manager();
}

//-----------------------------------------------------------------------------

void
Hierarchical_Problem_Manager::increment_num_problems(unsigned hier)
{
  return get_top_manager().increment_num_problems(my_hierarchy_);
}

//-----------------------------------------------------------------------------

bool
Hierarchical_Problem_Manager::get_next_problem_id(unsigned hier, ID_TYPE & value) const
{
  return get_top_manager().get_next_problem_id(hier, value);
}

//-----------------------------------------------------------------------------

bool 
Hierarchical_Problem_Manager::solve_standalone(double & dt) 
{ 
  if( m_verbose )
    cout << "Hierarchical_Problem_Manager::solve_standalone() ...." << endl; 

  return solve_nonlinear(dt);
}

//-----------------------------------------------------------------------------

bool 
Hierarchical_Problem_Manager::is_converged() 
{
  if( m_verbose )
    cout << "Hierarchical_Problem_Manager::is_converged() ....";

  bool is_conv = true;
  map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::const_iterator p_iter =  problems_.begin(),
    p_iter_end = problems_.end();
  for( ; p_iter_end != p_iter; ++p_iter )
    is_conv &= p_iter->second->is_converged();

  if( m_verbose )
    cout << is_conv << endl;

  return is_conv;
}

//-----------------------------------------------------------------------------

void 
Hierarchical_Problem_Manager::reinitialize() 
{
  //cout << "Hierarchical_Problem_Manager::reinitialize() ....";
  map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::const_iterator p_iter =  problems_.begin(),
    p_iter_end = problems_.end();
  for( ; p_iter_end != p_iter; ++p_iter )
    p_iter->second->reinitialize();
}

//-----------------------------------------------------------------------------

ID_TYPE 
Hierarchical_Problem_Manager::add_problem(Teuchos::RCP<Model_Evaluator> prob)
{
  ID_TYPE id = Problem_Manager::add_problem(prob);

  // We define children as only leaf problems, ie not Hierarchical_Problem_Manager's
  if( prob->num_my_children() )
    std::copy(prob->get_children().begin(), prob->get_children().end(), std::inserter(my_children_, my_children_.end()));
  else
    my_children_.insert(id);

  return id;
}

//-----------------------------------------------------------------------------

void 
Hierarchical_Problem_Manager::copy_incoming_transfers(const Problem_Manager & top_pm)
{
  typedef multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >::const_iterator MMciter;

  const multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> > & all_xfers = top_pm.get_all_registered_transfers();

  vector<ID_TYPE> descendants = owning_manager_.get_all_descendants(); // includes my children
  //cout << "Problem " << my_id_ << " has " << descendants.size() << " siblings and descendants :\n";
  for( size_t i = 0; i < descendants.size(); ++i )
  {
    if( my_children_.end() != my_children_.find(descendants[i]))
    {
      //cout << "\t" << descendants[i] << " - mine.  ";
      for( size_t j = 0; j < descendants.size(); ++j )
      {
        if( my_children_.end() == my_children_.find(descendants[j]))
        {
          //cout << " [found descdendat not mine " << descendants[j] << "] ";
          pair<ID_TYPE, ID_TYPE> id_key(descendants[j], descendants[i]);
          if( all_xfers.count(id_key) )
          {
            pair<MMciter, MMciter> xfers = all_xfers.equal_range(id_key);
            for( MMciter iter = xfers.first; iter != xfers.second; ++iter )
            {
              incoming_transfers_.insert(pair<ID_TYPE, Teuchos::RCP<Data_Transfer_Operator> >(iter->first.first, iter->second) );
              //cout << "Added transfer <" << iter->first.first << ", " << iter->first.second << ">.";
            }
          }
        }
      }
    }
    else
    {
      //cout << "\t" << descendants[i] << " - sibling's.  ";
    }
    //cout << endl;
  }

  map<ID_TYPE, Teuchos::RCP<Model_Evaluator> >::iterator p_iter =  problems_.begin(),
    p_iter_end = problems_.end();
  for( ; p_iter_end != p_iter; ++p_iter )
  {
    Model_Evaluator& problem  = *(p_iter->second);
    problem.copy_incoming_transfers(top_pm);
  }
}

} // namespace LIME
