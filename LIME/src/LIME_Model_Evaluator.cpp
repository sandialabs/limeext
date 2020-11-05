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
// LIME_Model_Evaluator.cpp
//
// Description:
// 
// ********************************************************************************************

//-----------------------------------------------------------------------------
// Include Directives

#include <stdexcept>
#include <vector>

#include "LIME_Model_Evaluator.hpp"
#include "LIME_Problem_Manager.hpp"

using std::logic_error;
using std::runtime_error;
using std::vector;

using LIME::Problem_Manager;

//-----------------------------------------------------------------------------

namespace LIME {

Model_Evaluator::Model_Evaluator(Problem_Manager & pm, string name) :
  my_name_(name), 
  m_computed_residual(false), 
  owning_manager_(pm),
  my_children_()
{
  ID_TYPE new_id;
  if( !pm.get_next_problem_id(pm.get_hierarchy(), new_id) )
    throw runtime_error("Could not obtain a valid provisional problem id for this hierarchy level.");

  //string problemName = problem->createInArgs().modelEvalDescription();
  //problems_[new_id] = problem;
  //id_to_name_[new_id] = problemName;
  //name_to_id_[problemName] = new_id;
  my_id_ = new_id;
  pm.add_provisional_problem_id(new_id);
  pm.increment_num_problems(pm.get_hierarchy());
}

//-----------------------------------------------------------------------------

void 
Model_Evaluator::copy_incoming_transfers(const Problem_Manager & top_pm)
{
  if( 0 != num_my_children() )
    throw logic_error("copy_incoming_transfers should have been handled by Hierarchical_Problem_Manager.");

  typedef multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> >::const_iterator MMciter;

  const multimap<pair<ID_TYPE, ID_TYPE>, Teuchos::RCP<Data_Transfer_Operator> > & all_xfers = top_pm.get_all_registered_transfers();
  vector<ID_TYPE> descendants = owning_manager_.get_all_descendants(); // includes me
  //cout << "Problem " << my_id_ << " has " << descendants.size() << " siblings and descendants :\n";
  for( size_t i = 0; i < descendants.size(); ++ i)
  {
    if( my_id_ == descendants[i] )
    {
      //cout << "\t" << descendants[i] << " - me.  ";
    }
    else
    {
      //cout << "\t" << descendants[i] << " - sibling's.  ";
      pair<ID_TYPE, ID_TYPE> id_key(descendants[i], my_id_);
      if( all_xfers.count(id_key) )
      {
        pair<MMciter, MMciter> xfers = all_xfers.equal_range(id_key);
        for( MMciter iter = xfers.first; iter != xfers.second; ++iter )
        {
          incoming_transfers_.insert(pair<ID_TYPE, Teuchos::RCP<Data_Transfer_Operator> >(iter->first.first, iter->second) );
          //cout << "Added transfer <" << iter->first.first << ", " << my_id_ << ">.";
        }
      }
    }
    //cout << endl;
  }
}

//-----------------------------------------------------------------------------

vector<Teuchos::RCP<Data_Transfer_Operator> >
Model_Evaluator::get_incoming_transfers() const
{ 
  typedef multimap<ID_TYPE, Teuchos::RCP<Data_Transfer_Operator> >::const_iterator MMciter;

  vector<Teuchos::RCP<Data_Transfer_Operator> > xfers;
  for( MMciter iter = incoming_transfers_.begin(); iter != incoming_transfers_.end(); ++iter )
    xfers.push_back(iter->second);

  return xfers;
}

} // namespace LIME

