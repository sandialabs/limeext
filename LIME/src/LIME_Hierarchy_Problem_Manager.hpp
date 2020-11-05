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
// LIME_Hierarchy_Problem_Manager.hpp
//
// Description: 
// 
// ********************************************************************************************

#ifndef LIME_HIERARCHY_PROBLEM_MANAGER_H
#define LIME_HIERARCHY_PROBLEM_MANAGER_H

#include <utility>

#include "LIME_Model_Evaluator.hpp"
#include "LIME_Problem_Manager.hpp"

//-----------------------------------------------------------------------------
// 

namespace LIME {

  class Hierarchical_Problem_Manager : public Problem_Manager,
                                       public LIME::Model_Evaluator
  {
    public:

      Hierarchical_Problem_Manager(LIME::Problem_Manager & pm, bool verbose=false);

      virtual bool supports_standalone_solve() const { return true; }
      virtual bool solve_standalone(double & dt);
      virtual bool is_converged();
      virtual void reinitialize();

      virtual ~Hierarchical_Problem_Manager() {} 

      virtual const Problem_Manager & get_top_manager() const;
      virtual Problem_Manager & get_top_manager();

      // Get num problems for a hierarchy level
      virtual void increment_num_problems(unsigned hier);

      // Get unique problem id
      virtual bool get_next_problem_id(unsigned hier, ID_TYPE & new_id) const;

      virtual ID_TYPE add_problem(Teuchos::RCP<LIME::Model_Evaluator> prob);

      virtual void copy_incoming_transfers(const LIME::Problem_Manager & top_pm);
  };

}

#endif
