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
// LIME_Problem_As_Operator.hpp
//
// Description:
// 
// ********************************************************************************************

#ifndef problem_as_operator_hpp
#define problem_as_operator_hpp

#include <iostream>
//#include <limits>

#include "Epetra_Operator.h" // base class
#include "LIME_Model_Evaluator.hpp" // wrapped class

//-----------------------------------------------------------------------------

namespace LIME {

class problem_as_operator : public Epetra_Operator
{
 public:

  problem_as_operator(LIME::Model_Evaluator* prob) :
    wrapped_problem_(prob) { }
  virtual ~problem_as_operator() {}

  // Inherited from Epetra_Operator 
  // ....  and allows us to handle the ApplyInverse method where we perform 
  // ....  physics-based preconditioning

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { 
      // The two incoming vectors are actually the same and so a copy needs
      // to be done to avoid inadvertant data corruption.
      Epetra_Vector x(*X(0));
      return wrapped_problem_->apply_preconditioner(x, *Y(0)) ? 0 : -1; 
    }

  // Boilerplate to fulfill pure virtual methods of Epetra_Operator
  virtual int SetUseTranspose(bool UseTranspose)
    { return -1; }

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { std::cout << "problem_as_operator::Apply. Should't happen." << std::endl; return -1; /* We don't support this. */}

  virtual double NormInf() const
    { return std::numeric_limits<double>::max(); }

  virtual const char * Label() const
    { string name = wrapped_problem_->createInArgs().modelEvalDescription(); return name.c_str(); }

  virtual bool UseTranspose() const
    { return false; }

  virtual bool HasNormInf() const
    { return false; }

  virtual const Epetra_Comm & Comm() const
    { return wrapped_problem_->get_x_map()->Comm(); }

  virtual const Epetra_Map & OperatorDomainMap() const
    {  return *wrapped_problem_->get_x_map(); }

  virtual const Epetra_Map & OperatorRangeMap() const
    { return *wrapped_problem_->get_f_map(); }


 protected:
  LIME::Model_Evaluator * wrapped_problem_;
};
 
}

#endif
