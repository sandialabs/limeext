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
// LIME_PrePostOperator.hpp
//
// Description:
// 
// ********************************************************************************************

#ifdef HAVE_LIME_JFNK

#ifndef USER_PREPOSTOPERATOR_H
#define USER_PREPOSTOPERATOR_H

#include "NOX_Abstract_PrePostOperator.H"

//-----------------------------------------------------------------------------

namespace LIME {

class PrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  PrePostOperator() : NOX::Abstract::PrePostOperator() { }

  ~PrePostOperator() { }

  void runPreIterate(const NOX::Solver::Generic& solver);
  void runPostIterate(const NOX::Solver::Generic& solver);
  void runPreSolve(const NOX::Solver::Generic& solver);
  void runPostSolve(const NOX::Solver::Generic& solver);

  void outputX(const string &, const NOX::Solver::Generic& solver);
  void outputF(const string &, const NOX::Solver::Generic& solver);
  void outputX_and_F(const string &, const NOX::Solver::Generic& solver);
  void outputX_and_dX(const string &, const NOX::Solver::Generic& solver);
};

} // namespace LIME

#endif
#endif
