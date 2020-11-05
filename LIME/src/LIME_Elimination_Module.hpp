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
// LIME_Elimination_Module.hpp
//
// Description: 
// 
// ********************************************************************************************

#ifndef ELIMINATION_MODULE_H
#define ELIMINATION_MODULE_H

#include "Teuchos_RefCountPtr.hpp"
#include "EpetraExt_ModelEvaluator.h"

//-----------------------------------------------------------------------------
// 

class Epetra_Vector;

class Elimination_Module
{
 public:

  Elimination_Module() { }

  virtual ~Elimination_Module() { }

  virtual Teuchos::RCP<Epetra_Vector> get_input_vector() const;
  virtual Teuchos::RCP<const Epetra_Vector> get_output_vector() const;
  virtual bool perform_elimination() = 0;

protected:

  Teuchos::RCP<Epetra_Vector> inputs_;
  Teuchos::RCP<Epetra_Vector> outputs_;
};

inline
Teuchos::RCP<Epetra_Vector> Elimination_Module::get_input_vector() const
{ return inputs_; }

inline
Teuchos::RCP<const Epetra_Vector> Elimination_Module::get_output_vector() const
{ return outputs_; }

#endif
