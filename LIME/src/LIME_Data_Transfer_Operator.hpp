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
// LIME_Data_Transfer_Operator.hpp
//
// Description: Users should inherit from this, Data_Transfer_Operator, class and 
// implement the perform_data_transfer method.  The class is constructed with 
// two arguments representing the physics collaborating in the transfer operation.
// The first constructor argument is the source of the transfer data and is
// conceptually const while the second constructor argument is the target of the
// transfer and is clearly non-const.  The constructor will perform a runtime
// check to ensure that the source and target are not the same.
//
// As there is no encapsulation in this class, it would be cleaner if we used a 
// struct and store references to the source and target problem manager.  If the
// transfer operators are stored in a list, we may need a default constructor that
// would require using pointers that can be default initialized rather than
// references which can not be default initialized.  An additional improvement
// would be to follow a naming convention where typedef's are not upper case.
// Upper-case is usually reserved for macros whereas aliases (typedef) are lower
// or camel case.
//
// 
// ********************************************************************************************

#ifndef DATA_TRANSFER_OPERATOR_H
#define DATA_TRANSFER_OPERATOR_H

#include <stdexcept>

//-----------------------------------------------------------------------------
// 

namespace LIME {

typedef unsigned ID_TYPE;

class Model_Evaluator;

class Data_Transfer_Operator
{
 public:

  //Data_Transfer_Operator(ID_TYPE s, ID_TYPE t) : source_id_(s), target_id_(t) 
  //{
  //  // ensure user never codes a data transfer to itself
  //  if (source_id_ == target_id_) {
  //    throw std::runtime_error("Error: User attempted to register a transfer operator with itself");
  //  }
  //}

  Data_Transfer_Operator(Teuchos::RCP<Model_Evaluator> s, Teuchos::RCP<Model_Evaluator> t) : source_(s), target_(t) 
  {
    // ensure user never codes a data transfer to itself
    if (source_ == target_) {
      throw std::runtime_error("Error: User attempted to register a transfer operator with itself");
    }
  }

  virtual ~Data_Transfer_Operator() {}

  virtual bool perform_data_transfer() const = 0;

  ID_TYPE source_id() const { return source_->get_id(); }
  ID_TYPE target_id() const { return target_->get_id(); }

  const Teuchos::RCP<Model_Evaluator> & source() const { return source_; }
  Teuchos::RCP<Model_Evaluator> & target() const { return target_; }

private:

  ID_TYPE source_id_, target_id_;
  Teuchos::RCP<Model_Evaluator> source_;
  mutable Teuchos::RCP<Model_Evaluator> target_;
};

}

#endif
