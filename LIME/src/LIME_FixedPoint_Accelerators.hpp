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
// LIME_FixedPoint_Accelerators.hpp
//
// Description: Users should inherit from this, FixedPoint_Accelerator, class and 
// implement the compute_update method.
// 
// ********************************************************************************************

#ifndef LIME_FIXEDPOINT_ACCELERATORS_H
#define LIME_FIXEDPOINT_ACCELERATORS_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Epetra_MultiVector.h"

namespace LIME {

  class FixedPoint_Accelerator
  {
    public:

      FixedPoint_Accelerator() {}

      virtual ~FixedPoint_Accelerator() {}

      virtual bool compute_update(Epetra_Vector &) = 0;

    private:

  };


// For now, we implement the Anderson Acceleration using Teuchos BLAS/LAPACK
// and associated matrix/vector data.  This implies a serial implementation
// with no support for parallel, eg mpi.  The main component needed to make this
// parallel is a distributed QR factorization algorithm.

  class Anderson_Acceleration : public FixedPoint_Accelerator
  {

    public:

      Anderson_Acceleration(int, bool verbose = false);

      virtual ~Anderson_Acceleration() {}

      virtual void initialize(const Epetra_Vector &);

      virtual bool compute_update(Epetra_Vector &);

    private:

      bool verbose_;

      int mix_dim_;
      int compute_count_;

      Teuchos::RCP<Epetra_MultiVector> Dmatrix_;
      Teuchos::RCP<Epetra_MultiVector> Kmatrix_;

      Teuchos::LAPACK<int,double> lapack_;

      Teuchos::SerialDenseVector<int,double> x_vec_;
      Teuchos::SerialDenseVector<int,double> x_prev_;
      Teuchos::SerialDenseVector<int,double> x_diff_;
      Teuchos::SerialDenseVector<int,double> x_diff_prev_;
      Teuchos::SerialDenseVector<int,double> alpha_;
      Teuchos::SerialDenseVector<int,double> x_err_;

  };
}

#endif
