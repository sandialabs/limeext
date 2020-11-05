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
// LIME_FixedPoint_Accelerators.cpp
//
// Description: 
// 
// ********************************************************************************************

// Relevant headers
#include <LIME_FixedPoint_Accelerators.hpp>

// Trilinos Objects
#include <Epetra_Vector.h>

typedef double ScalarType;

using namespace LIME;

Anderson_Acceleration::Anderson_Acceleration(int mix_dim, bool verbose) :
  FixedPoint_Accelerator(),
  verbose_(verbose),
  mix_dim_(mix_dim)
{
}

void 
Anderson_Acceleration::initialize(const Epetra_Vector & vec)
{
  const Epetra_BlockMap & map = vec.Map();

  Dmatrix_ = Teuchos::rcp( new Epetra_MultiVector(map, mix_dim_) );
  Kmatrix_ = Teuchos::rcp( new Epetra_MultiVector(map, mix_dim_) );

  Kmatrix_->PutScalar(0.0);
  Dmatrix_->PutScalar(0.0);

  int xlen = vec.MyLength();

  Teuchos::LAPACK<int,double> lapack_;

  x_vec_.resize(xlen);
  x_prev_.resize(xlen);
  x_diff_.resize(xlen);
  x_diff_prev_.resize(xlen);
  alpha_.resize(mix_dim_);

  compute_count_ = 0;
}

bool 
Anderson_Acceleration::compute_update(Epetra_Vector & xVec)
{
  if( verbose_)
    cout << endl;

  int vec_index = compute_count_ % mix_dim_;
  int num_vecs = std::min(compute_count_+1, mix_dim_);
  int vec_dim = xVec.MyLength();

  Teuchos::SerialDenseVector<int,ScalarType> x_vec(Teuchos::Copy, (double*)&xVec[0], vec_dim);

  x_diff_ = x_vec;
  x_diff_ -= x_prev_;
  for( int j = 0; j < vec_dim; ++j )
    (*Kmatrix_)[vec_index][j] = x_vec(j);
  if( compute_count_ < 1 )
    x_diff_prev_ = x_diff_;
  else
  {
    int diff_vec_index = (compute_count_-1) % (mix_dim_-1);
    int num_diff_vecs = std::min(compute_count_, mix_dim_-1);

    for( int j = 0; j < vec_dim; ++j )
      (*Dmatrix_)[diff_vec_index][j] = x_diff_(j) - x_diff_prev_(j);
    x_diff_prev_ = x_diff_;

    if( verbose_ )
    {
      cout << "\nDmat :" << endl;
      Dmatrix_->Print(cout);
    }

    // Copy Epetra_MultiVector into SerialDenseMatrix
    int numCols = num_diff_vecs;
    Teuchos::SerialDenseMatrix<int,ScalarType> testMat(Teuchos::Copy, (*Dmatrix_)[0], vec_dim, vec_dim, numCols);

    int info = 0;
    int lwork = -1;
    std::vector<ScalarType> tau(std::min(testMat.numRows(), testMat.numCols()));
    std::vector<ScalarType> work(1);
    lapack_.GEQRF(testMat.numRows(), testMat.numCols(), testMat.values(), testMat.stride(), &tau[0], &work[0], lwork, &info);

    // Compute QR factorization of testMat
    lwork = (int)work[0];
    work.resize(lwork);
    lapack_.GEQRF(testMat.numRows(), testMat.numCols(), testMat.values(), testMat.stride(), &tau[0], &work[0], lwork, &info);

    if( verbose_ )
    {
      cout << "After call to GEQRF, testMat : " << endl;
      testMat.print(cout);
    }

    // Explicitly construct Q and R factors 
    // NOTE:  The upper triangular part of testMat is copied into R and testMat becomes Q.
    Teuchos::SerialDenseMatrix<int,ScalarType> Rtmp( testMat.numRows(), numCols );
    for(int ii = 0; ii < numCols; ++ii)
      for(int jj = ii; jj < numCols; ++jj)
        Rtmp(ii,jj) = testMat(ii,jj); 

    if( verbose_ )
    {
      cout << "\nR matrix :" << endl;
      Rtmp.print(cout);
    }

    // Update RHS 
    Teuchos::SerialDenseVector<int,ScalarType> rhs(x_diff_);
    int num_reflectors = (int)tau.size();
    lapack_.ORMQR('L', 'T', rhs.length(), 1, num_reflectors, testMat.values(), testMat.stride(), &tau[0], rhs.values(), rhs.length(), &work[0], lwork, &info);
    rhs.scale(-1.0);

    Teuchos::SerialDenseMatrix<int,ScalarType> working_mat(num_diff_vecs,num_diff_vecs);
    Teuchos::SerialDenseVector<int,ScalarType> working_rhs(num_diff_vecs);
    for( int row = 0; row < num_diff_vecs; ++row )
    {
      working_rhs(row) = rhs(row);
      for( int col = row; col < num_diff_vecs; ++col )
        working_mat(row,col) = Rtmp(row,col);
    }

    if( verbose_ )
    {
      cout << "\nnew RHS :" << endl;
      working_rhs.print(cout);
      cout << "\nnew LHS :" << endl;
      working_mat.print(cout);
    }

    lapack_.TRTRS('U', 'N', 'N', working_mat.numRows(), 1, working_mat.values(), working_mat.numRows(), working_rhs.values(), working_mat.numRows(), &info);

    if( verbose_ )
    {
      cout << "\nAlpha_bar from QR solve (i=" << compute_count_ << ") :" << endl;
      working_rhs.print(cout);
    }
    int first_diff_index = ( num_diff_vecs < mix_dim_-1 ? 0 : compute_count_ % (mix_dim_-1) );
    alpha_(0) = -working_rhs(first_diff_index);
    int prev_index = first_diff_index;
    for( int j = 1; j < num_diff_vecs; ++j )
    {
      int diff_vec_index = (first_diff_index+j) % (mix_dim_-1);
      alpha_(j) = working_rhs(prev_index) - working_rhs(diff_vec_index);
      prev_index = diff_vec_index;
    }
    alpha_(num_diff_vecs) = 1.0 + working_rhs(prev_index);

    if( verbose_ )
    {
      cout << "\nAlpha from QR solve (i=" << compute_count_ << ") :" << endl;
      alpha_.print(cout);

      cout << "\nKmat :" << endl;
      Kmatrix_->Print(cout);
    }

    x_vec.putScalar(0.0);
    int first_vec_index = ( num_vecs < mix_dim_ ? 0 : (compute_count_+1) % mix_dim_ );
    for( int k = 0; k < num_vecs; ++k )
    {
      int vec_index = (first_vec_index + k) % mix_dim_;
      for( int row = 0; row < vec_dim; ++row )
        x_vec(row) += alpha_(k)*(*Kmatrix_)[vec_index][row];
    }
  }
  if( verbose_ )
  {
    cout << "index = " << (compute_count_+1) << endl;
    cout << endl << "-----------------\nx_prev :";
    x_prev_.print(cout);
    cout << endl << "-----------------\nx_vec :";
    x_vec.print(cout);
    cout << endl << "-----------------\nx_diff :";
    x_diff_.print(cout);
    cout << endl << "-----------------\nx_diff_prev :";
    x_diff_prev_.print(cout);
  }

  x_prev_ = x_vec;
  for( int i = 0; i < vec_dim; ++i )
    xVec[i] = x_vec(i);

  compute_count_++;

  return true;
}
