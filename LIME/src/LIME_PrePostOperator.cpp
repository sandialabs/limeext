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
// LIME_PrePostOperator.cpp
//
// Description:
// 
// ********************************************************************************************

#include "LIME_config.hpp"

#ifdef HAVE_LIME_JFNK

//-----------------------------------------------------------------------------
// Include Directives

#include "LIME_PrePostOperator.hpp"
#include "NOX_Epetra.H"

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::runPreIterate(const NOX::Solver::Generic& solver)
{
  return outputX("runPreIterate", solver);
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::runPostIterate(const NOX::Solver::Generic& solver)
{
  outputX_and_F("runPostIterate", solver);
  return outputX_and_dX("runPostIterate", solver);
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::runPreSolve(const NOX::Solver::Generic& solver)
{
  return outputX("runPreSolve", solver);
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::runPostSolve(const NOX::Solver::Generic& solver)
{
  return outputX_and_F("runPostSolve", solver);
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::outputX(const string & label, const NOX::Solver::Generic& solver)
{
  const NOX::Epetra::Group * solnGrp = dynamic_cast<const NOX::Epetra::Group*>(&solver.getSolutionGroup());
  const Epetra_Vector & solution = dynamic_cast<const NOX::Epetra::Vector&>(solnGrp->getX()).getEpetraVector();

  cout << "\n\t-----------------------------------------------------";
  cout << "\n\tUserPrePostOperator::" << label << " : solution vector:";
  cout << "\n\t-----------------------------------------------------" << endl;
  for( int i = 0; i < solution.MyLength(); ++i )
    cout << std::setprecision(16) << "x[" << i << "] = " << solution[i] << endl;
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::outputF(const string & label, const NOX::Solver::Generic& solver)
{
  const NOX::Epetra::Group * solnGrp = dynamic_cast<const NOX::Epetra::Group*>(&solver.getSolutionGroup());
  const Epetra_Vector & residual = dynamic_cast<const NOX::Epetra::Vector&>(solnGrp->getF()).getEpetraVector();

  cout << "\n\t-----------------------------------------------------";
  cout << "\n\tUserPrePostOperator::" << label << " : residual vector:";
  cout << "\n\t-----------------------------------------------------" << endl;
  for( int i = 0; i < residual.MyLength(); ++i )
    cout << std::setprecision(16) << "F[" << i << "] = " << residual[i] << endl;
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::outputX_and_F(const string & label, const NOX::Solver::Generic& solver)
{
  const NOX::Epetra::Group * solnGrp = dynamic_cast<const NOX::Epetra::Group*>(&solver.getSolutionGroup());
  const Epetra_Vector & solution = dynamic_cast<const NOX::Epetra::Vector&>(solnGrp->getX()).getEpetraVector();
  const Epetra_Vector & residual = dynamic_cast<const NOX::Epetra::Vector&>(solnGrp->getF()).getEpetraVector();

  cout << "\n\t-----------------------------------------------------";
  cout << "\n\tUserPrePostOperator::" << label << " : solution and residual vectors:";
  cout << "\n\t-----------------------------------------------------" << endl;
  for( int i = 0; i < solution.MyLength(); ++i )
    cout << std::setprecision(16) << "x[" << i << "] = " << solution[i] << ",\tF[" << i << "] = " << residual[i] << endl;
}

//-----------------------------------------------------------------------------

void 
LIME::PrePostOperator::outputX_and_dX(const string & label, const NOX::Solver::Generic& solver)
{
  const NOX::Epetra::Group * oldSolnGrp = dynamic_cast<const NOX::Epetra::Group*>(&solver.getPreviousSolutionGroup());
  //const NOX::Epetra::Group * solnGrp = dynamic_cast<const NOX::Epetra::Group*>(&solver.getSolutionGroup());
  const Epetra_Vector & oldSolution = dynamic_cast<const NOX::Epetra::Vector&>(oldSolnGrp->getX()).getEpetraVector();
  const Epetra_Vector & deltaX   = dynamic_cast<const NOX::Epetra::Vector&>(oldSolnGrp->getNewton()).getEpetraVector();

  cout << "\n\t-----------------------------------------------------";
  cout << "\n\tUserPrePostOperator::" << label << " : old solution and deltaX vectors:";
  cout << "\n\t-----------------------------------------------------" << endl;
  for( int i = 0; i < oldSolution.MyLength(); ++i )
    cout << std::setprecision(16) << "x[" << i << "] = " << oldSolution[i] << ",\tdX[" << i << "] = " << deltaX[i] << endl;
}
#endif
