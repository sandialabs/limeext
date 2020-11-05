
#include "Brusselator_Temperature.hpp"

#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.hpp"

#include "LIME_Problem_Manager.hpp"

Brusselator_Temperature::Brusselator_Temperature(LIME::Problem_Manager & pm, int numGlobalNodes, const string & outputDir, bool verbose) :
  Brusselator_Equation(pm, "Brusselator_Temperature", numGlobalNodes, outputDir, verbose)
{
  m_g_map = Teuchos::rcp(new Epetra_LocalMap(numResponses, 0, pm.Comm()));
  m_g_vec = Teuchos::rcp(new Epetra_Vector(*m_g_map), false);

  EpetraExt::ModelEvaluator::OutArgs outArgs = createOutArgs();
  const ModelEvaluator::Evaluation<Epetra_Vector> g_vec(m_g_vec, EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT);
  outArgs.set_g(0, g_vec);

  pm.post_problem_response("BrusselatorResponse", Teuchos::rcp(this, false), m_g_vec);
}

//-----------------------------------------------------------------------------

Brusselator_Temperature::~Brusselator_Temperature()
{
}

//-----------------------------------------------------------------------------

EpetraExt::ModelEvaluator::OutArgs
Brusselator_Temperature::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription(my_name_);

  // For now, just support residual fills
  outArgs.setSupports(OUT_ARG_f, true);

  outArgs.set_Np_Ng(3, 1);

  return outArgs;
}

//-----------------------------------------------------------------------------

void 
Brusselator_Temperature::initializeSolution()
{
  Epetra_Vector & soln = *m_solution;
  Epetra_Vector & x = *m_mesh_x;

  // Here is a constant initial condition
  //soln.PutScalar(0.6);

  // ... or

  // Here we do a sinusoidal perturbation of the unstable
  // steady state.
  double pi = 4.*atan(1.0);

  double alpha  = (*m_p_vec)[0];
  for (int i=0; i<x.MyLength(); i++)
    soln[i] = alpha + 1.e-1*sin(1.0*pi*x[i]);

  *m_old_solution = soln;
} 

//-----------------------------------------------------------------------------

void 
Brusselator_Temperature::finalize_solution(const Epetra_Vector & x)
{
  Brusselator_Equation::finalize_solution(x);

  double max_temp = 0.0;
  x.MaxValue(&max_temp);

  if( max_temp > m_max_temp )
    m_max_temp = max_temp;

  cout << " \t----  Max Temp = " << m_max_temp << endl;
}

//-----------------------------------------------------------------------------

void 
Brusselator_Temperature::reinitialize()
{
  Brusselator_Equation::reinitialize();
  m_max_temp = 0.0;
}

//-----------------------------------------------------------------------------

// Matrix and Residual Fills

bool 
Brusselator_Temperature::evaluate(
    NOX::Epetra::Interface::Required::FillType fillType,
    const Epetra_Vector & soln, 
    Epetra_Vector & rhs) const 
{
  bool fillRes = false;
  bool fillJac = false;

  // We currently make no provision for combined Residual and Jacobian fills
  // i.e. it is either or

  if( NOX::Epetra::Interface::Required::Jac == fillType )
    fillJac = true;
   else
   {
    fillRes = true;
   }

  int numDep = 1;
  if( 1 != numDep )
    throw "Brusselator_Temperature::evaluate : number of auxiliary data vecors != 1.";

  Teuchos::RCP<const Epetra_Map> OverlapMap = get_x_map();

  // TODO - Create the overlapped solution and position vectors 
  Epetra_Vector u(*OverlapMap);
  Epetra_Vector uold(*OverlapMap);
  vector<Epetra_Vector*> dep(numDep);
  for( int i = 0; i<numDep; i++)
    dep[i] = new Epetra_Vector(*OverlapMap);
  Epetra_Vector xvec(*OverlapMap);

  // TODO - Export Solution to Overlap vector
  // If the vector to be used in the fill is already in the Overlap form,
  // we simply need to map on-processor from column-space indices to
  // OverlapMap indices. Note that the old solution is simply fixed data that
  // needs to be sent to an OverlapMap (ghosted) vector.  The conditional
  // treatment for the current soution vector arises from use of
  // FD coloring in parallel.
  uold = *m_old_solution;
  //uold.Import(*m_old_solution, *Importer, Insert);
  //for( int i = 0; i<numDep; i++ )
  //  dep[i]->Import(*( (*(depSolutions.find(depProblems[i]))).second ), *Importer, Insert);
  *(dep[0]) = *(m_aux_data);

  //xvec.Import(*xptr, *Importer, Insert);
  xvec = *m_mesh_x;

  //if( NOX::Epetra::Interface::Required::FD_Res == fillType )
  //  // Overlap vector for solution received from FD coloring, so simply reorder
  //  // on processor
  //  u.Export(*soln, *ColumnToOverlapImporter, Insert);
  //else // Communication to Overlap vector is needed
  //  u.Import(*soln, *Importer, Insert);
  u = soln;

  //int OverlapNumMyNodes = OverlapMap->NumMyElements();
  int numOverlapNodes = m_x_epetra_map->NumMyElements();

  int row;
  double alpha  = (*m_p_vec)[0];
  double beta   = (*m_p_vec)[1];
  double Dcoeff = (*m_p_vec)[2];
  //cout << "In Brusselator_Temperature::evaluate, using Dcoeff = " << Dcoeff << endl;
  //double jac;
  double xx[2];
  double uu[2]; 
  double uuold[2];
  vector<double*> ddep(numDep);
  for( int i = 0; i<numDep; i++)
    ddep[i] = new double[2];
  Basis basis;

  double convection = 0.0;

  // Zero out the objects that will be filled
  //if ( fillJac ) A->PutScalar(0.0);
  if ( fillRes ) rhs.PutScalar(0.0);

  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < numOverlapNodes-1; ne++) 
  {
    
    // Loop Over Gauss Points
    for(int gp=0; gp < 2; gp++) 
    {
      // Get the solution and coordinates at the nodes 
      xx[0]=xvec[ne];
      xx[1]=xvec[ne+1];
      uu[0] = u[ne];
      uu[1] = u[ne+1];
      uuold[0] = uold[ne];
      uuold[1] = uold[ne+1];
      for( int i = 0; i<numDep; i++ ) 
      {
        ddep[i][0] = (*dep[i])[ne];
        ddep[i][1] = (*dep[i])[ne+1];
      }
      // Calculate the basis function at the gauss point
      basis.getBasis(gp, xx, uu, uuold, ddep);

      // Loop over Nodes in Element
      for (int i=0; i< 2; i++) 
      {
        row=OverlapMap->GID(ne+i);
        if( m_x_epetra_map->MyGID(row) ) 
        {
          if ( fillRes ) 
          {
            convection = 0.0;
            //if( useConvection )
            //  convection = basis.ddep[id_vel]*basis.duu/basis.dx;

            (rhs)[m_x_epetra_map->LID(OverlapMap->GID(ne+i))]+=
              + basis.wt*basis.dx
              * ((basis.uu - basis.uuold)/m_dt * basis.phi[i] 
              + convection * basis.phi[i] 
              + (1.0/(basis.dx*basis.dx))*Dcoeff*basis.duu*basis.dphide[i]
              + basis.phi[i] * ( -alpha + (beta+1.0)*basis.uu
              - basis.uu*basis.uu*basis.ddep[0]) );
          }
        }
        // Loop over Trial Functions
        //if ( fillJac ) 
        //{
        //  for( int j = 0; j < 2; ++j ) 
        //  {
        //    if (m_x_epetra_map->MyGID(row)) 
        //    {
        //      column=OverlapMap->GID(ne+j);
        //      jac=basis.wt*basis.dx*(
        //              basis.phi[j]/m_dt*basis.phi[i] 
        //              +(1.0/(basis.dx*basis.dx))*Dcoeff*basis.dphide[j]*
        //                                                basis.dphide[i]
        //              + basis.phi[i] * ( (beta+1.0)*basis.phi[j]
        //              - 2.0*basis.uu*basis.phi[j]*basis.ddep[id_spec]) );  
        //      A->SumIntoGlobalValues(row, 1, &jac, &column);
        //    }
        //  }
        //}
      }
    }
  } 

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(0)=1
  //if (MyPID==0) {
    if ( fillRes )
      (rhs)[0]= soln[0] - alpha;
    //if ( fillJac ) {
    //  int column=0;
    //  double jac=1.0;
    //  A->ReplaceGlobalValues(0, 1, &jac, &column);
    //  column=1;
    //  jac=0.0;
    //  A->ReplaceGlobalValues(0, 1, &jac, &column);
    //}
  //}
  // U(1)=1
  if ( m_x_epetra_map->LID(m_x_epetra_map->MaxAllGID()) >= 0 ) {
    int lastDof = m_x_epetra_map->LID(m_x_epetra_map->MaxAllGID());
    if ( fillRes )
      (rhs)[lastDof] = soln[lastDof] - alpha;
    //if ( fillJac ) {
    //  int row=m_x_epetra_map->MaxAllGID();
    //  int column = row;
    //  double jac = 1.0;
    //  A->ReplaceGlobalValues(row, 1, &jac, &column);
    //  jac=0.0;
    //  column--;
    //  A->ReplaceGlobalValues(row, 1, &jac, &column);
    //}
  }

  // Sync up processors to be safe
  //Comm->Barrier();
 
  //A->FillComplete();

  //A->Print(cout);

  //if( fillRes )
  //  cout << "For residual fill :" << endl << rhs << endl;

  //if( fillJac ) {
  //  cout << "For jacobian fill :" << endl;
  //  A->Print(cout);
  //}

  // Cleanup
  for( int i = 0; i<numDep; i++)
  {
    delete [] ddep[i];
    delete    dep[i];
  }

  return true;
}

//-----------------------------------------------------------------------------

/* A representative response function

   g = 1.5 - T_max over all time steps

*/

void 
Brusselator_Temperature::compute_problem_response(const Epetra_Vector & x) const
{
  double max_temp = 0.0;
  x.MaxValue(&max_temp);

  if( max_temp > m_max_temp )
    throw std::logic_error("Brusselator_Temperature::compute_problem_response : incoming solution state has max Temp > previously computed value. Shouldn't happen.");

  (*m_g_vec)[0] = 1.5 - m_max_temp;

  cout << " ----  Max Temp = " << m_max_temp << ", and response = " << (*m_g_vec)[0] << endl;

  return;
}
