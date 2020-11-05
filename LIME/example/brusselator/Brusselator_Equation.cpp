
#include "Brusselator_Equation.hpp"

#include "NOX_Common.H"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Basis.hpp"

#include "LIME_Problem_Manager.hpp"

using LIME::Model_Evaluator;

Brusselator_Equation::Brusselator_Equation(LIME::Problem_Manager & pm, const string & name, int numGlobalNodes, const string & outputDir, bool verbose) :
  Model_Evaluator(pm, name),
  m_output_dir(outputDir),
  verbose_(verbose),
  m_xmin(0.0),
  m_xmax(1.0),
  m_dt(1.0e-2),
  m_time(0.0),
  m_max_time_steps(250),
  m_time_step(0),
  useConvection(false)
{
  // allocate state vector
  m_f_epetra_map = Teuchos::rcp(new Epetra_Map(numGlobalNodes, 0, pm.Comm()));
  m_x_epetra_map = Teuchos::rcp(new Epetra_Map(numGlobalNodes, 0, pm.Comm()));
  m_solution = Teuchos::rcp(new Epetra_Vector(*m_x_epetra_map));

  // Create mesh and solution vectors

  // We first initialize the mesh and then the solution since the latter
  // can depend on the mesh.
  m_mesh_x = Teuchos::rcp( new Epetra_Vector(*m_x_epetra_map) );
  double Length= m_xmax - m_xmin;
  m_dx=Length/((double) numGlobalNodes-1);
  for (int i=0; i < m_x_epetra_map->NumMyElements(); i++)
    (*m_mesh_x)[i]=m_xmin + m_dx*((double) m_x_epetra_map->MinMyGID()+i);

  // Create extra vector needed for this transient problem
  m_old_solution = Teuchos::rcp( new Epetra_Vector(*m_x_epetra_map) );

  // Create parameters and responses vectors
  m_p_map = Teuchos::rcp(new Epetra_LocalMap(numParams, 0, pm.Comm()));
  m_p_vec = Teuchos::rcp(new Epetra_Vector(*m_p_map), false);

  m_param_names = Teuchos::rcp( new Teuchos::Array<string>(numParams) );
  (*m_param_names)[0] = "alpha" ; (*m_p_vec)[0]= 0.6;
  (*m_param_names)[1] = "beta"  ; (*m_p_vec)[1]= 2.0;
  (*m_param_names)[2] = "Dcoeff"; (*m_p_vec)[2]= 0.025;
}

//-----------------------------------------------------------------------------

Brusselator_Equation::~Brusselator_Equation()
{
}

//-----------------------------------------------------------------------------

EpetraExt::ModelEvaluator::InArgs 
Brusselator_Equation::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription(my_name_);

  inArgs.setSupports(IN_ARG_x, true);

  inArgs.set_Np(3);

  return inArgs;
}

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  // Get the incoming solution state
  const Epetra_Vector x(*(inArgs.get_x().get()));

  if( outArgs.get_f().get() ) // Signals a residual fill request, eg computeF from NOX
  {
    Epetra_Vector & f = *(outArgs.get_f().get());

    // compute nonlinear residual
    evaluate(NOX::Epetra::Interface::Required::Residual, x, f);
  }

  if( outArgs.get_W().get() ) // Signals either a computeJacobian or a computePreconditioner
  {
    throw "Brusselator_Equation::evalModel : Jacobian matrix support not available.";
  }

  if( 0 < outArgs.Ng() )
  {
    Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(0); 

    if( outArgs.get_g(0).get() ) // Signals a response "g" eval
      compute_problem_response(x); // could pass in parameters if needed
  }
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Map> 
Brusselator_Equation::get_x_map() const 
{ 
  return m_x_epetra_map; 
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Map> 
Brusselator_Equation::get_f_map() const 
{ 
  return m_f_epetra_map; 
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Map> 
Brusselator_Equation::get_g_map(int) const 
{ 
  return Teuchos::null; 
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Vector> 
Brusselator_Equation::get_x_init() const 
{
  return m_solution; 
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Epetra_Vector> 
Brusselator_Equation::get_p_init(int l) const 
{
  if( 0 != l )
    throw "Brusselator_Equation::get_p_init : Only single parameter vector (index = 0) exists.";

  return m_p_vec; 
}

//-----------------------------------------------------------------------------

Teuchos::RCP<const Teuchos::Array<std::string> > 
Brusselator_Equation::get_p_names(int l) const
{
  if( 0 != l )
    throw "Brusselator_Equation::get_p_names : Only single parameter vector (index = 0) exists.";

  return m_param_names; 
}

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::initialize()
{
  initializeSolution();
  if(verbose_)
  {
    string filename(my_name_ + "_init");
    string full_filename = m_output_dir + filename;
    LIME::dump_vector_to_file(*m_solution, full_filename);
    filename += ".txt";
    //write_mesh_and_vector_to_text(*m_solution, filename);
  }
}

//-----------------------------------------------------------------------------
void 
Brusselator_Equation::reinitialize()
{
  m_time = 0.0;
  initializeSolution();
}

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::finalize_solution(const Epetra_Vector& x)
{
  *m_old_solution = x;

  if(verbose_)
  {
    std::ostringstream suffix;
    suffix << m_time_step << std::flush;
    string suffixed_name = createInArgs().modelEvalDescription() + "_" + suffix.str();
    string full_suffixed_name = m_output_dir + suffixed_name;

    LIME::dump_vector_to_file(x, full_suffixed_name);

    suffixed_name += ".txt";

    write_mesh_and_vector_to_text(x, suffixed_name);
  }
}

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::write_mesh_and_vector_to_text(const Epetra_Vector& x, const string & filename)
{
  string full_filename = m_output_dir + filename;
  int numNodes = x.Map().NumMyElements();
  FILE *ifp;
  ifp = fopen(full_filename.c_str(), "w");
  for (int i = 0; i < numNodes; i++)
    fprintf(ifp, "%d  %E  %E \n", i, (*m_mesh_x)[i], x[i]);
  fclose(ifp);
}

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::set_xfer_data(const Teuchos::RCP<Epetra_Vector> & x)
{
  if( Teuchos::is_null(m_aux_data) )
    m_aux_data = Teuchos::rcp(new Epetra_Vector(*x));
  else
    *m_aux_data = *x;
} 

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::initializeSolution()
{
  m_solution->PutScalar(0.0);
  *m_old_solution = *m_solution;
} 

//-----------------------------------------------------------------------------

void 
Brusselator_Equation::generateGraph()
{
  //// Declare required variables
  //int i;
  //int row, column;
  //int numOverlapNodes = OverlapMap->NumMyElements();
  //int OverlapMinMyNodeGID;
  //if (MyPID==0) OverlapMinMyNodeGID = StandardMap->MinMyGID();
  //else OverlapMinMyNodeGID = StandardMap->MinMyGID()-1;
  //
  //// Loop Over # of Finite Elements on Processor
  //for (int ne=0; ne < numOverlapNodes-1; ne++) {
  //        
  //  // Loop over Nodes in Element
  //  for (i=0; i<2; i++) {

  //    // If this node is owned by current processor, add indices
  //    if (StandardMap->MyGID(OverlapMap->GID(ne+i))) {

  //      // Loop over unknowns in Node
  //      row=OverlapMap->GID(ne+i);

  //      // Loop over supporting nodes
  //      for( int j = 0; j < 2; ++j) 
  //      {
  //        // Loop over unknowns at supporting nodes
  //        column=OverlapMap->GID(ne+j);
  //        AA->InsertGlobalIndices(row, 1, &column);
  //      }
  //    }
  //  }
  //}
  //AA->FillComplete();
  //
  //return;
}
