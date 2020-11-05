
#include "LIME_Problem_Manager.hpp"
#include "ConvDiff_EQ.hpp"

// A simple 1D linear basis class used for quadrature
namespace {

  class Basis {

    public:

      Basis() :
        eta(0.0),
        wt(0.0)
      { 
        phi = new double[2];
        dphide = new double[2];
      }

      // Destructor
      ~Basis() {
        delete [] phi;
        delete [] dphide;
      }

      // Calculate linear 1D basis at the specified gauss point
      void getBasis(int gp, double *x, double *u)
      {
        int N = 2;
        if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
        if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

        // Calculate basis function and derivatives at nodel pts
        phi[0]=(1.0-eta)/2.0;
        phi[1]=(1.0+eta)/2.0;
        dphide[0]=-0.5;
        dphide[1]=0.5;

        // Caculate basis function and derivative at GP.
        dx=0.5*(x[1]-x[0]);
        xx=0.0;
        uu=0.0;
        duu=0.0;
        for( int i = 0; i < N; ++i ) 
        {
          xx += x[i] * phi[i];
          uu += u[i] * phi[i];
          duu += u[i] * dphide[i];
        }
      }

    public:
      // Variables that are calculated at the gauss point
      double *phi, *dphide; 
      double xx, eta, wt;
      double uu, duu;
      double dx;
  }; // end Basis
}

// Constructor - creates the Epetra objects (maps and vectors) 
ConvDiff_EQ::ConvDiff_EQ(
            LIME::Problem_Manager & pm  ,
            double peclet               ,
            double radiation            ,
            double kappa                ,
            double bcWeight             ,
            INTERFACE_LOCATION iLoc     ,
            double xmin                 ,
            double xmax                 ,
            double Tleft                ,
            double Tright               ,
            int numDofs                 ,
            string name                   ) :
  Model_Evaluator(pm, name)          ,
  xmin_            ( xmin      )   ,
  xmax_            ( xmax      )   ,
  Tleft_           ( Tleft     )   ,
  Tright_          ( Tright    )   ,
  peclet_          ( peclet    )   ,
  radiation_       ( radiation )   ,
  kappa_           ( kappa     )   ,
  iLocation_       ( iLoc      )   ,
  bcWeight_        ( bcWeight  )   ,
  interface_T_from_xfer_(0.0)      ,
  interface_Flux_from_xfer_(0.0)
{
  // Create mesh and solution vectors
  epetra_map_ = Teuchos::rcp(new Epetra_Map(numDofs, 0, pm.Comm()));

  // We first initialize the mesh and then the solution since the latter can depend on the mesh.

  mesh_ = Teuchos::rcp( new Epetra_Vector(*epetra_map_) );
  double dx   = (xmax_ - xmin_) / ( (double) numDofs - 1 );

  for( int i = 0; i < numDofs; ++i ) 
    (*mesh_)[i]=xmin_ + dx*((double) epetra_map_->MinMyGID()+i);

  // Create and initialize (using default provided) the solution vector
  me_interface_soln_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));

  initialize();
}

//-----------------------------------------------------------------------------

EpetraExt::ModelEvaluator::InArgs 
ConvDiff_EQ::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  // Register our identify
  inArgs.setModelEvalDescription(my_name_);

  // Signal that we can do calculations with incoming state x
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}

//-----------------------------------------------------------------------------

EpetraExt::ModelEvaluator::OutArgs 
ConvDiff_EQ::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  // Register our identify - consistent with createInArgs
  outArgs.setModelEvalDescription(my_name_);

  // Signal that we can compute a residual vector
  outArgs.setSupports(OUT_ARG_f, true);

  return outArgs;
}

//-----------------------------------------------------------------------------

// Initialize method
void 
ConvDiff_EQ::initialize()
{

  // Set exact interface temperature accordingly
  if( LEFT == iLocation_ )
    T1_exact_ = Tleft_ ;
  else
    T1_exact_ = Tright_;

  // Choose appropriate element and nodes for interfacial BC enforcement
  int OverlapNumMyNodes = epetra_map_->NumMyElements(); // this may break in parallel
  if( LEFT == iLocation_ )
  {
    interface_elem_ = 0                     ;
    local_node_     = 0                     ;
    interface_node_ = 0                     ;
    opposite_node_  = OverlapNumMyNodes - 1 ;
    dirScale_       = 1.0                   ;
  }
  else
  {
    interface_elem_ = OverlapNumMyNodes - 2 ;
    local_node_     = 1                     ;
    interface_node_ = OverlapNumMyNodes - 1 ;
    opposite_node_  = 0                     ;
    dirScale_       = -1.0                  ;
  }

  // Use this method to compute and output the analytic solution and its first derivative
  exactSolution_     = Teuchos::rcp( new Epetra_Vector(*me_interface_soln_) );
  dTdx_              = Teuchos::rcp( new Epetra_Vector(*me_interface_soln_) );
  Epetra_Vector & x = *mesh_;

  if( 1.e-20 < fabs(peclet_) )
  {
    for( int i = 0; i < x.MyLength(); ++i ) 
    {
      (*exactSolution_)[i] = (T1_exact_ - Tleft_*exp(peclet_) + (Tleft_ - T1_exact_)*exp(peclet_*(x[i]-xmin_))) /
                            ( 1.0 - exp(peclet_) );
      (*dTdx_         )[i] = peclet_ * (Tleft_ - T1_exact_)*exp(peclet_*(x[i]-xmin_)) /
                            ( 1.0 - exp(peclet_) );
    }
  }
  else
  {
    for( int i = 0; i < x.MyLength(); ++i ) 
    {
      (*exactSolution_)[i] = (Tright_ - T1_exact_)*(x[i] - xmax_) + Tright_  ;
      (*dTdx_         )[i] =  Tright_ - T1_exact_                          ;
    }
  }

  ostringstream sval; 
  std::string fileName1 = my_name_+"_analytic_" + sval.str(); 
  std::string fileName2 = my_name_+"_dTdx_"     + sval.str(); 

  ofstream outFile1(fileName1.c_str());
  ofstream outFile2(fileName2.c_str());
  if( !outFile1 || !outFile2 )
  {
    std::string msg = "ERROR: Could not open one of files \"" + fileName1 + "\"" + " or \""
                      + fileName2 + "\"";
    throw msg;
  }

  for( int i = 0; i < x.MyLength(); ++i ) 
  {
    outFile1 << i << "  " << x[i] << "  " << (*exactSolution_)[i] << endl;
    outFile2 << i << "  " << x[i] << "  " << (*dTdx_         )[i] << endl;
  }
}

//-----------------------------------------------------------------------------

// Set me_interface_soln_ to desired initial condition
void 
ConvDiff_EQ::initializeSolution()
{
  me_interface_soln_->PutScalar(0.5);

  //Epetra_Vector perturb(soln);
  //perturb.Random();
  //perturb.Scale(0.20);

  //soln.Update(1.0, perturb, 1.0);
} 

//-----------------------------------------------------------------------------

void 
ConvDiff_EQ::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  // Create a "View" of the incoming solution state, x
  //  This avoids a copy of data and becomes important for real-sized problems
  const Epetra_Vector x(View, *(inArgs.get_x().get()), 0);

  //cout << "CON1D_ModelEval_w_Resid::evalModel called. state:\n";

  if( outArgs.get_f().get() ) // A non-NULL f-vector signals a residual fill request
  {
    // Get a reference to the vector we will populate
    Epetra_Vector & f = *(outArgs.get_f().get());

    // compute our application's nonlinear residual
    evaluate(NOX::Epetra::Interface::Required::Residual, x, f);
  }

  // We don't yet support this and respond to such a request with an error
  if( outArgs.get_W().get() ) // Signals either a computeJacobian or a computePreconditioner
    throw std::runtime_error("ConvDiff_EQ::evalModel : Jacobian matrix support not available yet.");
}

//-----------------------------------------------------------------------------

// Matrix and Residual Fills
bool 
ConvDiff_EQ::evaluate( NOX::Epetra::Interface::Required::FillType flag,
		    const Epetra_Vector & soln, 
		    Epetra_Vector & rhs) const
{

  // Zero out the objects that will be filled
  rhs.PutScalar(0.0);

  int row;
  double * xx    = new double[2];
  double * uu    = new double[2]; 

  Basis basis;

  // Loop Over # of Finite Elements on Processor
  for( int ne = 0; ne < epetra_map_->NumMyElements()-1; ++ne )
  {
    // Loop Over Gauss Points
    for( int gp = 0; gp < 2; ++gp ) 
    {
      // Get the solution and coordinates at the nodes 
      xx[0]=(*mesh_)[ne];
      xx[1]=(*mesh_)[ne+1];
      uu[0] = soln[ne];
      uu[1] = soln[ne+1];
      // Calculate the basis function and variables at the gauss points
      basis.getBasis(gp, xx, uu);

      // Loop over Nodes in Element
      for( int i = 0; i < 2; ++i )
      {
	row = epetra_map_->GID(ne+i);
	if( epetra_map_->MyGID(row) ) 
        {
          rhs[epetra_map_->LID(epetra_map_->GID(ne+i))] +=
            + basis.wt * basis.dx
            * ( peclet_ * (basis.duu / basis.dx) * basis.phi[i] 
            +   kappa_ * (1.0/(basis.dx*basis.dx)) * basis.duu * basis.dphide[i] );
	}
      }
    }
  } 

  //if( NOX::Epetra::Interface::Required::Residual == flag )
  //{
  //  int lastDof = epetra_map_->LID(epetra_map_->MaxAllGID());
  //  cout << "\t\"" << myName << "\" u[0] = " << (*soln)[0] 
  //       << "\tu[N] = " << (*soln)[lastDof] << endl;
  //  cout << "\t\"" << myName << "\" RHS[0] = " << (*rhs)[0] 
  //       << "\tRHS[N] = " << (*rhs)[lastDof] << endl << endl;
  //}


  // Apply BCs

  double heatFlux = compute_heat_flux(soln);

  double bcResidual = bcWeight_        * (heatFlux - interface_Flux_from_xfer_ ) -
                      (1.0 - bcWeight_)* (soln[interface_node_] - interface_T_from_xfer_ );

  int lastDof = epetra_map_->LID(epetra_map_->MaxAllGID());

  // "Left" boundary
  if( LEFT == iLocation_ ) // this may break in parallel
  {
    rhs[0]       = bcResidual;
    rhs[lastDof] = soln[lastDof] - Tright_;
  }
  // "Right" boundary
  else
  {
    rhs[0]       = soln[0] - Tleft_;
    rhs[lastDof] = bcResidual;
  }

  //if( flag == NOX::Epetra::Interface::Required::Residual )
  //{
  //  cout << setprecision(10) << "ConvDiff_EQ::evaluate : \"" << my_name_ << "\" x_Vector :" << endl << soln << endl;
  //  cout << setprecision(10) << "ConvDiff_EQ::evaluate : \"" << my_name_ << "\" rhs_Vector :" << endl << rhs << endl;
  //}

  delete [] xx    ;
  delete [] uu    ;

  return true;
}

//-----------------------------------------------------------------------------

// A fill specialized to the single node at the coupling interface
double 
ConvDiff_EQ::compute_heat_flux(const Epetra_Vector & soln) const
{
  // Declare required variables
  int row;
  double * xx = new double[2];
  double * uu = new double[2]; 

  Basis basis;

  double flux = 0.0;

  // Loop Over Gauss Points
  for( int gp = 0; gp < 2; ++gp ) 
  {
    // Get the solution and coordinates at the nodes 
    xx[0]=(*mesh_)[interface_elem_];
    xx[1]=(*mesh_)[interface_elem_+1];
    uu[0] = soln[interface_elem_];
    uu[1] = soln[interface_elem_+1];

    // Calculate the basis function and variables at the gauss points
    basis.getBasis(gp, xx, uu);

    row = epetra_map_->GID( interface_elem_ + local_node_ );

    if( epetra_map_->MyGID(row) ) 
    {
      flux += 
        + basis.wt * basis.dx
        * ( peclet_ * (basis.duu / basis.dx) * basis.phi[local_node_] 
        +   kappa_ * (1.0/(basis.dx*basis.dx)) * basis.duu * basis.dphide[local_node_] );
    }
  }

  delete [] xx    ;
  delete [] uu    ;

  //int lastDof = epetra_map_->LID(epetra_map_->MaxAllGID());
  //cout << "\t\"" << myName << "\" u[0] = " << u[0] 
  //     << "\tu[N] = " << u[lastDof] << endl;
  //cout << u << endl;
  //cout << "\t\"" << myName << "\" flux = " << flux << endl << endl;

  // Scale domain integration according to interface position
  flux *= dirScale_;

  // Now add radiation contribution to flux
  flux += radiation_ * ( pow(soln[interface_node_], 4) - pow(soln[opposite_node_], 4) );

  return flux;
}

//-----------------------------------------------------------------------------

double 
ConvDiff_EQ::get_interface_temp() const
{
  
  if( LEFT == iLocation_ )
    return (*me_interface_soln_)[0];
  else
  {
    int lastDof = epetra_map_->LID(epetra_map_->MaxAllGID());
    return (*me_interface_soln_)[lastDof];
  }
}

//-----------------------------------------------------------------------------

void 
ConvDiff_EQ::outputStatus( ostream & os ) 
{
  
  //std::string location = ( iLocation_ == LEFT ) ? "Left" : "Right";

  //os << "\"" << myName << "\" couples at an interface with Problem \"" << depProbPtr->getName()
  //   << "\" on the " << location << endl;

  return;
}

//-----------------------------------------------------------------------------

double 
ConvDiff_EQ::computeAnalyticInterfaceTemp(
                  double radiation      ,
                  double T_left         ,
                  double T_right        ,
                  double kappa          ,
                  double peclet           )
{
  
  int    nIters         = 0;            // iteration counter
  double T_int          = 1.0;          // initial guess
  double tol            = 1.e-12;       // solve tolerance
  double residual       = radiation*( pow(T_int, 4) - pow(T_right, 4) ) 
                          + (T_left - T_int)*peclet*exp(peclet)/(1.0 - exp(peclet))
                          - kappa*(T_right - T_int);
  double dfdT           = 0.0;

  // Simple Newton loop
  while( tol < fabs(residual) )
  {
    ++nIters;

    dfdT = radiation*4.0*pow(T_int,3)
           - peclet*exp(peclet)/(1.0 - exp(peclet))
           + kappa;

    if( 1.e-15 > fabs(dfdT) )
    {
      cout << "ConvDiff_EQ::computeAnalyticInterfaceTemp:\n"
              "Warning: Obtained near-zero derivtative. Aborting calculation." << endl;
      return(T_int);
    }

    T_int = T_int - residual/dfdT;

    residual = radiation*( pow(T_int, 4) - pow(T_right, 4) ) 
               + (T_left - T_int)*peclet*exp(peclet)/(1.0 - exp(peclet))
               - kappa*(T_right - T_int);
  }

  cout << "Analytic interfacial temperature = " << setprecision(8) << T_int << endl;
  cout << "Residual = " << residual << " in " << nIters << " iterations." << endl;

  return T_int;
}

//-----------------------------------------------------------------------------
