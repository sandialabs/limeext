
#include "ThreeNonlinEQsRes.hpp"

using namespace LIME;

// -----------------------------------------------------------------------------------

double 
Three_Nonlinear_Equations_w_Resid::residual_1(double x1, double x2, double x3, bool verbose)
{
  if(verbose)
    cout << "In residual_1 x1 = " << x1 << ", x2 = " << x2 << " and x3 = " << x3 << endl;
  double resid = 2.0*x1 - x2 + pow(x3,3) - 31.0;
  if(verbose)
    cout << "\treturning resid1 = " << resid << endl;
  return resid;
}

// -----------------------------------------------------------------------------------

double 
Three_Nonlinear_Equations_w_Resid::residual_2(double x1, double x2, double x3, bool verbose)
{
  if(verbose)
    cout << "In residual_2 with x2 = " << x2 << ", x1 = " << x1 << " and x3 = " << x3 << endl;
  double resid = x1 + 2.0*x2 - 2.0*x3*x3 + 21.0;
  if(verbose)
    cout << "\treturning resid2 = " << resid << endl;
  return resid;
}

// -----------------------------------------------------------------------------------

double 
Three_Nonlinear_Equations_w_Resid::residual_3(double x1, double x2, double x3, bool verbose)
{
  if(verbose)
    cout << "In residual_3 with x3 = " << x3 << ", x1 = " << x1 << ", x2 = " << x2 << endl;
  double resid = x1*x1 + x2 +2.0*x3 - 5.0;
  if(verbose)
    cout << "\treturning resid3 = " << resid << endl;
  return resid;
}

// -----------------------------------------------------------------------------------

Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::EQ_ModelEval(LIME::Problem_Manager & pm, double (*f_ptr)(double, double, double, bool), EQUATION_ID my_id, bool verbose) :
  LIME::Model_Evaluator(pm, "EQ_ModelEval"),
  verbose_(verbose),
  my_eq_id_(my_id),
  my_func_(f_ptr)
{
  epetra_map_ = Teuchos::rcp( new Epetra_Map(1, 0, pm.Comm()) );
  soln_vec_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
  initializeSolution();
}

// -----------------------------------------------------------------------------------

const vector<vector<double> > &
Three_Nonlinear_Equations_w_Resid::get_gold_jacobian_coeffs()
{
  static vector< vector<double> > gold_coeffs;
  if( gold_coeffs.empty() )
  {
    gold_coeffs.resize(3);
    for( size_t i = 0; i < gold_coeffs.size(); ++i )
      gold_coeffs[i].resize(3);
    gold_coeffs[0][0] = 2.0; gold_coeffs[0][1] = -1.0; gold_coeffs[0][2] =  27.0;
    gold_coeffs[1][0] = 1.0; gold_coeffs[1][1] =  2.0; gold_coeffs[1][2] = -12.0;
    gold_coeffs[2][0] = 2.0; gold_coeffs[2][1] =  1.0; gold_coeffs[2][2] =   2.0;
  }

  return gold_coeffs;
}

// -----------------------------------------------------------------------------------


EpetraExt::ModelEvaluator::InArgs 
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  // Register our identify
  inArgs.setModelEvalDescription(my_name_);

  // Signal that we can do calculations with incoming state x
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}

// -----------------------------------------------------------------------------------

EpetraExt::ModelEvaluator::OutArgs 
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  // Register our identify - consistent with createInArgs
  outArgs.setModelEvalDescription(my_name_);

  // Signal that we can compute a residual vector
  outArgs.setSupports(OUT_ARG_f, true);

  return outArgs;
}

// -----------------------------------------------------------------------------------

void 
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const
{
  // Copy incoming solution state into our own solution vector
  (*soln_vec_) = *(inArgs.get_x());
  const Epetra_Vector & x = *soln_vec_;

  if( outArgs.get_f().get() ) // A non-NULL f-vector signals a residual fill request
  {
    // Get a reference to the vector we will populate
    Epetra_Vector & f = *(outArgs.get_f().get());

    switch (my_eq_id_)
    {
      case EQ_ONE:
        f[0] = (*my_func_)(x[0], other_x_[EQ_TWO], other_x_[EQ_THREE], verbose_);
        break;
      case EQ_TWO:
        f[0] = (*my_func_)(other_x_[EQ_ONE], x[0], other_x_[EQ_THREE], verbose_);
        break;
      case EQ_THREE:
        f[0] = (*my_func_)(other_x_[EQ_ONE], other_x_[EQ_TWO], x[0], verbose_);
        break;
      default:
        throw std::runtime_error("Invalid equation id in EQ_ModelEval::solve_standalone.");
    }
  }
  else
    // We don't support any other callbacks and so should never get
    // called for anything other than a residual fill.
    throw std::runtime_error("ResEQ_ModelEval::evalModel : Only residual fills are supported.");
}

// -----------------------------------------------------------------------------------

Teuchos::RCP<Three_Nonlinear_Equations_w_Resid::EQ_ModelEval>
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::create_instance(LIME::Problem_Manager & pm, EQUATION_ID eq_id, bool verbose)
{
  Teuchos::RCP<Three_Nonlinear_Equations_w_Resid::EQ_ModelEval> eq = 
    Teuchos::rcp(new EQ_ModelEval(pm, get_fn(eq_id), eq_id, verbose));

  return eq;
}

// -----------------------------------------------------------------------------------

double (*Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::get_fn(EQUATION_ID eq_id))(double, double, double, bool)
{
  switch ( eq_id )
  {
    case EQ_ONE:
      return &residual_1;
      break;
    case EQ_TWO:
      return &residual_2;
      break;
    case EQ_THREE:
      return &residual_3;
      break;
    default:
      throw std::runtime_error("Invalid equation id in EQ_ModelEval get_fn.");
  }
}

// -----------------------------------------------------------------------------------

std::vector<Teuchos::RCP<LIME::Data_Transfer_Operator> >
Three_Nonlinear_Equations_w_Resid::get_all_transfers(const std::vector<Teuchos::RCP<EQ_ModelEval> > & MEs, bool verbose)
{
  if( 3 != MEs.size() )
    throw std::runtime_error("Incorrect number of EQ_ModelEvals passed to get_all_transfers(...).");
  for( int i = 0; i < 3; ++i )
    if( (int) MEs[i]->get_my_eq_id() != i )
      throw std::runtime_error("Incorrect ordering or content EQ_ModelEvals passed to get_all_transfers(...).");

  std::vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > all_problem_xfers;

  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[1], MEs[0], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[0], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[0], MEs[1], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[1], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[0], MEs[2], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[1], MEs[2], verbose) ) );

  return all_problem_xfers;
}

// -----------------------------------------------------------------------------------

void 
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::initializeSolution() 
{
  (*soln_vec_)[0] = double(my_eq_id_+1)+double(my_eq_id_+1)/10.0;
  switch (my_eq_id_) {

    case EQ_ONE:
      my_known_soln_ = 1.0;
      other_x_[EQ_TWO] = -1.e20;
      other_x_[EQ_THREE] = -1.e20;
      break;

    case EQ_TWO:
      my_known_soln_ = -2.0;
      other_x_[EQ_ONE] = -1.e20;
      other_x_[EQ_THREE] = -1.e20;
      break;

    case EQ_THREE:
      my_known_soln_ = 3.0;
      other_x_[EQ_ONE] = -1.e20;
      other_x_[EQ_TWO] = -1.e20;
      break;

    default:
      throw std::runtime_error("Invalid equation id in EQ_ModelEval constructor.");
  }
  if(verbose_)
    cout << "Initialized eq " << my_eq_id_ << endl;
}

// -----------------------------------------------------------------------------------

void 
Three_Nonlinear_Equations_w_Resid::EQ_ModelEval::reinitialize() 
{
  iterate_vals_.clear();
  initializeSolution();
  if(verbose_)
    cout << "Reinitialized eq " << my_eq_id_ << endl;
}

// -----------------------------------------------------------------------------------

bool 
Three_Nonlinear_Equations_w_Resid::EQ_Transfer::perform_data_transfer() const
{
  const Teuchos::RCP<EQ_ModelEval> fromME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(source());
  const Teuchos::RCP<EQ_ModelEval>   toME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(target());
  toME->set_other_x( fromME->get_my_eq_id(), fromME->get_my_x() );
  if(verbose_)
    cout << "Transferred value " << fromME->get_my_x() << " from eq " << fromME->get_my_eq_id() 
         << " (" << fromME->get_id() << ") to eq " << toME->get_my_eq_id() << " (" << toME->get_id() 
         << ")" << endl;

  return true;
}
