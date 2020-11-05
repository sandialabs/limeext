
#include "FiveLinearEQs.hpp"

using namespace LIME;

double 
Five_Linear_Equations::solve_x1(double x2, double x3, double, bool verbose)
{
  if(verbose)
    cout << "In solve_x1 with x2 = " << x2 << " and x3 = " << x3 << endl;
  double x1 = (7.0 + x2 - x3)/2.0;
  if(verbose)
    cout << "\treturning x1 = " << x1 << endl;
  return x1;
}

double 
Five_Linear_Equations::solve_x2(double x1, double x3, double, bool verbose)
  {
    if(verbose)
      cout << "In solve_x2 with x1 = " << x1 << " and x3 = " << x3 << endl;
    double x2 = (-6.0 -x1 + x3)/2.0;
    if(verbose)
      cout << "\treturning x2 = " << x2 << endl;
    return x2;
  }

double 
Five_Linear_Equations::solve_x3(double x1, double x2, double x5, bool verbose)
  {
    if(verbose)
      cout << "In solve_x3 with x1 = " << x1 << ", x2 = " << x2 << " and x5 = " << x5 << endl;
    double x3 = (6.0 - x1 - x2 + x5)/2.0;
    if(verbose)
      cout << "\treturning x3 = " << x3 << endl;
    return x3;
  }

double 
Five_Linear_Equations::solve_x4(double x1, double x3, double, bool verbose)
  {
    if(verbose)
      cout << "In solve_x4 with x1 = " << x1 << " and x3 = " << x3 << endl;
    double x4 = -x1 + x3;
    if(verbose)
      cout << "\treturning x4 = " << x4 << endl;
    return x4;
  }

double 
Five_Linear_Equations::solve_x5(double x2, double x3, double x4, bool verbose)
  {
    if(verbose)
      cout << "In solve_x5 with x1 = " << x2 << ", x3 = " << x3 << " and x4 = " << x4 << endl;
    double x5 = -10.0 -2.0*x2 + x3 + x4;
    if(verbose)
      cout << "\treturning x5 = " << x5 << endl;
    return x5;
  }

Teuchos::RCP<Five_Linear_Equations::EQ_ModelEval>
Five_Linear_Equations::EQ_ModelEval::create_instance(LIME::Problem_Manager & pm, EQUATION_ID eq_id, bool verbose)
{
  Teuchos::RCP<Five_Linear_Equations::EQ_ModelEval> eq = 
    Teuchos::rcp(new EQ_ModelEval(pm, get_fn(eq_id), eq_id, verbose));

  return eq;
}

double (*Five_Linear_Equations::EQ_ModelEval::get_fn(EQUATION_ID eq_id))(double, double, double, bool)
{
  switch ( eq_id )
  {
    case EQ_ONE:
      return &solve_x1;
      break;
    case EQ_TWO:
      return &solve_x2;
      break;
    case EQ_THREE:
      return &solve_x3;
      break;
    case EQ_FOUR:
      return &solve_x4;
      break;
    case EQ_FIVE:
      return &solve_x5;
      break;
    default:
      throw std::runtime_error("Invalid equation id in EQ_ModelEval get_fn.");
  }
}

std::vector<Teuchos::RCP<LIME::Data_Transfer_Operator> >
Five_Linear_Equations::get_all_transfers(const std::vector<Teuchos::RCP<EQ_ModelEval> > & MEs, bool verbose)
{
  if( 5 != MEs.size() )
    throw std::runtime_error("Incorrect number of EQ_ModelEvals passed to get_all_transfers(...).");
  for( int i = 0; i < 5; ++i )
    if( (int) MEs[i]->get_my_eq_id() != i )
      throw std::runtime_error("Incorrect ordering or content EQ_ModelEvals passed to get_all_transfers(...).");

  std::vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > all_problem_xfers;

  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[1], MEs[0], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[0], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[0], MEs[1], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[1], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[0], MEs[2], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[1], MEs[2], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[4], MEs[2], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[0], MEs[3], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[3], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[1], MEs[4], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[2], MEs[4], verbose) ) );
  all_problem_xfers.push_back( Teuchos::rcp( new EQ_Transfer(MEs[3], MEs[4], verbose) ) );

  return all_problem_xfers;
}

void 
Five_Linear_Equations::EQ_ModelEval::initializeSolution() 
{
  my_x_ = 0.0;
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

    case EQ_FOUR:
      my_known_soln_ = 2.0;
      other_x_[EQ_ONE] = -1.e20;
      other_x_[EQ_THREE] = -1.e20;
      break;

    case EQ_FIVE:
      my_known_soln_ = -1.0;
      other_x_[EQ_ONE] = -1.e20;
      other_x_[EQ_THREE] = -1.e20;
      break;

    default:
      throw std::runtime_error("Invalid equation id in EQ_ModelEval constructor.");
  }
  if(verbose_)
    cout << "Initialized eq " << my_eq_id_ << endl;
}

void 
Five_Linear_Equations::EQ_ModelEval::reinitialize() 
{
  iterate_vals_.clear();
  initializeSolution();
  if(verbose_)
    cout << "Reinitialized eq " << my_eq_id_ << endl;
}

void 
Five_Linear_Equations::EQ_ModelEval::solve_standalone()
{
  double dummy_double = 0.0;

  switch (my_eq_id_)
  {
    case EQ_ONE:
      my_x_ = (*my_func_)(other_x_[EQ_TWO], other_x_[EQ_THREE], dummy_double, verbose_);
      break;
    case EQ_TWO:
      my_x_ = (*my_func_)(other_x_[EQ_ONE], other_x_[EQ_THREE], dummy_double, verbose_);
      break;
    case EQ_THREE:
      my_x_ = (*my_func_)(other_x_[EQ_ONE], other_x_[EQ_TWO], other_x_[EQ_FIVE], verbose_);
      break;
    case EQ_FOUR:
      my_x_ = (*my_func_)(other_x_[EQ_ONE], other_x_[EQ_THREE], dummy_double, verbose_);
      break;
    case EQ_FIVE:
      my_x_ = (*my_func_)(other_x_[EQ_TWO], other_x_[EQ_THREE], other_x_[EQ_FOUR], verbose_);
      break;
    default:
      throw std::runtime_error("Invalid equation id in EQ_ModelEval::solve_standalone.");
  }
  if(verbose_)
    cout << "solve_standalone eq " << my_eq_id_ << endl;
}

bool 
Five_Linear_Equations::EQ_ModelEval::is_converged()
{
  // Store our solution states for later retrieval and comparison
  iterate_vals_.push_back(my_x_);
  bool is_conv =  fabs(my_known_soln_ - my_x_) < FP_UNIT_TEST_TOL;
  if(verbose_)
    cout << "Problem " << my_eq_id_ << " (" << my_id_ << ") is_converged: " << is_conv << endl;
  return is_conv;
}


bool 
Five_Linear_Equations::EQ_Transfer::perform_data_transfer() const
{
  const Teuchos::RCP<EQ_ModelEval> fromME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(source());
  const Teuchos::RCP<EQ_ModelEval>   toME = Teuchos::rcp_dynamic_cast<EQ_ModelEval>(target());
  toME->set_other_x( fromME->get_my_eq_id(), fromME->get_my_x() );
  if(verbose_)
    cout << "Transferred value " << fromME->get_my_x() << " from eq " << fromME->get_my_eq_id() << " (" << fromME->get_id() << ") to eq " << toME->get_my_eq_id() << " (" << toME->get_id() << ")" << endl;

  return true;
}
