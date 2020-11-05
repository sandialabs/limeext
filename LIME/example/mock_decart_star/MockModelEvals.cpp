
#include "MockModelEvals.hpp"

using namespace LIME;

// -----------------------------------------------------------------------------------

Mock_DeCART_APP::Mock_DeCART_APP(double H_rxn, double alpha, double beta, double a, double b) :
  H_rxn_(H_rxn),
  alpha_(alpha),
  beta_(beta),
  a_(a),
  b_(b)
{
}

double
Mock_DeCART_APP::solve_temp(double cdot, double rho, bool verbose) const
{
  double qdot = H_rxn_ * cdot;

  return (-beta_*rho + sqrt(beta_*beta_*rho*rho+ 4.0*alpha_*qdot))/(2.0*alpha_);
}

double
Mock_DeCART_APP::compute_depletion_rate(double krxn, double cval, bool verbose) const
{
  return -krxn*cval;
}

double
Mock_DeCART_APP::compute_rxn_coeff(double temp, bool verbose) const
{
  if( verbose )
    cout << "\t\t\tMock_DeCART_APP::compute_rxn_coeff using temp = " << temp;

  double krxn = a_ * (1.0 + b_*temp);

  if( verbose )
    cout << ", krxn = " << krxn << endl;

  return krxn;
}

double
Mock_DeCART_APP::compute_energy_error(double cdot, double rho, double temp, bool verbose) const
{
  double lhs = H_rxn_ * cdot;
  double rhs = alpha_*temp*temp + beta_*rho*temp;

  return (lhs - rhs);
}

double
Mock_DeCART_APP::compute_depletion_error(double cdot, double krxn, double cval, bool verbose) const
{
  double lhs = cdot;
  double rhs = -krxn*cval;

  return (lhs - rhs);
}

// -----------------------------------------------------------------------------------

Mock_DeCART_Energy_ME::Mock_DeCART_Energy_ME(LIME::Problem_Manager & pm, Mock_DeCART_APP & app, double temp0, bool verbose) : 
  Model_Evaluator(pm),
  decart_app_(app),
  verbose_(verbose),
  temp_(temp0)
{
}

void 
Mock_DeCART_Energy_ME::solve_standalone()
{
  if(verbose_)
    cout << "In Mock_DeCART_Energy_ME::solve_standalone, cdot = " << cdot_xfer_ << ", rho = " << rho_xfer_ << endl;

  temp_ = decart_app_.solve_temp(cdot_xfer_, rho_xfer_, verbose_);

  if(verbose_)
    cout << "Found temp = " << temp_ << endl;
}

bool 
Mock_DeCART_Energy_ME::is_converged()
{
  double error = decart_app_.compute_energy_error(cdot_xfer_, rho_xfer_, temp_, verbose_);
  bool converged = (abs(error) < 1.e-6);

  if(verbose_)
    cout << "\nMock_DeCART_Energy_ME error: " << error << endl;

  return converged;
}

// -----------------------------------------------------------------------------------

Mock_DeCART_Depletion_ME::Mock_DeCART_Depletion_ME(LIME::Problem_Manager & pm, Mock_DeCART_APP & app, double cval0, bool verbose) : 
  Model_Evaluator(pm),
  decart_app_(app),
  verbose_(verbose),
  cval_(cval0),
  cval_old_(cval0),
  time_(0.0),
  dt_(0.1),
  time_step_(0)
{
}

void 
Mock_DeCART_Depletion_ME::update_time()
{
  //// We could do something more elaborate here but shoose for now to do simple EEuler step
  //// and reset the old cval to the current
  //cdot_ = decart_app_.compute_depletion_rate(krxn_xfer_, cval_old_, verbose_);
  cval_ = cval_old_ + dt_*cdot_;
  cval_old_ = cval_;
  time_ += dt_;
  time_step_++;

  if(verbose_)
    cout << "In Mock_DeCART_Depletion_ME::update_time, advancing cval to next time: cval = " << cval_old_ << endl;
}

void 
Mock_DeCART_Depletion_ME::solve_standalone()
{
  if(verbose_)
    cout << "In Mock_DeCART_Depletion_ME::solve_standalone, krxn = " << krxn_xfer_ << endl;

  // Simple Explicit Euler method to get solution at next time step
  cdot_ = decart_app_.compute_depletion_rate(krxn_xfer_, cval_old_, verbose_);
  //cval_ = cval_old_ + dt_*cdot_;


  if(verbose_)
    cout << "Found cval = " << cval_ << ", (cdot = " << cdot_ << ")" << endl;
}

bool 
Mock_DeCART_Depletion_ME::is_converged()
{
  double error = decart_app_.compute_depletion_error(cdot_, krxn_xfer_, cval_, verbose_);
  bool converged = (abs(error) < 1.e-6);

  if(verbose_)
    cout << "\nMock_DeCART_Depletion_ME error: " << error << endl;

  return converged;
}

// -----------------------------------------------------------------------------------

Mock_Star_APP::Mock_Star_APP(double rho0, double gamma) :
  rho0_(rho0),
  gamma_(gamma)
{
}

double
Mock_Star_APP::solve_density(double temp, bool verbose)
{
  return rho0_ * (1.0 + gamma_*temp);
}

double
Mock_Star_APP::compute_density_error(double temp, double rho, bool verbose) const
{
  double lhs = rho;
  double rhs = rho0_ * (1.0 + gamma_*temp);

  return (lhs - rhs);
}

// -----------------------------------------------------------------------------------

Mock_Star_ME::Mock_Star_ME(LIME::Problem_Manager & pm, Mock_Star_APP & app, bool verbose) :
  Model_Evaluator(pm),
  star_app_(app),
  verbose_(verbose)
{
  rho_ = star_app_.get_rho0();
}

void 
Mock_Star_ME::solve_standalone()
{
  if(verbose_)
    cout << "In Mock_Star_ME::solve_standalone using temp = " << temp_xfer_ << endl;

  rho_ = star_app_.solve_density(temp_xfer_, verbose_);

  if(verbose_)
    cout << "Found rho = " << rho_ << endl;
}

bool 
Mock_Star_ME::is_converged()
{
  double error = star_app_.compute_density_error(temp_xfer_, rho_, verbose_);
  bool converged = (abs(error) < 1.e-6);

  if(verbose_)
    cout << "\nMock_Star_ME error: " << error << endl;

  return converged;
}

// -----------------------------------------------------------------------------------

bool
Transfer_Mock_DeCART_To_Star::perform_data_transfer() const
{ 
  const Teuchos::RCP<Mock_DeCART_Energy_ME> & decart_energy_me = Teuchos::rcp_dynamic_cast<Mock_DeCART_Energy_ME>(source());
  const Teuchos::RCP<Mock_Star_ME>          & star_me          = Teuchos::rcp_dynamic_cast<Mock_Star_ME>(target());

  double temp = decart_energy_me->get_temp();

  if(verbose_)
    cout << "\tIn Transfer_Mock_DeCART_To_Star::perform_data_transfer, sending temp = " << temp << " to Mock_Star." << endl;

  star_me->set_temp_xfer(temp);
   
  return true;
}

// -----------------------------------------------------------------------------------

bool
Transfer_Mock_Star_To_DeCART::perform_data_transfer() const
{ 
  const Teuchos::RCP<Mock_Star_ME>          & star_me          = Teuchos::rcp_dynamic_cast<Mock_Star_ME>(source());
  const Teuchos::RCP<Mock_DeCART_Energy_ME> & decart_energy_me = Teuchos::rcp_dynamic_cast<Mock_DeCART_Energy_ME>(target());

  double rho = star_me->get_rho();

  if(verbose_)
    cout << "\tIn Transfer_Mock_Star_To_DeCART::perform_data_transfer, sending rho = " << rho << " to Mock_DeCART_Energy_ME." << endl;
   
  decart_energy_me->set_rho_xfer(rho);

  return true;
}

// -----------------------------------------------------------------------------------

bool
Transfer_Mock_DeCART_Energy_To_Depletion::perform_data_transfer() const
{ 
  const Teuchos::RCP<Mock_DeCART_Energy_ME>    & decart_energy_me    = Teuchos::rcp_dynamic_cast<Mock_DeCART_Energy_ME>(source());
  const Teuchos::RCP<Mock_DeCART_Depletion_ME> & decart_depletion_me = Teuchos::rcp_dynamic_cast<Mock_DeCART_Depletion_ME>(target());

  // Compute the rate coefficient using temperature from Mock_DeCART_Energy_ME
  double krxn = decart_energy_me->get_decart_app().compute_rxn_coeff(decart_energy_me->get_temp(), verbose_);

  if(verbose_)
    cout << "\tIn Transfer_Mock_DeCART_Energy_To_Depletion::perform_data_transfer, sending krxn = " << krxn << " to Mock_DeCART_Depletion_ME." << endl;
   
  decart_depletion_me->set_krxn_xfer(krxn);

  return true;
}

// -----------------------------------------------------------------------------------

bool
Transfer_Mock_DeCART_Depletion_To_Energy::perform_data_transfer() const
{ 
  const Teuchos::RCP<Mock_DeCART_Depletion_ME> & decart_depletion_me = Teuchos::rcp_dynamic_cast<Mock_DeCART_Depletion_ME>(source());
  const Teuchos::RCP<Mock_DeCART_Energy_ME>    & decart_energy_me    = Teuchos::rcp_dynamic_cast<Mock_DeCART_Energy_ME>(target());

  double cdot = decart_depletion_me->get_cdot();

  if(verbose_)
    cout << "\tIn Transfer_Mock_DeCART_Depletion_To_Energy::perform_data_transfer, sending cdot = " << cdot << " to Mock_DeCART_Energy_ME." << endl;
   
  decart_energy_me->set_cdot_xfer(cdot);

  return true;
}

