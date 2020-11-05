#ifndef MOCK_DESTAR_MODEL_EVALS_HPP
#define MOCK_DESTAR_MODEL_EVALS_HPP

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

namespace LIME 
{
  class Mock_DeCART_APP
  {
    public:

      Mock_DeCART_APP(double H_rxn, double alpha, double beta, double a, double b);
      ~Mock_DeCART_APP() {}

    // -----------------------------------------------------------------------------------
    // Contrive and describe ...
    //
    //
    // The system has ???
    // -----------------------------------------------------------------------------------

    double solve_temp(double qdot, double rho, bool verbose) const;
    double compute_depletion_rate(double krxn, double cval, bool verbose) const;
    double compute_rxn_coeff(double temp, bool verbose) const;
    double compute_energy_error(double cdot, double rho, double temp, bool verbose) const;
    double compute_depletion_error(double cdot, double krxn, double cval, bool verbose) const;

    private:

      double H_rxn_;
      double alpha_;
      double beta_;
      double a_;
      double b_;
  };

  class Mock_Star_APP
  {
    public:

      Mock_Star_APP(double rho0, double gamma);
      ~Mock_Star_APP() {}

    // -----------------------------------------------------------------------------------
    // Contrive and describe ...
    //
    //
    // The system has ???
    // -----------------------------------------------------------------------------------

    double get_rho0() const
      { return rho0_; }

    double solve_density(double temp, bool verbose);
    double compute_density_error(double temp, double rho, bool verbose) const;

    private:

      double rho0_;
      double gamma_;
  };

  // -------------------------------------------------------------
  // A Simple Model Evaluator for Mock DeCART
  // -------------------------------------------------------------
  class Mock_DeCART_Energy_ME : public LIME::Model_Evaluator
  {
    public:

      Mock_DeCART_Energy_ME(LIME::Problem_Manager & pm, Mock_DeCART_APP &, double temp0, bool verbose);

      ~Mock_DeCART_Energy_ME() {}

      bool supports_standalone_solve() const { return true; }

      bool is_transient() const { return false; }

      void solve_standalone();
      bool is_converged();

      // Data access
      double get_temp() const
        { return temp_; }

      void set_rho_xfer( double val )
        { rho_xfer_ = val; }

      void set_cdot_xfer( double val )
        { cdot_xfer_ = val; }

      const Mock_DeCART_APP & get_decart_app() const
        { return decart_app_; }

    private:

      Mock_DeCART_APP & decart_app_;
      bool verbose_;
      bool converged_;

      // my dof
      double temp_;

      // Data to be populated from Xfer from Mock_Star
      double cdot_xfer_;
      double rho_xfer_;
  };


  class Mock_DeCART_Depletion_ME : public LIME::Model_Evaluator
  {
    public:

      Mock_DeCART_Depletion_ME(LIME::Problem_Manager & pm, Mock_DeCART_APP &, double cval0, bool verbose);

      ~Mock_DeCART_Depletion_ME() {}

      bool supports_standalone_solve() const { return true; }

      bool is_transient() const { return true; }

      unsigned int get_max_steps() const { return 1; }
      double get_time_step() const { return dt_; }
      double get_max_time() const { return 0.5; }
      void update_time();

      void solve_standalone();
      bool is_converged();

      double get_cdot() const
        { return cdot_; }

      void set_krxn_xfer( double val )
        { krxn_xfer_ = val; }

      const Mock_DeCART_APP & get_decart_app() const
        { return decart_app_; }

    private:

      Mock_DeCART_APP & decart_app_;
      bool verbose_;
      bool converged_;

      // my dof
      double cval_;

      // previous state
      double cval_old_;

      // auxiliary data
      double cdot_;

      // time info
      double time_;
      double dt_;
      int time_step_;

      // Data to be populated from Xfer from Mock_Star
      double krxn_xfer_;
  };


  namespace Mock_Star_EQs
  {
    // -----------------------------------------------------------------------------------
    // Contrive and describe ...
    //
    //
    // The system has ???
    // -----------------------------------------------------------------------------------

    double solve_temp_and_density(double power, bool verbose);
  }

  class Mock_Star_ME : public LIME::Model_Evaluator
  {
    public:

      Mock_Star_ME(LIME::Problem_Manager & pm, Mock_Star_APP &, bool verbose = false);

      ~Mock_Star_ME() {}

      bool supports_standalone_solve() const { return true; }

      bool is_transient() const { return false; }

      void solve_standalone();
      bool is_converged();

      // Data access
      double get_rho() const
        { return rho_; }

      void set_temp_xfer( double val )
        { temp_xfer_ = val; }

    private:

      Mock_Star_APP & star_app_;
      bool verbose_;
      bool converged_;

      // my dof
      double rho_;

      // Data to be populated from Xfer from Mock_Star
      double temp_xfer_;
  };

  // -------------------------------------------------------------------------
  // A Simple Transfer operator for moving data between Mock_DeCART and Mock_Star Model Evaluators
  // -------------------------------------------------------------------------
  class Transfer_Mock_DeCART_To_Star : public LIME::Data_Transfer_Operator 
  {
    public:

      Transfer_Mock_DeCART_To_Star(const Teuchos::RCP<Mock_DeCART_Energy_ME> from, Teuchos::RCP<Mock_Star_ME> to, bool verbose = false ) :
        LIME::Data_Transfer_Operator(from, to),
        verbose_(verbose)
    { }

      bool perform_data_transfer() const;

    private:
      bool verbose_;
  };

  // -------------------------------------------------------------------------
  // A Simple Transfer operator for moving data between Mock_DeCART and Mock_Star Model Evaluators
  // -------------------------------------------------------------------------
  class Transfer_Mock_Star_To_DeCART : public LIME::Data_Transfer_Operator 
  {
    public:

      Transfer_Mock_Star_To_DeCART(const Teuchos::RCP<Mock_Star_ME> from, Teuchos::RCP<Mock_DeCART_Energy_ME> to, bool verbose = false ) :
        LIME::Data_Transfer_Operator(from, to),
        verbose_(verbose)
    { }

      bool perform_data_transfer() const;

    private:
      bool verbose_;
  };

  // -------------------------------------------------------------------------
  // A Simple Transfer operator for moving data between Mock_DeCART Energy and Mock_DeCART Depletion Model Evaluators
  // -------------------------------------------------------------------------
  class Transfer_Mock_DeCART_Energy_To_Depletion: public LIME::Data_Transfer_Operator 
  {
    public:

      Transfer_Mock_DeCART_Energy_To_Depletion(const Teuchos::RCP<Mock_DeCART_Energy_ME> from, Teuchos::RCP<Mock_DeCART_Depletion_ME> to, bool verbose = false ) :
        LIME::Data_Transfer_Operator(from, to),
        verbose_(verbose)
    { }

      bool perform_data_transfer() const;

    private:
      bool verbose_;
  };

  // -------------------------------------------------------------------------
  // A Simple Transfer operator for moving data between Mock_DeCART Energy and Mock_DeCART Depletion Model Evaluators
  // -------------------------------------------------------------------------
  class Transfer_Mock_DeCART_Depletion_To_Energy: public LIME::Data_Transfer_Operator 
  {
    public:

      Transfer_Mock_DeCART_Depletion_To_Energy(const Teuchos::RCP<Mock_DeCART_Depletion_ME> from, Teuchos::RCP<Mock_DeCART_Energy_ME> to, bool verbose = false ) :
        LIME::Data_Transfer_Operator(from, to),
        verbose_(verbose)
    { }

      bool perform_data_transfer() const;

    private:
      bool verbose_;
  };
}

#endif
