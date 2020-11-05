#ifndef THREE_NONLINEAR_EQUATIONS_W_RESIDUALS_HPP
#define THREE_NONLINEAR_EQUATIONS_W_RESIDUALS_HPP

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

namespace LIME 
{
  namespace Three_Nonlinear_Equations_w_Resid
  {
    // -----------------------------------------------------------------------------------
    // A simple three-equation nonlinear system to test exposed residuals
    //
    //               2*x1   -   x2 + x3^3   - 31 = 0 
    //                 x1   + 2*x2 - 2*x3^2 + 21 = 0 
    //                 x1^2 +   x2 + 2*x3   - 5  = 0 
    //
    // The system has an exact result of x1 = 1, x2 = -2, x3 = 3
    // -----------------------------------------------------------------------------------

    double residual_1(double x1, double x2, double x3, bool verbose);
    double residual_2(double x1, double x2, double x3, bool verbose);
    double residual_3(double x1, double x2, double x3, bool verbose);

    const vector<vector<double> > & get_gold_jacobian_coeffs();

    // -------------------------------------------------------------
    // A Simple Model Evaluator for each equation in a linear system
    // -------------------------------------------------------------
    class EQ_ModelEval : public LIME::Model_Evaluator
    {
      public:

        enum  EQUATION_ID { 
          EQ_ONE = 0,
          EQ_TWO,
          EQ_THREE
        };

        static Teuchos::RCP<Three_Nonlinear_Equations_w_Resid::EQ_ModelEval> create_instance(LIME::Problem_Manager &, EQUATION_ID, bool verbose = false);
        static double (*get_fn(EQUATION_ID))(double, double, double, bool);

        EQ_ModelEval(LIME::Problem_Manager & pm, double (*f_ptr)(double, double, double, bool), EQUATION_ID my_id, bool verbose = false);

        ~EQ_ModelEval() {}

        bool supports_standalone_solve() const { return false; }

        bool supports_residual() const { return true; }

        // Methods needed to support residual fills
        EpetraExt::ModelEvaluator::InArgs createInArgs() const;
        EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
        void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

        Teuchos::RCP<const Epetra_Map> get_x_map() const
        { return epetra_map_; }

        Teuchos::RCP<const Epetra_Map> get_f_map() const
        { return epetra_map_; }

        Teuchos::RCP<const Epetra_Vector> get_x_init() const
        { return soln_vec_; }

        Teuchos::RCP<Epetra_Vector> get_x_state()
        { return soln_vec_; }

        Teuchos::RCP<Epetra_Vector> get_x_state() const
        { return soln_vec_; }

        Teuchos::RCP<Epetra_Map> epetra_map_;
        Teuchos::RCP<Epetra_Vector> soln_vec_;
        void initializeSolution();

        void reinitialize();

        const double & get_my_x() { return (*soln_vec_)[0]; }
        double get_known_x() { return my_known_soln_; }
        EQUATION_ID get_my_eq_id() { return my_eq_id_; }
        void set_other_x(EQUATION_ID eq_id, double x) { other_x_[eq_id] = x; }
        const vector<double> & get_my_iterates() { return iterate_vals_; }

      private:

        bool verbose_;
        EQUATION_ID my_eq_id_;
        mutable std::map<EQUATION_ID, double> other_x_;
        double my_known_soln_;
        double (*my_func_)(double, double, double, bool);
        vector<double> iterate_vals_;
    };

    // -------------------------------------------------------------------------
    // A Simple Transfer operator for moving data between EQ_ModelEval instances
    // -------------------------------------------------------------------------
    class EQ_Transfer : public LIME::Data_Transfer_Operator 
    {
      public:

        EQ_Transfer(Teuchos::RCP<EQ_ModelEval> from, Teuchos::RCP<EQ_ModelEval> to, bool verbose = false) :
          LIME::Data_Transfer_Operator(from, to),
          verbose_(verbose)
      { }

        bool perform_data_transfer() const;

      private:
        bool verbose_;
    };

    std::vector<Teuchos::RCP<LIME::Data_Transfer_Operator> > get_all_transfers(const std::vector<Teuchos::RCP<EQ_ModelEval> > &, bool verbose = false);
  }
}

#endif
