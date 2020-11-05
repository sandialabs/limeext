#ifndef FIVE_LINEAR_EQUATIONS_FIXTURE_HPP
#define FIVE_LINEAR_EQUATIONS_FIXTURE_HPP

// Relevant headers
#include <LIME_Hierarchy_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_Data_Transfer_Operator.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

const double FP_UNIT_TEST_TOL = 1.e-4;

namespace LIME 
{
  namespace Five_Linear_Equations 
  {
    // -----------------------------------------------------------------------------------
    // A Five-equation linear system to various algorithms
    // 
    //    This is a contrived linear system broken into five standalone equations (ME s):
    // 
    //        2.0*x1 -     x2 +     x3                   =   7.0
    //            x1 + 2.0*x2 -     x3                   = - 6.0
    //            x1 +     x2 + 2.0*x3 -              x5 =   6.0
    //            x1          -     x3 +     x4          =   0.0
    //                 2.0*x2 -     x3 -     x4 +     x5 = -10.0
    // 
    //  which should produce the solution (1.0, -2.0, 3.0, 2.0, -1.0) to within tolerance.
    // -----------------------------------------------------------------------------------

    double solve_x1(double x2, double x3, double,    bool verbose);
    double solve_x2(double x1, double x3, double,    bool verbose);
    double solve_x3(double x1, double x2, double x5, bool verbose);
    double solve_x4(double x1, double x3, double,    bool verbose);
    double solve_x5(double x2, double x3, double x4, bool verbose);

    // -------------------------------------------------------------
    // A Simple Model Evaluator for each equation in a linear system
    // -------------------------------------------------------------
    class EQ_ModelEval : public LIME::Model_Evaluator
    {
      public:

        enum  EQUATION_ID { 
          EQ_ONE = 0,
          EQ_TWO,
          EQ_THREE,
          EQ_FOUR,
          EQ_FIVE
        };

        static Teuchos::RCP<Five_Linear_Equations::EQ_ModelEval> create_instance(LIME::Problem_Manager &, EQUATION_ID, bool verbose = false);
        static double (*get_fn(EQUATION_ID))(double, double, double, bool);

        EQ_ModelEval(LIME::Problem_Manager & pm, double (*f_ptr)(double, double, double, bool), EQUATION_ID my_id, bool verbose = false) :
          LIME::Model_Evaluator(pm, "EQ_ModelEval"),
          verbose_(verbose),
          my_eq_id_(my_id),
          my_func_(f_ptr)
      {
        initializeSolution();
      }

        ~EQ_ModelEval() {}
        bool supports_standalone_solve() const { return true; }

        void initializeSolution();

        void reinitialize();

        void solve_standalone();

        bool is_converged();

        const double & get_my_x() { return my_x_; }
        EQUATION_ID get_my_eq_id() { return my_eq_id_; }
        void set_other_x(EQUATION_ID eq_id, double x) { other_x_[eq_id] = x; }
        const vector<double> & get_my_iterates() { return iterate_vals_; }

      private:

        bool verbose_;
        EQUATION_ID my_eq_id_;
        double my_x_;
        std::map<EQUATION_ID, double> other_x_;
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
