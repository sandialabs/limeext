
// Unit test headers
#include <Teuchos_UnitTestHarness.hpp>

// Relevant headers
#include <LIME_Problem_Manager.hpp>
#include <LIME_Model_Evaluator.hpp>
#include <LIME_FixedPoint_Accelerators.hpp>

// Trilinos Objects
#include <Epetra_SerialComm.h>

#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ScalarTraits.hpp"

typedef double ScalarType;

using namespace LIME;

namespace {

  // -------------------------------------------------------------
  // A Simple Model Evaluator for each equation in a linear system
  // -------------------------------------------------------------

  class Test_Accel_ModelEval : public Model_Evaluator 
  {
    public:

      static const int MIXDIM = 6;
      static const int SOLDIM = 5;

      Test_Accel_ModelEval(Problem_Manager & pm, bool verbose = false) :
        Model_Evaluator(pm),
        verbose_(verbose)
      {
        matrix_.reshape(SOLDIM, SOLDIM);
        b_.resize(SOLDIM);
        soln_.resize(SOLDIM);
        work_.resize(SOLDIM);

        matrix_(0,0) = 0.823295 ;      matrix_(0,1) = 0.10794   ;     matrix_(0,2) = 0.904459 ;      matrix_(0,3) = 0.213938 ;      matrix_(0,4) = -0.686642;
        matrix_(1,0) = -0.604897;      matrix_(1,1) = -0.0452059;     matrix_(1,2) = 0.83239  ;      matrix_(1,3) = -0.967399;      matrix_(1,4) = -0.198111;
        matrix_(2,0) = -0.329554;      matrix_(2,1) = 0.257742  ;     matrix_(2,2) = 0.271423 ;      matrix_(2,3) = -0.514226;      matrix_(2,4) = -0.740419;
        matrix_(3,0) = 0.536459 ;      matrix_(3,1) = -0.270431 ;     matrix_(3,2) = 0.434594 ;      matrix_(3,3) = -0.725537;      matrix_(3,4) = -0.782382;
        matrix_(4,0) = -0.444451;      matrix_(4,1) = 0.0268018 ;     matrix_(4,2) = -0.716795;      matrix_(4,3) = 0.608354 ;      matrix_(4,4) = 0.997849 ;

        if( verbose_ )
        {
          cout << "----------------------------------------------\nGenerated matrix :" << endl;
          matrix_.print(cout);
        }

        soln_(0) = 0.823295;
        soln_(1) = 0.10794;
        soln_(2) = 0.904459;
        soln_(3) = 0.213938;
        soln_(4) = -0.686642;
        if( verbose_ )
        {
          cout << "----------------------------------------------\nGold solution :" << endl;
          soln_.print(cout);
        }

        b_ = soln_;
        b_.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -1.0, matrix_, soln_, 1.0);
        if( verbose_ )
        {
          cout << "----------------------------------------------\nRHS :" << endl;
          b_.print(cout);
        }
      }

      ~Test_Accel_ModelEval() {}

      bool solve_standalone( Teuchos::SerialDenseVector<int,ScalarType> & vec )
      {
        work_.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, matrix_, vec, 0.0);
        vec = b_;
        vec += work_;

        return true;
      }

      const Teuchos::SerialDenseVector<int,ScalarType> & get_soln() const
      { return soln_; }

    private:

      bool verbose_;

      Teuchos::SerialDenseMatrix<int,ScalarType> matrix_;
      Teuchos::SerialDenseVector<int,ScalarType> b_;
      Teuchos::SerialDenseVector<int,ScalarType> soln_;
      Teuchos::SerialDenseVector<int,ScalarType> work_;
  };
}


TEUCHOS_UNIT_TEST(UnitTestAndersonAcceleration, testUnit2)
{
  Epetra_SerialComm comm;
  bool verbose = false;
  Problem_Manager pm(comm, verbose);

  Test_Accel_ModelEval model_eval(pm, verbose);

  const int soldim = Test_Accel_ModelEval::SOLDIM;
  const int mixdim = Test_Accel_ModelEval::MIXDIM;

  LIME::Anderson_Acceleration accelerator(mixdim, verbose);

    const Teuchos::SerialDenseVector<int,ScalarType> & x_true = model_eval.get_soln();
    Teuchos::SerialDenseVector<int,ScalarType> x_err(soldim);

    Epetra_BlockMap map(soldim, 1, 0, comm);
    Epetra_Vector xVec(map);
    Teuchos::SerialDenseVector<int,ScalarType> x_vec(soldim);

    accelerator.initialize(xVec);

    cout << endl;
    for( int i = 0; i < 6; ++i )
    {
      model_eval.solve_standalone(x_vec);
      for( int j = 0; j < soldim; ++j )
        xVec[j] = x_vec(j);
      accelerator.compute_update(xVec);
      for( int j = 0; j < soldim; ++j )
        x_vec(j) = xVec[j];
      if( verbose )
      {
        cout << "index = " << (i+1) << endl;
        cout << endl << "-----------------\nx_vec : " << xVec << endl;
      }
      x_err = x_vec;
      x_err.scale(-1.0);
      for( int j = 0; j < soldim; ++j )
        x_err(j) += x_true(j);
      cout << "Error norm = " << sqrt(x_err.dot(x_err)) << endl;
    }

  // Test that we make it to the end
  TEST_ASSERT(true);
}
