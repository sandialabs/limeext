#ifndef LIME_EXAMPLE_BRUSSELATOR_TEMPERATURE_HPP
#define LIME_EXAMPLE_BRUSSELATOR_TEMPERATURE_HPP

#include "Brusselator_Equation.hpp"
#include "NOX_Epetra.H"

class Brusselator_Temperature : public Brusselator_Equation
{

  public:

    Brusselator_Temperature(LIME::Problem_Manager &, int NumGlobalUnknowns, const string & outputDir = "./", bool verbose=false);

    ~Brusselator_Temperature();

    // Inherited from EpetraExt::ModelEvaluator.
    virtual EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;

    // Inherited from EpetraExt::ModelEvaluator.
    Teuchos::RCP<const Epetra_Map> get_g_map(int) const 
    { return m_g_map; }

    // Set initial condition for solution vector
    virtual void initializeSolution();

    // This allows us to reset our running max_temp if needed
    virtual void finalize_solution(const Epetra_Vector &);

    // A reset function performed at the start of each transient solve (call to integrate)
    virtual void reinitialize();

    // Evaluates the function (F) and/or the Jacobian using the solution 
    // values in solnVector.
    virtual bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
        const Epetra_Vector & solnVector, 
        Epetra_Vector & rhsVector) const ;

    virtual Teuchos::RCP<Epetra_Vector> get_response_vec() const
    { return m_g_vec; }

    // Evaluates a response function for use with TriKota or other
    virtual void compute_problem_response(const Epetra_Vector &) const;

  protected:

    static const int numResponses = 1;
    Teuchos::RCP<Epetra_LocalMap> m_g_map;
    Teuchos::RCP<Epetra_Vector> m_g_vec;
    double m_max_temp;
};
#endif
