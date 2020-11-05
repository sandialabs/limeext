#ifndef _NOX_EXAMPLE_BRUSSELATOR_SPECIES_B_H
#define _NOX_EXAMPLE_BRUSSELATOR_SPECIES_B_H

#include "Brusselator_Equation.hpp"
#include "NOX_Epetra.H"

class Brusselator_Species : public Brusselator_Equation
{

  public:

    Brusselator_Species(LIME::Problem_Manager &, int NumGlobalUnknowns, const string & outputDir = "./", bool verbose=false);

    virtual ~Brusselator_Species();

    // Inherited from EpetraExt::ModelEvaluator.
    virtual EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;

    // Initialize based on registrations
    virtual void initialize();

    // Set initial condition for solution vector
    void initializeSolution();

    // Evaluates the function (F) and/or the Jacobian using the solution 
    // values in solnVector.
    bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
        const Epetra_Vector & solnVector, 
        Epetra_Vector & rhsVector) const;

};
#endif
