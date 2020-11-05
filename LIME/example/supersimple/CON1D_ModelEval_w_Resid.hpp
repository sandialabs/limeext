#ifndef LIME_EXAMPLE_CON1D_W_RESID_MODELEVAL_HPP
#define LIME_EXAMPLE_CON1D_W_RESID_MODELEVAL_HPP

#include "CON1D_ModelEval.hpp"

class CON1D_ModelEval_w_Resid : public CON1D_ModelEval
{
  public:

    CON1D_ModelEval_w_Resid(LIME::Problem_Manager & pm, const string & name);

    virtual ~CON1D_ModelEval_w_Resid() {};

    // Methods in addition to those in CON1D_ModelEval that support residual fill callbacks
    virtual EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    virtual EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    virtual void evalModel(const InArgs&, const OutArgs&) const;
    virtual Teuchos::RCP<const Epetra_Map> get_x_map() const
      { return epetra_map_; }
    virtual Teuchos::RCP<const Epetra_Map> get_f_map() const
      { return epetra_map_; }
    virtual Teuchos::RCP<const Epetra_Vector> get_x_init() const;

    virtual Teuchos::RCP<Epetra_Vector> get_x_state();
    virtual Teuchos::RCP<Epetra_Vector> get_x_state() const;

    void initializeSolution();

    // Disable use of standalone solve using this...
    virtual bool supports_standalone_solve() const { return true; }

    // This should be removed in favor of querying via createInArgs
    virtual bool supports_residual() const { return true; }

  protected:

    // A data layour map needed to create/configure solver objects
    Teuchos::RCP<Epetra_Map> epetra_map_;

    // The ME state vector conformal to the data layout map and that is
    // also used to create/configure solver objects
    mutable Teuchos::RCP<Epetra_Vector> me_interface_soln_;

    // needed to get residual array values
    float* r; 
};
#endif

