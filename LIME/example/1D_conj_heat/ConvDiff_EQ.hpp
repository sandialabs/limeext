#ifndef LIME_EXAMPLE_CONVDIFF_EQ_HPP
#define LIME_EXAMPLE_CONVDIFF_EQ_HPP

#include "LIME_Model_Evaluator.hpp"

class ConvDiff_EQ : public LIME::Model_Evaluator
{
  public:

    enum INTERFACE_LOCATION {
      LEFT    ,
      RIGHT   
    };

    // Constructor
    ConvDiff_EQ(
        LIME::Problem_Manager & pm            ,
        double        peclet                  ,
        double        radiation               ,
        double        kappa                   ,
        double        bcWeight                ,
        INTERFACE_LOCATION                    ,
        double        xmin                    ,
        double        xmax                    ,
        double        Tleft                   ,
        double        Tright                  ,
        int           NumGlobalUnknowns = 0   ,
        string        name = ""                 ) ;

    // Destructor
    ~ConvDiff_EQ() { }

    // Method to compute exact interface temperature
    static double computeAnalyticInterfaceTemp( 
        double Radiation      ,
        double T_left  = 0.0  ,
        double T_right = 1.0  ,
        double kappa   = 0.1  ,
        double peclet  = 9.0      );


    // Compute and output analytic slution here
    void initialize();

    // Evaluates the function (F) and/or the Jacobian using the solution 
    // values in solnVector.
    bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
        const Epetra_Vector & soln, Epetra_Vector & rhs) const;

    // Compute the appropriate end heat flux
    double compute_heat_flux(const Epetra_Vector &) const;

    // Accessor
    double get_interface_temp() const;

    // Accessor
    Teuchos::RCP<Epetra_Vector> getExactSolution()
    { return exactSolution_; }

    // Additional status output
    virtual void outputStatus( ostream & os );

    virtual bool is_transient() const { return false; }

    virtual bool supports_standalone_solve() const { return false; }
    virtual bool supports_residual() const { return true; }


    virtual unsigned int get_max_steps() const { return 50; }

    virtual EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    virtual EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    virtual void evalModel(const InArgs&, const OutArgs&) const;
    virtual Teuchos::RCP<const Epetra_Map> get_x_map() const
    { return epetra_map_; }
    virtual Teuchos::RCP<const Epetra_Map> get_f_map() const
    { return epetra_map_; }
    virtual Teuchos::RCP<const Epetra_Vector> get_x_init() const
    { return me_interface_soln_; }

    virtual Teuchos::RCP<Epetra_Vector> get_x_state()
    { return me_interface_soln_; }
    virtual Teuchos::RCP<Epetra_Vector> get_x_state() const
    { return me_interface_soln_; }

    void initializeSolution();

    void set_xfer_T_val(double val)
    { interface_T_from_xfer_ = val; }

    void set_xfer_Flux_val(double val)
    { interface_Flux_from_xfer_ = val; }

  protected:

    // A data layour map needed to create/configure solver objects
    Teuchos::RCP<Epetra_Map> epetra_map_;

    // The ME state vector conformal to the data layout map and that is
    // also used to create/configure solver objects
    mutable Teuchos::RCP<Epetra_Vector> me_interface_soln_;

  protected:

    double xmin_                  ;
    double xmax_                  ;
    double Tleft_                 ;
    double Tright_                ;
    double T1_exact_              ;

    double peclet_                ;
    double radiation_             ;
    double kappa_                 ;

    INTERFACE_LOCATION iLocation_ ;
    double dirScale_              ;
    double bcWeight_              ;

    Teuchos::RCP<Epetra_Vector> exactSolution_ ;
    Teuchos::RCP<Epetra_Vector> dTdx_          ;

    int interface_elem_;
    int local_node_    ;
    int interface_node_;
    int opposite_node_ ;
    double interface_T_from_xfer_;
    double interface_Flux_from_xfer_;

    Teuchos::RCP<Epetra_Vector> mesh_ ;

};

#include "LIME_Data_Transfer_Operator.hpp"

class ConvDiff_Xfer : public LIME::Data_Transfer_Operator 
{
  public:
    ConvDiff_Xfer(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to) : 
      LIME::Data_Transfer_Operator(from, to)
  { }

    ~ConvDiff_Xfer() { }

    virtual bool perform_data_transfer() const
    {
      ConvDiff_EQ * s =  dynamic_cast<ConvDiff_EQ *>(&(*source()));
      ConvDiff_EQ * t =  dynamic_cast<ConvDiff_EQ *>(&(*target()));

      t->set_xfer_T_val(s->get_interface_temp());
      t->set_xfer_Flux_val(s->compute_heat_flux(*s->get_x_state()));
      //std::cout << "Transferred from \"" << source->get_my_name() << "\" to \"" << target->get_my_name() << "\" : "
      //  << "\tT_int = " << source->get_interface_temp() << ", Flux_int = " << source->compute_heat_flux(*source->get_x_state()) << endl;

      return true;
    }
};

#endif
