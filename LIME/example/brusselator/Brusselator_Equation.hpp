#ifndef LIME_EXAMPLE_BRUSSELATOR_EQUATION_HPP
#define LIME_EXAMPLE_BRUSSELATOR_EQUATION_HPP

#include "LIME_Model_Evaluator.hpp"
#include "NOX_Epetra.H"
#include "Epetra_LocalMap.h"
#include "Teuchos_Array.hpp"

class Brusselator_Equation : public LIME::Model_Evaluator
{

  public:

    Brusselator_Equation(LIME::Problem_Manager &, const string & name, int NumGlobalUnknowns, const string & outputDir, bool verbose=false);

    virtual ~Brusselator_Equation();

    // Inherited from EpetraExt::ModelEvaluator.
    virtual EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    virtual void evalModel(const InArgs&, const OutArgs&) const;
    virtual Teuchos::RCP<const Epetra_Map> get_x_map() const;
    virtual Teuchos::RCP<const Epetra_Map> get_f_map() const;
    virtual Teuchos::RCP<const Epetra_Map> get_g_map(int) const;
    virtual Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    virtual Teuchos::RCP<const Epetra_Vector> get_p_init(int) const;
    virtual Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;

    virtual Teuchos::RCP<Epetra_Vector> get_x_state()
    { return m_solution; }
    virtual Teuchos::RCP<Epetra_Vector> get_x_state() const
    { return m_solution; }

    virtual Teuchos::RCP<Epetra_Vector> get_p_state()
    { return m_p_vec; }
    virtual Teuchos::RCP<Epetra_Vector> get_p_state() const
    { return m_p_vec; }

    virtual void set_xfer_data(const Teuchos::RCP<Epetra_Vector> &);

    // Signal that we can provide a residaul for various uses
    virtual bool supports_residual() const { return true; }

    // Set initial condition for solution vector
    virtual void initializeSolution();

    // Initialize based on registrations
    virtual void initialize();

    // A reset function performed at the start of each transient solve (call to integrate)
    virtual void reinitialize();

    // Evaluates the function (F) and/or the Jacobian using the solution 
    // values in solnVector.
    virtual bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
        const Epetra_Vector & solnVector, 
        Epetra_Vector & rhsVector) const = 0;

    // Evaluates a response function for use with TriKota or other
    virtual void compute_problem_response(const Epetra_Vector & ) const
      { return; }

    virtual bool is_transient() const { return true; }
    virtual double get_max_time() const { return 1.e12; }
    virtual double get_time_step() const { return m_dt; }
    virtual double get_current_time() const { return m_time; }
    virtual unsigned int get_max_steps() const { return m_max_time_steps; }

    // set time
    virtual void update_time() { m_time += m_dt; m_time_step++; }

    // set time step
    virtual void set_time_step(double dt) { m_dt = dt;}

    // Reset problem for next parameter (time) step.
    // For now, this simply updates oldsoln with the given Epetra_Vector

    virtual void finalize_solution(const Epetra_Vector &);

  protected:

    // inserts the global column indices into the Graph
    virtual void generateGraph();

    // writes a vector with its mesh to a text file
    virtual void write_mesh_and_vector_to_text(const Epetra_Vector& x, const string & filename);

  protected:

    const string m_output_dir;
    bool verbose_;

    double m_xmin;
    double m_xmax;
    double m_dx;
    double m_dt;
    double m_time;
    int m_max_time_steps;
    int m_time_step;

    Epetra_Vector *rhs;

    bool useConvection;

    Teuchos::RCP<Epetra_Map> m_f_epetra_map, m_x_epetra_map;
    Teuchos::RCP<Epetra_Vector> m_mesh_x;
    mutable Teuchos::RCP<Epetra_Vector> m_solution;
    mutable Teuchos::RCP<Epetra_Vector> m_old_solution;
    Teuchos::RCP<Epetra_Vector> m_aux_data;

    // All (both) Brusselator equations share the same input params
    static const int numParams = 3;
    Teuchos::RCP<Epetra_LocalMap> m_p_map;
    Teuchos::RCP<Epetra_Vector> m_p_vec;
    Teuchos::RCP<Teuchos::Array<std::string> > m_param_names;
};

#include "LIME_Problem_Manager.hpp"
#include "LIME_Data_Transfer_Operator.hpp"

class Brusselator_EQ_Xfer : public LIME::Data_Transfer_Operator 
{
  public:
    Brusselator_EQ_Xfer(Teuchos::RCP<LIME::Model_Evaluator> from, Teuchos::RCP<LIME::Model_Evaluator> to) : 
      LIME::Data_Transfer_Operator(from, to)
  { }

    ~Brusselator_EQ_Xfer() { }

    virtual bool perform_data_transfer() const
    {
      Brusselator_Equation * s =  dynamic_cast<Brusselator_Equation *>(&(*source()));
      Brusselator_Equation * t =  dynamic_cast<Brusselator_Equation *>(&(*target()));

      t->set_xfer_data(s->get_x_state());

      return true;
    }
};

#endif
