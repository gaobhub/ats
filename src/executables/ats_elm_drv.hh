/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* Simulation controller and top-level driver For ONE E3SM Land Model (ELM) Timestep

 (1) The interface imitates ats 'main.cc', 'simulation_driver.cc', 'coordinator.cc' very much, WITHOUT heavy print-out info
 So the spec. are exactly same.

 (2) The building system will compile and link into a package of dynamic libraries,
 with libats_elm_interface.dylib as hooking class to ELM's EMI-EM_ATS interface.

 (3) ALL inputs, inc. *.xml, required by ATS would be in E3SM's input data directory under 'lnd/clm2/ats'.
     NOTE: all parameter list from *.xml ARE merely a data placeholder, by which ELM may reset or override.
          for an example, not matter what values of 'start time' and 'end time' in cycle_drvier,
          ELM will over-ride them each ELM timestep.

 Authors: F.-M. Yuan (yuanf@ornl.gov), Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

In the `"cycle driver`" sublist, the user specifies global control of the
simulation, including starting and ending times and restart options.

.. _coordinator-spec:
.. admonition:: coordinator-spec

    * `"start time`" ``[double]`` **0.** Specifies the start of time in model time.
    * `"start time units`" ``[string]`` **"s"** One of "s", "d", or "yr"

    ONE OF

    * `"end time`" ``[double]`` Specifies the end of the simulation in model time.
    * `"end time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

    OR

    * `"end cycle`" ``[int]`` **optional** If provided, specifies the end of the
      simulation in timestep cycles.

      END

    * `"subcycled timestep`" ``[bool]`` **false**  If true, this coordinator creates
      a third State object to store intermediate solutions, allowing for failed
      steps.
    * `"restart from checkpoint file`" ``[string]`` **optional** If provided,
      specifies a path to the checkpoint file to continue a stopped simulation.
    * `"wallclock duration [hrs]`" ``[double]`` **optional** After this time, the
      simulation will checkpoint and end.
    * `"required times`" ``[io-event-spec]`` **optional** An IOEvent_ spec that
      sets a collection of times/cycles at which the simulation is guaranteed to
      hit exactly.  This is useful for situations such as where data is provided at
      a regular interval, and interpolation error related to that data is to be
      minimized.
    * `"PK tree`" ``[pk-typed-spec-list]`` List of length one, the top level
      PK_ spec.

Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

Example:

.. code-block:: xml

   <ParameterList name="cycle driver">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
     <ParameterList name="PK tree">
       <ParameterList name="my richards pk">
         <Parameter name="PK type" type="string" value="richards" />
       </ParameterList>
     </ParameterList>
   </ParameterList>

*/
/*
 ELM-EMI calling subroutines include (NOT YET finalized):
  subroutine EM_ATS_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
  subroutine set_material_properties(this, l2e_init_list, bounds_clump)
  subroutine set_initial_conditions(this, l2e_init_list, bounds_clump)
  subroutine set_boundary_conditions(this, l2e_init_list, bounds_clump)
  subroutine EM_ATS_OneStep(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, bounds_clump)
  subroutine get_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
  subroutine EM_ATS_Finalize(this, dt, nstep, clump_rank, l2e_list, e2l_list, bounds_clump)

*/

#ifndef ATS_ELM_HH_
#define ATS_ELM_HH_

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "VerboseObject.hh"


namespace Amanzi {
  class TimeStepManager;
  class Visualization;
  class Checkpoint;
  class State;
  class TreeVector;
  class PK;
  class PK_ATS;
  class UnstructuredObservations;
};


namespace ATS {

class ats_elm_drv {

public:
  ats_elm_drv();
  ~ats_elm_drv();

  // some variables in common
  int comm_size;
  int comm_rank;
  std::string input_filename_;
  Teuchos::RCP<Teuchos::ParameterList> parameter_list_;
  Teuchos::RCP<Teuchos::ParameterList> ats_elm_drv_list_;

  // data from ELM
  double* elm_surf_gridsX;  // elm surface grids X coord in meters converted from lon, size of at least 2 (for 1 grid)
  double* elm_surf_gridsY;  // elm surface grids Y coord in meters converted from lat, size of at least 2 (for 1 grid)
  double* elm_surf_gridsZ;  // elm surface grids center elevation in meters, size of at least 1 (for 1 grid)
  double* elm_col_nodes;    // elm soil column nodes in meters, size of 16 (for 15 layers) from top to bottom with upward positive;
  int length_gridsX;
  int length_gridsY;
  int length_nodes;

  double* patm;             // elm atm. air pressure, unit: Pa, 1-D (surf-grids)
  double* soilp;            // elm soil column hydraulic pressure, unit: Pa, 2-D (surf-grids, soil col. layers)
  double* wtd;              // elm soil column water table depth, unit: m, 1-D (surf-grids)

  // PK container and factory
  Teuchos::RCP<Amanzi::PK> pk_;

  // states
  Teuchos::RCP<Amanzi::State> S_;

  // time step manager
  Teuchos::RCP<Amanzi::TimeStepManager> tsm_;

  // primary functions
  int drv_init(std::string input_filename); //,  // input .xml file passing from ELM
		  //const Teuchos::RCP<const Amanzi::Comm_type>& comm); // communicator passing from ELM

  int drv_setup(const double start_ts, const bool initialoutput=false);

  void ic_reset();

  void bc_reset();

  void ss_reset();

  void cycle_driver(const double start_ts, const double end_ts,
		  const bool resetIC);

  void finalize();

private:
  //
  void mesh_parameter_reset(const bool elm_matched=false);  // elm_matched: exactly matched with elm-domain; otherwise ranges only
  void mesh_vertices_reset();
  void cycle_driver_read_parameter();

  // PK methods to be called
  void StatePKsetup();
  double initialize();

  bool advance(double t_old, double t_new, double& dt_next);
  double get_dt(bool after_fail=false);
  void visualize(bool force=false);
  void checkpoint(double dt, bool force=false);
  Teuchos::RCP<Amanzi::State> get_next_state() { return S_next_; }

  Teuchos::RCP<Amanzi::State> S_inter_;
  Teuchos::RCP<Amanzi::State> S_next_;

  Teuchos::RCP<Amanzi::TreeVector> soln_;

  double t0_, t1_;
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;

  // Epetra communicator
  Amanzi::Comm_ptr_type comm_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Amanzi::Visualization> > visualization_;
  std::vector<Teuchos::RCP<Amanzi::Visualization> > failed_visualization_;
  Teuchos::RCP<Amanzi::Checkpoint> checkpoint_;
  bool restart_;
  std::string restart_filename_;

  // observations
  std::vector<Teuchos::RCP<Amanzi::UnstructuredObservations>> observations_;

  // timers
  Teuchos::RCP<Teuchos::Time> timer_;
  double duration_;
  bool subcycled_ts_;

  // fancy OS
  Teuchos::RCP<Amanzi::VerboseObject> vo_;

}; // close of class ats_elm_interface

} //close of namespace ATS

#endif
