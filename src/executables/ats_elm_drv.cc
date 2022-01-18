/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS driver class for ATS-ELM interface

License: see $ATS_DIR/COPYRIGHT
Authors: Ethan Coon, F.-M. Yuan @ ORNL

Implementation (driver) for the ats_elm_interface.  The interface driver is a class to:
(1) read-in inputs by filename passing from ELM, and communicator synchorizing.
(2) pass data from ELM to ATS and setup/initialize ATS.
(3) in ONE single ELM-timestep, run cycle driver, which runs the overall, top level timestep loop.
    - inputs from ELM: starting/ending time in seconds and will reset "end time" in cycle_driver parameter list;
    - option from ELM: reset initial conditions or states from ELM (by default NO)
(4) return data to ELM.

NOTE: Since ELM ususually runs at half-hourly timestep for hundred/thousand years,
      it's a must to turn off all screen/std print-out to avoid large amounts of ascii log file.
      Therefore, set global verbosity level to NONE by default.
------------------------------------------------------------------------- */

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"
#include "exceptions.hh"

#include "boost/filesystem.hpp"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

//
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

//
#include "GeometricModel.hh"
#include "State.hh"
#include "ats_mesh_factory.hh"

//
#include "InputAnalysis.hh"
#include "Units.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"

// registration files
#include "state_evaluators_registration.hh"

#include "ats_relations_registration.hh"
#include "ats_transport_registration.hh"
#include "ats_energy_pks_registration.hh"
#include "ats_energy_relations_registration.hh"
#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_deformation_registration.hh"
#include "ats_bgc_registration.hh"
#include "ats_surface_balance_registration.hh"
#include "ats_mpc_registration.hh"
#include "ats_sediment_transport_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#ifdef ALQUIMIA_ENABLED
#include "pks_chemistry_registration.hh"
#endif

#define DEBUG_MODE 1


// -----------------------------------------------------------------------------
//
// ats_elm interface main class
//
// -----------------------------------------------------------------------------

#include "ats_elm_drv.hh"

namespace ATS {

ats_elm_drv::ats_elm_drv() {

  // create and start the global timer
  timer_ = Teuchos::rcp(new Teuchos::Time("wallclock_monitor",true));

};

ats_elm_drv::~ats_elm_drv() {
   //destructor (TODO)

};

// -----------------------------------------------------------------------------
// setup & initializing
// -----------------------------------------------------------------------------

int ats_elm_drv::drv_init(std::string input_filename) {
		//const Teuchos::RCP<const Amanzi::Comm_type>& comm) {

  //
  input_filename_ = input_filename;

  // STEP 1: configuring MPI, mesh (domain), states, ... from input file passing by ELM

  // -- create/set communicator (TODO)
  comm_ = Amanzi::getDefaultComm();  // TODO: should be sync. with 'comm'
  comm_size = comm_->NumProc();
  comm_rank = comm_->MyPID();

  // STEP 2: read-in a generic .xml from input file passing by ELM

  if (!boost::filesystem::exists(input_filename_)) {
    if (comm_rank == 0) {
	   std::cerr << "ERROR: input file \"" << input_filename_ << "\" does not exist." << std::endl;
	}
	return 1;
  }

  // parse input file
  parameter_list_ = Teuchos::getParametersFromXmlFile(input_filename_);

  // STEP 3: set default verbosity level to NONE
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

  return 0;

}


int ats_elm_drv::drv_setup(const double start_ts, const bool initialoutput) {

  // (1) create meshes, regions, state from input .xml, AND over-ride them from elm if any.

  Teuchos::ParameterList& plist = *parameter_list_;

  // (1a) ELM surface grids/soil columns passing into ats mesh parameter lists
  mesh_parameter_reset();

  // (1b) geometric model and regions
  Teuchos::ParameterList reg_params = plist.sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm_) );

  // (1c) state
  Teuchos::ParameterList state_plist = plist.sublist("state");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));

  // (1d) create and register meshes
  ATS::Mesh::createMeshes(plist, comm_, gm, *S_);

  //mesh_vertices_reset(); //TODO - not yet works

  // (1e) 'cycle driver' or 'coordinator' parameter-list
  ats_elm_drv_list_ = Teuchos::sublist(parameter_list_, "cycle driver");
  cycle_driver_read_parameter();

  // create the top level PK
  Teuchos::ParameterList pk_tree_list = ats_elm_drv_list_->sublist("PK tree");
  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string &pk_name = pk_tree_list.name(pk_item);

  // (2) model setup

  // create the solution
  soln_ = Teuchos::rcp(new Amanzi::TreeVector());

  // create the pk
  Amanzi::PKFactory pk_factory;
  pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, parameter_list_, S_, soln_);

  // create the checkpointing
  Teuchos::ParameterList& chkp_plist = parameter_list_->sublist("checkpoint");
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, *S_));

  // create the observations
  Teuchos::ParameterList& observation_plist = parameter_list_->sublist("observations");
  for (auto& sublist : observation_plist) {
    if (observation_plist.isSublist(sublist.first)) {
      observations_.emplace_back(Teuchos::rcp(new Amanzi::UnstructuredObservations(
                observation_plist.sublist(sublist.first))));
    } else {
      Errors::Message msg("\"observations\" list must only include sublists.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // check whether meshes are deformable, and if so require a nodal position
  for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {

    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      std::string node_key;
      if (mesh->first != "domain") node_key= mesh->first+std::string("-vertex_coordinate");
      else node_key = std::string("vertex_coordinate");

      S_->RequireField(node_key)->SetMesh(mesh->second.first)->SetGhosted()
          ->AddComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension());
    }

    // writes region information
    if (parameter_list_->isSublist("analysis")) {
      Amanzi::InputAnalysis analysis(mesh->second.first, mesh->first);
      analysis.Init(parameter_list_->sublist("analysis").sublist(mesh->first));
      analysis.RegionAnalysis();
      analysis.OutputBCs();
    }
  }

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());

  // verbose object (NOTE: already set option level to NONE)
  vo_ = Teuchos::rcp(new Amanzi::VerboseObject("ats_elm_driver", *parameter_list_));

  // start at time t = t0 and initialize the state.
  double dt_restart = -1;
  {
    //Teuchos::TimeMonitor monitor(*setup_timer_);
    StatePKsetup();
    dt_restart = initialize();

    //
    // ic_reset();

  }

  // if output initial states
  if (initialoutput) {
    double dt = dt_restart > 0 ? dt_restart : get_dt(false);

    // visualization at IC
    visualize(false);
    checkpoint(dt);
  }

  // have to sync starting time with ELM
  t0_ = start_ts;
  ats_elm_drv_list_->set("start time",t0_,"s");

  return 0;
}

// ELM domain, surface-grids/soil profile, passing into ats via 'plist'
void ats_elm_drv::mesh_parameter_reset(const bool elm_matched) {

  /*---------------------------------------------------------------------------------*/
  // (1) mesh, including 'surface' and 'domain'

  Teuchos::RCP<Teuchos::ParameterList> mesh_plist_ = Teuchos::sublist(parameter_list_, "mesh");
  Teuchos::RCP<Teuchos::ParameterList> domain_plist_ = Teuchos::sublist(mesh_plist_, "domain");
  std::cout << "mesh type: " << domain_plist_->get<Teuchos::string>("mesh type") << std:: endl;

  if (domain_plist_->get<Teuchos::string>("mesh type") == "generate mesh") {

	  Teuchos::RCP<Teuchos::ParameterList>coords_plist_ = Teuchos::sublist(domain_plist_, "generate mesh parameters");

	  Teuchos::Array<double> coords_low;
	  Teuchos::Array<double> coords_high;
	  Teuchos::Array<double> coords_point;

	  //surface grids over-ride (TODO)

	  // modify vertical range
	  coords_low = coords_plist_->get<Teuchos::Array<double>>("domain low coordinate");
	  coords_high = coords_plist_->get<Teuchos::Array<double>>("domain high coordinate");

      // let's ats determine mesh, except for box size (ranges)
      coords_low[coords_low.size()-1] = elm_col_nodes[length_nodes-1];
      coords_high[coords_high.size()-1] = elm_col_nodes[0];
      Teuchos::Array<int> ncells = coords_plist_->get<Teuchos::Array<int>>("number of cells");
      ncells[0] = length_gridsX-1;
      ncells[1] = length_gridsY-1;
      ncells[2] = length_nodes-1;
      coords_plist_->set("number of cells", ncells);

      if (elm_matched) {
		//coords = elm_col_nodes;  //incorrect and not here
	  }

      coords_plist_->set("domain low coordinate", coords_low, "m");
	  coords_plist_->set("domain high coordinate", coords_high, "m");

	  /*---------------------------------------------------------------------------------*/
	  // (2) regions, including 'computational domain', 'surface domain', 'surface', 'bottom face', etc.
	  Teuchos::RCP<Teuchos::ParameterList> region_plist_ = Teuchos::sublist(parameter_list_, "regions");

	  Teuchos::RCP<Teuchos::ParameterList> compu_domain_ = Teuchos::sublist(region_plist_, "computational domain");
	  coords_plist_ = Teuchos::sublist(compu_domain_, "region: box");
	  coords_plist_->set("low coordinate", coords_low, "m");    // exactly same as 'mesh->domain'
	  coords_plist_->set("high coordinate", coords_high, "m");  // exactly same as 'mesh->domain'

	  Teuchos::RCP<Teuchos::ParameterList> sideset_surface_ = Teuchos::sublist(region_plist_, "surface");
	  coords_plist_ = Teuchos::sublist(sideset_surface_, "region: plane");
	  coords_point = coords_plist_->get<Teuchos::Array<double>>("point");
	  //coords_point[coords_point.size()-1] = elm_surf_gridsZ[0]; //TODO - this is what really needed
	  coords_point[coords_point.size()-1] = elm_col_nodes[0];
	  coords_plist_->set("point", coords_point, "m");

  }; // 'mesh type' is 'generate mesh'

}

void ats_elm_drv::mesh_vertices_reset(){

  //for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
  //     mesh!=S_->mesh_end(); ++mesh) {

      if (S_->HasMesh("domain")) {
        auto mesh_ = S_->GetMesh("domain");
        int dim = 3;
        Amanzi::AmanziGeometry::Point coords(dim);

        // number of vertices
        int nV = mesh_ -> num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);  // or, 'ALL'?

        // collect coordinates and override Z coords (TODO: X, Y)

        for (int iV=0; iV<nV; iV++) {
          // get the coords of the node
          mesh_ -> node_get_coordinates(iV,&coords);
          auto old = coords;

          // need to known Z indices for current node
          int n_nabvid = 0;
          int nabvid = mesh_->node_get_node_above(iV);
          while (nabvid != -1) {
        	  n_nabvid = n_nabvid + 1;
        	  nabvid = mesh_->node_get_node_above(nabvid);
          }
          coords[dim-1] = elm_col_nodes[n_nabvid];  // NOT [length_nodes-nextid-1] ??? (further checking)

          //mesh_ -> node_set_coordinates(iV, coords);
          // S_->GetMesh("domain") NON-changeable?? - (TODO) needs a new thought here


          std::cout<< "checking vertice " << iV <<" - iZ "<<n_nabvid<<" ; ";
          std::cout << "old coords: "<< old << " - reset: " << coords;
          std::cout<< std::endl;

        }
      }

   //}

}

void ats_elm_drv::cycle_driver_read_parameter() {
  Amanzi::Utils::Units units;
  t0_ = ats_elm_drv_list_->get<double>("start time");
  std::string t0_units = ats_elm_drv_list_->get<std::string>("start time units", "s");
  if (!units.IsValidTime(t0_units)) {
    Errors::Message msg;
    msg << "ats_elm_drv start time: unknown time units type: \"" << t0_units << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  bool success;
  t0_ = units.ConvertTime(t0_, t0_units, "s", success);

  t1_ = ats_elm_drv_list_->get<double>("end time");
  std::string t1_units = ats_elm_drv_list_->get<std::string>("end time units", "s");
  if (!units.IsValidTime(t1_units)) {
    Errors::Message msg;
    msg << "ats_elm_drv end time: unknown time units type: \"" << t1_units << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  t1_ = units.ConvertTime(t1_, t1_units, "s", success);

  max_dt_ = ats_elm_drv_list_->get<double>("max time step size [s]", 1.0e99);
  min_dt_ = ats_elm_drv_list_->get<double>("min time step size [s]", 1.0e-12);
  cycle0_ = ats_elm_drv_list_->get<int>("start cycle",0);
  cycle1_ = ats_elm_drv_list_->get<int>("end cycle",-1);
  duration_ = ats_elm_drv_list_->get<double>("wallclock duration [hrs]", -1.0);
  subcycled_ts_ = ats_elm_drv_list_->get<bool>("subcycled timestep", false);

  // restart control
  restart_ = ats_elm_drv_list_->isParameter("restart from checkpoint file");
  if (restart_) restart_filename_ = ats_elm_drv_list_->get<std::string>("restart from checkpoint file");
}



void ats_elm_drv::StatePKsetup() {
  // Set up the states, creating all data structures.
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  S_->RequireScalar("dt", "ats_elm_drv");

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->Setup(S_.ptr());
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();
}

double ats_elm_drv::initialize() {
  Teuchos::OSTab tab = vo_->getOSTab();

  // Restart from checkpoint part 1:
  //  - get the time prior to initializing anything else
  if (restart_) {
    S_->set_time(Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_));
  }

  // Initialize the state
  *S_->GetScalarData("dt", "ats_elm_drv") = 0.;
  S_->GetField("dt","ats_elm_drv")->set_initialized();
  S_->InitializeFields();

  // Initialize the process kernels
  pk_->Initialize(S_.ptr());

  // Restart from checkpoint part 2:
  // -- load all other data
  double dt_restart = -1;
  if (restart_) {
    dt_restart = Amanzi::ReadCheckpoint(*S_, restart_filename_);
    t0_ = S_->time();
    cycle0_ = S_->cycle();

    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first)) {
        Amanzi::DeformCheckpointMesh(*S_, mesh->first);
      }
    }
  }

  // Final checks.
  S_->CheckNotEvaluatedFieldsInitialized();
  S_->InitializeEvaluators();
  S_->InitializeFieldCopies();
  S_->CheckAllFieldsInitialized();

  // commit the initial conditions.
  pk_->CommitStep(S_->time(), S_->time(), S_);

  // Write dependency graph.
  S_->WriteDependencyGraph();
  S_->InitializeIOFlags();

  // Check final initialization
  WriteStateStatistics(*S_, *vo_);

  // Set up visualization
  auto vis_list = Teuchos::sublist(parameter_list_,"visualization");
  for (auto& entry : *vis_list) {
    std::string domain_name = entry.first;

    if (S_->HasMesh(domain_name)) {
      // visualize standard domain
      auto mesh_p = S_->GetMesh(domain_name);
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);
      if (!sublist_p->isParameter("file name base")) {
        if (domain_name.empty() || domain_name == "domain") {
          sublist_p->set<std::string>("file name base", std::string("ats_vis"));
        } else {
          sublist_p->set<std::string>("file name base", std::string("ats_vis_")+domain_name);
        }
      }

      if (S_->HasMesh(domain_name+"_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
        mesh_p = S_->GetMesh(domain_name+"_3d");

      // vis successful timesteps
      auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
      vis->set_name(domain_name);
      vis->set_mesh(mesh_p);
      vis->CreateFiles(false);

      visualization_.push_back(vis);

    } else if (Amanzi::Keys::isDomainSet(domain_name)) {
      // visualize domain set
      const auto& dset = S_->GetDomainSet(Amanzi::Keys::getDomainSetName(domain_name));
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);

      if (sublist_p->get("visualize individually", false)) {
        // visualize each subdomain
        for (const auto& subdomain : *dset) {
          Teuchos::ParameterList sublist = vis_list->sublist(subdomain);
          sublist.set<std::string>("file name base", std::string("ats_vis_")+subdomain);
          auto vis = Teuchos::rcp(new Amanzi::Visualization(sublist));
          vis->set_name(subdomain);
          vis->set_mesh(S_->GetMesh(subdomain));
          vis->CreateFiles(false);
          visualization_.push_back(vis);
        }
      } else {
        // visualize collectively
        auto domain_name_base = Amanzi::Keys::getDomainSetName(domain_name);
        if (!sublist_p->isParameter("file name base"))
          sublist_p->set("file name base", std::string("ats_vis_")+domain_name_base);
        auto vis = Teuchos::rcp(new Amanzi::VisualizationDomainSet(*sublist_p));
        vis->set_name(domain_name_base);
        vis->set_mesh(dset->get_referencing_parent());
        for (const auto& subdomain : *dset) {
          vis->set_subdomain_mesh(subdomain, S_->GetMesh(subdomain));
        }
        vis->CreateFiles(false);
        visualization_.push_back(vis);
      }
    }
  }

  // make observations at time 0
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(cycle0_);

  // set up the TSM
  // -- register visualization times
  for (const auto& vis : visualization_) vis->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register observation times
  for (const auto& obs : observations_) obs->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

  // -- register any intermediate requested times
  if (ats_elm_drv_list_->isSublist("required times")) {
    Teuchos::ParameterList& sublist = ats_elm_drv_list_->sublist("required times");
    Amanzi::IOEvent pause_times(sublist);
    pause_times.RegisterWithTimeStepManager(tsm_.ptr());
  }

  // Create an intermediate state that will store the updated solution until
  // we know it has succeeded.
  S_next_ = Teuchos::rcp(new Amanzi::State(*S_));
  *S_next_ = *S_;
  if (subcycled_ts_) {
    S_inter_ = Teuchos::rcp(new Amanzi::State(*S_));
    *S_inter_ = *S_;
  } else {
    S_inter_ = S_;
  }

  // set the states in the PKs Passing null for S_ allows for safer subcycling
  // -- PKs can't use it, so it is guaranteed to be pristinely the old
  // timestep.  This code is useful for testing that PKs don't use S.
  // pk_->set_states(Teuchos::null, S_inter_, S_next_);

  // That said, S is required to be valid, since it is valid for all time
  // (e.g. from construction), whereas S_inter and S_next are only valid after
  // set_states() is called.  This allows for standard interfaces.
  pk_->set_states(S_, S_inter_, S_next_);

  return dt_restart;
}

// reset ats initial conditions (IC)
void ats_elm_drv::ic_reset() {
  // Three (3) types of ICs
  // (1) constants, in State
  // (2)

  // (3) real IC, i.e. primary variable in PKs
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  std::string pk_name_ = "flow";
  std::string pv_key = flow_plist_->get<std::string>("primary variable key");
  std::cout<<"flow pv_key: "<< pv_key <<std::endl;

  if (S_->HasField(pv_key)){
	Teuchos::RCP<Amanzi::Field> field = S_->GetField(pv_key, pk_name_);
	std::cout <<"state field: "<< field->fieldname()<<" - type: "<< field->type()<<std::endl;
    std::cout << "data: " << *(S_->GetFieldData(pv_key)->ViewComponent("cell")) <<std::endl;

    auto mesh_ = S_->GetMesh("domain");
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    auto &pc = *(S_->GetFieldData(pv_key)->ViewComponent("cell"));
    for (int c = 0; c < ncells; ++c) {
      const auto& xyzc = mesh_->cell_centroid(c);
      std::cout <<"coords: "<<xyzc<<" - "<<pc[0][c] << " - "<< soilp[c]<<std::endl;
      pc[0][c] = soilp[c];
    }

    std::cout << "data 2: " << *(S_->GetFieldData(pv_key)->ViewComponent("cell")) <<std::endl;

  }

}

// reset ats initial conditions (IC)
void ats_elm_drv::bc_reset() {

  // TODO
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  std::string pk_name_ = "flow";

}

// reset ats source-sink terms (SS)
void ats_elm_drv::ss_reset() {

  //
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  std::string pk_name_ = "flow";
  std::string pv_key = flow_plist_->get<std::string>("primary variable key");
  std::cout<<"flow pv_key: "<< pv_key <<std::endl;

  for (auto f_it = S_->field_begin(); f_it != S_->field_end(); ++f_it) {
    std::string name(f_it->first);
    std::cout<<"checking state field: " <<name <<std::endl;
  }

  if (S_->HasField(pv_key)){
	Teuchos::RCP<Amanzi::Field> field = S_->GetField(pv_key, pk_name_);
	std::cout <<"state field: "<< field->fieldname()<<" - type: "<< field->type()<<std::endl;
    std::cout << "data: " << *(S_->GetFieldData(pv_key)->ViewComponent("cell")) <<std::endl;

    auto mesh_ = S_->GetMesh("domain");
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    auto &pc = *(S_->GetFieldData(pv_key)->ViewComponent("cell"));
    for (int c = 0; c < ncells; ++c) {
      const auto& xyzc = mesh_->cell_centroid(c);
      std::cout <<"coords: "<<xyzc<<" - "<<pc[0][c] << " - "<< soilp[c]<<std::endl;
      //pc[0][c] = soilp[c];
    }

    std::cout << "data 2: " << *(S_->GetFieldData(pv_key)->ViewComponent("cell")) <<std::endl;

  }

}


// -----------------------------------------------------------------------------
// ONE single ELM-timestep
// -----------------------------------------------------------------------------

void ats_elm_drv::cycle_driver(const double start_ts, const double end_ts, const bool resetIC_from_elm) {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // staring/ending time (seconds) from ELM, usually one ELM-timestep
  double ts = end_ts - start_ts;
  t0_ = start_ts;
  t1_ = end_ts;
  max_dt_ = ts;
  ats_elm_drv_list_->set("end time",t1_,"s");

  //
  if (resetIC_from_elm) {

	  // reset initial states
	  // TODO

  }

  // get the intial timestep
  double dt = get_dt(false);

  // iterate process kernels
  {

    bool fail = false;

    std::cout << "INFO: ATS runs start at " << S_->time() << std::endl;

    while ((S_->time() < t1_) &&
           ((cycle1_ == -1) || (S_->cycle() <= cycle1_)) &&
           (duration_ < 0 || timer_->totalElapsedTime(true) < duration) &&
           dt > 0.) {
      *S_->GetScalarData("dt", "ats_elm_drv") = dt;
      *S_inter_->GetScalarData("dt", "ats_elm_drv") = dt;
      *S_next_->GetScalarData("dt", "ats_elm_drv") = dt;

      S_->set_initial_time(S_->time());
      S_->set_final_time(S_->time() + dt);
      S_->set_intermediate_time(S_->time());

      fail = advance(S_->time(), S_->time() + dt, dt);

    } // while not finished

    if (not fail) {std::cout << "INFO: ATS runs end at " << S_->time() << std::endl;}

  }

  //
  //visualize(false);
  //checkpoint(dt);
  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics(S_next_);
  checkpoint_->Write(*S_next_, 0.0, true);

}


bool ats_elm_drv::advance(double t_old, double t_new, double& dt_next) {
  double dt = t_new - t_old;

  S_next_->advance_time(dt);
  bool fail = pk_->AdvanceStep(t_old, t_new, false);
  fail |= !pk_->ValidStep();

  // advance the iteration count and timestep size
  S_next_->advance_cycle();

  if (!fail) {
    // commit the state
    pk_->CommitStep(t_old, t_new, S_next_);

    // make observations, vis, and checkpoints
    for (const auto& obs : observations_) obs->MakeObservations(S_next_.ptr());
    visualize();
    dt_next = get_dt(fail);
    checkpoint(dt_next); // checkpoint with the new dt

    // we're done with this time step, copy the state
    *S_ = *S_next_;
    *S_inter_ = *S_next_;

  } else {
    // Failed the timestep.
    // Potentially write out failed timestep for debugging
    for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_next_);

    // The timestep sizes have been updated, so copy back old soln and try again.
    *S_next_ = *S_;
    *S_inter_ = *S_;

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        // collect the old coordinates
        std::string node_key;
        if (mesh->first != "domain") node_key= mesh->first+std::string("-vertex_coordinate");
        else node_key = std::string("vertex_coordinate");

        Teuchos::RCP<const Amanzi::CompositeVector> vc_vec = S_->GetFieldData(node_key);
        vc_vec->ScatterMasterToGhosted();
        const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);

        std::vector<int> node_ids(vc.MyLength());
        Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
        for (int n=0; n!=vc.MyLength(); ++n) {
          node_ids[n] = n;
          if (mesh->second.first->space_dimension() == 2) {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n]);
          } else {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n], vc[2][n]);
          }
        }

        // undeform the mesh
        Amanzi::AmanziGeometry::Point_List final_positions;
        mesh->second.first->deform(node_ids, old_positions, false, &final_positions);
      }
    }
    dt_next = get_dt(fail);
  }
  return fail;
}

double ats_elm_drv::get_dt(bool after_fail) {
  // get the physical step size
  double dt = pk_->get_dt();
  double dt_pk = dt;
  if (dt < 0.) return dt;

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("ats_elm_drv: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) {
    dt = max_dt_;
  }

  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_next_->time(), dt, after_fail);
  if (subcycled_ts_) dt = std::min(dt, dt_pk);
  return dt;
}

// -----------------------------------------------------------------------------

void ats_elm_drv::visualize(bool force) {
  // write visualization if requested
  bool dump = force;
  if (!dump) {
    for (const auto& vis : visualization_) {
      if (vis->DumpRequested(S_next_->cycle(), S_next_->time())) {
        dump = true;
      }
    }
  }

  if (dump) {
    pk_->CalculateDiagnostics(S_next_);
  }

  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(S_next_->cycle(), S_next_->time())) {
      WriteVis(*vis, *S_next_);
    }
  }
}

void ats_elm_drv::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    checkpoint_->Write(*S_next_, dt);
  }
}

// -----------------------------------------------------------------------------

void ats_elm_drv::finalize() {
  // finalizing simulation
  WriteStateStatistics(*S_, *vo_);

  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics(S_next_);
  checkpoint_->Write(*S_next_, 0.0, true);

  // flush observations to make sure they are saved
  for (const auto& obs : observations_) obs->Flush();
}

} // close namespace ATS

