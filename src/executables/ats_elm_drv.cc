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
#include "StateDefs.hh"
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

// passing some of ELM data by resetting 'parameter list' in *.xml
// but have to do this prior to model setup
void ats_elm_drv::plist_reset(const double start_ts) {

  // (0a) ELM surface grids/soil columns passing into ats mesh parameter lists
  //     if mesh type is 'generate mesh' from box-range
  plist_general_mesh_reset();

  // (0b) 'cycle driver' or 'coordinator' parameter-list
  ats_elm_drv_plist_ = Teuchos::sublist(parameter_list_, "cycle driver");
  plist_cycle_driver_reset(start_ts);

  // (0c) parameter lists of 'pks' and 'state' for resetting from ELM
  ats_elm_pks_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  ats_elm_state_plist_ = Teuchos::sublist(parameter_list_, "state");
  plist_material_reset();

}

// ELM domain, surface-grids/soil profile, passing into ats via 'plist'
void ats_elm_drv::plist_general_mesh_reset(const bool elm_matched) {

  /*---------------------------------------------------------------------------------*/
  // (1) mesh, including 'surface' and 'domain'

  Teuchos::RCP<Teuchos::ParameterList> mesh_plist_ = Teuchos::sublist(parameter_list_, "mesh");
  Teuchos::RCP<Teuchos::ParameterList> domain_plist_ = Teuchos::sublist(mesh_plist_, "domain");

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

// override whatever read-in from *.xml
void ats_elm_drv::plist_material_reset() {

  // porosity, viscosity, & permibility in 'state -> field evaluators'
  Teuchos::RCP<Teuchos::ParameterList> field_evals_plist_ = Teuchos::sublist(ats_elm_state_plist_, "field evaluators");
  Teuchos::RCP<Teuchos::ParameterList> state_constants_plist_ = Teuchos::sublist(ats_elm_state_plist_, "initial conditions");

  // (1a) porosity
  auto base_poro_plist_ = Teuchos::sublist(
		  Teuchos::sublist(field_evals_plist_, "base_porosity"), "function");
  int c=0;
  for (auto& entry : *base_poro_plist_) {
    std::string layer_name = entry.first;

    Teuchos::RCP<Teuchos::ParameterList> layer_plist_ = Teuchos::sublist(
		                  Teuchos::sublist(
		                    Teuchos::sublist(base_poro_plist_, layer_name),"function"),
					    "function-constant");
    double poro = layer_plist_->get<double>("value");
    poro = porosity[c];
    layer_plist_->set("value", poro);
    c++;
  }

  // (1b) permeability
  auto visco_plist_ = Teuchos::sublist(field_evals_plist_, "viscosity_liquid");
  double visco = visco_plist_->get<double>("value");
  //std::cout<<"viscosity: "<< visco<<std::endl;

  auto den_mass_plist_ = Teuchos::sublist(field_evals_plist_, "mass_density_liquid");
  double den_mass = den_mass_plist_->get<double>("value");
  //std::cout<<"water mass density: "<< den_mass<<std::endl;

  auto gravity_plist_ = Teuchos::sublist(state_constants_plist_, "gravity");
  double gravity = gravity_plist_->get<Teuchos::Array<double>>("value")[2];
  //std::cout<<"gravity: "<< gravity<<std::endl;

  auto perm_plist_ = Teuchos::sublist(
		  Teuchos::sublist(field_evals_plist_, "permeability"), "function");
  c=0;
  for (auto& entry : *perm_plist_) {
	std::string layer_name = entry.first;

    Teuchos::RCP<Teuchos::ParameterList> layer2_plist_ = Teuchos::sublist(
		                  Teuchos::sublist(
		                    Teuchos::sublist(perm_plist_, layer_name),"function"),
					    "function-constant");
    double perm = layer2_plist_->get<double>("value");
    //std::cout<<"permeablity checking: "<<perm <<" - "
    //		<<hksat[c]*visco/den_mass/(-gravity) <<std::endl;
    perm = hksat[c]*visco/den_mass/(-gravity);   // need to check unit coversion here
    layer2_plist_->set("value", perm);
    c++;
  }


  // water retention curve models, in 'pks -> flow'
  Teuchos::RCP<Teuchos::ParameterList> wrm_plist_ = Teuchos::sublist(
		  Teuchos::sublist(ats_elm_pks_plist_, "flow"), "water retention evaluator");

  std::string wrm_type = wrm_plist_->get<std::string>("WRM Type");
  auto wrm_constants_plist_ = Teuchos::sublist(wrm_plist_, "WRM parameters");
  c=0;
  for (auto& entry : *wrm_constants_plist_) {
    std::string layer_name = entry.first;

    auto layer3_plist_ = Teuchos::sublist(wrm_constants_plist_, layer_name);

    if(wrm_type == "van Genuchten"){
       double vG_alpha = layer3_plist_->get<double>("van Genuchten alpha [Pa^-1]");
       double vG_n = layer3_plist_->get<double>("van Genuchten n [-]");
       double vG_sr = layer3_plist_->get<double>("residual saturation [-]");
       double vG_s0 = layer3_plist_->get<double>("smoothing interval width [saturation]", 0.0);

       //std::cout<<"WRM -"<<layer_name<<" -vG_alpha "<<vG_alpha
       //		   <<" -vG_n "<<vG_n<<" -Sr "<<vG_sr<< " -S0 "<<vG_s0
       //		   <<std::endl;
       //(TODO) reset Parameters, if available from ELM
    } else if(wrm_type == "Clapp Hornberger"){
       double smpsat = layer3_plist_->get<double>("Clapp Hornberger smpsat [Pa]");
       double bsw = layer3_plist_->get<double>("Clapp Hornberger bsw [-]");
       double sr = layer3_plist_->get<double>("residual saturation [-]", 0.0);
       double s0 = layer3_plist_->get<double>("near-saturation inflection point interval [saturation]", 0.08);
       //double pcx = layer3_plist_->get<double>("dry-end smoothing starting point [Pa]", 1.0e10);

       smpsat = CH_smpsat[c];
       bsw = CH_bsw[c];
       sr = CH_sr[c];

       layer3_plist_->set("Clapp Hornberger smpsat [Pa]", smpsat);
       layer3_plist_->set("Clapp Hornberger bsw [-]", bsw);
       layer3_plist_->set("residual saturation [-]", sr);

    }
    c++;
    //std::cout<<"WRM - "<< layer_name <<": "<< *layer3_plist_<<std::endl;
  }
}

void ats_elm_drv::plist_cycle_driver_reset(const double t0) {
  bool success;
  Amanzi::Utils::Units units;

  // t0_
  ats_elm_drv_plist_->set("start time",t0,"s");
  ats_elm_drv_plist_->set("start time units","s");
  t0_ = ats_elm_drv_plist_->get<double>("start time");

  // t1_
  ats_elm_drv_plist_->set("end time",t0+1800.0,"s");
  ats_elm_drv_plist_->set("start time units","s");
  t1_ = ats_elm_drv_plist_->get<double>("end time");

  // dt_
  ats_elm_drv_plist_->set("max time step size [s]", 1800.00);
  max_dt_ = ats_elm_drv_plist_->get<double>("max time step size [s]", 1800.00);
  ats_elm_drv_plist_->set("min time step size [s]", 1.00);
  min_dt_ = ats_elm_drv_plist_->get<double>("min time step size [s]", 1.0e-12);
  cycle0_ = ats_elm_drv_plist_->get<int>("start cycle",0);
  cycle1_ = ats_elm_drv_plist_->get<int>("end cycle",-1);
  duration_ = ats_elm_drv_plist_->get<double>("wallclock duration [hrs]", -1.0);
  subcycled_ts_ = ats_elm_drv_plist_->get<bool>("subcycled timestep", false);

  // restart control
  restart_ = ats_elm_drv_plist_->isParameter("restart from checkpoint file");
  if (restart_) restart_filename_ = ats_elm_drv_plist_->get<std::string>("restart from checkpoint file");
}

int ats_elm_drv::drv_setup(const double start_ts, const bool initialoutput) {

  // -----------------------------------------------------------------------------------------
  // (1) create meshes, regions, state from input .xml, AND over-ride them from elm if any.

  Teuchos::ParameterList& plist = *parameter_list_;

  // (1a) geometric model and regions
  Teuchos::ParameterList reg_params = plist.sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm_) );

  // (1b) state
  Teuchos::ParameterList state_plist = plist.sublist("state");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));

  // (1c) create and register meshes
  ATS::Mesh::createMeshes(plist, comm_, gm, *S_);

  // reset mesh coordinates, if directly-readable 'mesh' file not available
  //TODO: to be a mapping of vertices from ATS to ELM, currently only do some checking
  //Amanzi::Key mesh_name = "domain";
  //mesh_vertices_reset(mesh_name, false);
  //mesh_name = "surface";
  //mesh_vertices_reset(mesh_name, false);

  // create the top level PKs
  Teuchos::ParameterList pk_tree_plist = ats_elm_drv_plist_->sublist("PK tree");
  if (pk_tree_plist.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_plist.begin();
  const std::string &pk_name = pk_tree_plist.name(pk_item);

  // (2) model setup

  // create the solution
  soln_ = Teuchos::rcp(new Amanzi::TreeVector());

  // create the pk
  Amanzi::PKFactory pk_factory;
  pk_ = pk_factory.CreatePK(pk_name, pk_tree_plist, parameter_list_, S_, soln_);

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
  dt_restart_ = -1;
  {
    //Teuchos::TimeMonitor monitor(*setup_timer_);
    StatePKsetup();
    dt_restart_ = initialize();

  }

  // appears ELM's state variable not yet ready during its initializing
  //ic_reset();

  // if output initial states
  //if (initialoutput) {
  //  double dt = dt_restart_ > 0 ? dt_restart_ : get_dt(false);

    // visualization at IC
  //  visualize(false);
  //  checkpoint(dt);
  //}

  // have to sync starting time with ELM
  t0_ = start_ts;
  ats_elm_drv_plist_->set("start time",t0_,"s");

  return 0;
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
  auto vis_plist = Teuchos::sublist(parameter_list_,"visualization");
  for (auto& entry : *vis_plist) {
    std::string domain_name = entry.first;

    if (S_->HasMesh(domain_name)) {
      // visualize standard domain
      auto mesh_p = S_->GetMesh(domain_name);
      auto sublist_p = Teuchos::sublist(vis_plist, domain_name);
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
      auto sublist_p = Teuchos::sublist(vis_plist, domain_name);

      if (sublist_p->get("visualize individually", false)) {
        // visualize each subdomain
        for (const auto& subdomain : *dset) {
          Teuchos::ParameterList sublist = vis_plist->sublist(subdomain);
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
  if (ats_elm_drv_plist_->isSublist("required times")) {
    Teuchos::ParameterList& sublist = ats_elm_drv_plist_->sublist("required times");
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

// -----------------------------------------------------------------------------
// ONE single ELM-timestep
// -----------------------------------------------------------------------------

void ats_elm_drv::cycle_driver(const double start_ts, const double end_ts, const bool resetIC_from_elm) {
  // wallclock duration -- in seconds
  const double duration(duration_);

  // staring/ending time (seconds) from ELM, usually one ELM-timestep
  double ts = end_ts - start_ts;
  t0_ = start_ts;
  t1_ = end_ts;
  max_dt_ = ts;
  ats_elm_drv_plist_->set("end time",t1_,"s");

  // -- re-register the final time
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  tsm_->RegisterTimeEvent(t1_);

  //
  if (resetIC_from_elm) {

	  // reset initial states
	  // TODO

  }

  // some checking
  std::cout<<std::endl<<"prior-ATS timestep ---------------------------------------"<<std::endl;
  std::cout<<"S_: "<<std::endl;
  get_data(S_);
  std::cout<<"S_next_: "<<std::endl;
  get_data(S_next_);


  // get the intial timestep
  double dt = get_dt(false);

  // iterate process kernels
  {

    bool fail = false;

    std::cout <<std::endl<< "INFO: ATS runs start at " << S_->time() << std::endl;

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

      std::cout << "INFO: ATS running ... " << S_->time() << " @ dt = "<< dt << std::endl;

      fail = advance(S_->time(), S_->time() + dt, dt);


    } // while not finished

    if (not fail) {std::cout << "INFO: ATS runs end at " << S_->time() << std::endl;}

  }

  // force to output for each time-step
  visualize(true);
  checkpoint(dt, true);

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

// -------------------------------------------------------------------------------------------
// the following currently only does checking mesh vertices
void ats_elm_drv::mesh_vertices_reset(std::string mesh_name, bool reset_from_elm){

  if (S_->HasMesh(mesh_name)) {
        auto mesh_ = S_->GetMesh(mesh_name);
        mesh_ ->build_columns();

        int dim = 3;
        if (mesh_name=="surface") {dim=2;}
        Amanzi::AmanziGeometry::Point coords(dim);

        // number of vertices
        int nV = mesh_ -> num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);  // or, 'ALL'?

        std::cout<< "checking vertice:  " << mesh_name <<" - Total Vertice no.: "<<nV<< std::endl;
        std::cout<< "iV - n_node_abv - Coordinates (X,Y[,Z]) "<<std::endl;
        for (int iV=0; iV<nV; iV++) {
          // get the coords of the node
          mesh_ -> node_get_coordinates(iV,&coords);

          // need to known Z indices for current node
          int n_nabvid = 0;
          int nabvid = mesh_->node_get_node_above(iV);
          while (nabvid != -1) {
        	  n_nabvid = n_nabvid + 1;
        	  nabvid = mesh_->node_get_node_above(nabvid);
          }

          // override Z coords (TODO: X, Y), if required
          if (reset_from_elm) {
        	  coords[dim-1] = elm_col_nodes[n_nabvid];  // NOT [length_nodes-nextid-1] ??? (further checking)
        	  //mesh_ -> node_set_coordinates(iV, coords);
        	  // S_->GetMesh("domain") NON-changeable?? - (TODO) needs a new thought here
          }

          std::cout<< iV <<" - "<<n_nabvid<<" - ("<< coords<<")"<< std::endl;

        }
   }

}

// reset ats initial conditions (IC)
void ats_elm_drv::ic_reset() {
  // Three (3) types of ICs
  // (1) constants, in State
  // (2)

  // (3) real IC, i.e. primary variable in PKs
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  Teuchos::RCP<Teuchos::ParameterList> surfflow_plist_ = Teuchos::sublist(pk_plist_, "overland flow");

  std::string pk_name_ = "flow";
  std::string pv_key = flow_plist_->get<std::string>("primary variable key");
  if (S_->HasField(pv_key)){
	Teuchos::RCP<Amanzi::Field> field = S_->GetField(pv_key, pk_name_);

    auto mesh_ = S_->GetMesh("domain");
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    auto &pc = *(S_->GetFieldData(pv_key)->ViewComponent("cell"));
    for (int c = 0; c < ncells; ++c) {
      pc[0][c] = soilp[c];
    }

    //
    Teuchos::RCP<Amanzi::FieldEvaluator> fm = S_->GetFieldEvaluator(pv_key);
    Teuchos::RCP<Amanzi::PrimaryVariableFieldEvaluator> solution_evaluator =
      Teuchos::rcp_dynamic_cast<Amanzi::PrimaryVariableFieldEvaluator>(fm);
    if (solution_evaluator != Teuchos::null){
      solution_evaluator->SetFieldAsChanged(S_.ptr());
      std::cout<<" Changed? "<< pv_key <<std::endl;
    }

  }

  //
  std::string pk2_name_ = "overland flow";
  std::string pv2_key = surfflow_plist_->get<std::string>("primary variable key");
  if (S_->HasField(pv2_key)){
	Teuchos::RCP<Amanzi::Field> field = S_->GetField(pv2_key, pk2_name_);

    auto mesh2_ = S_->GetMesh("surface");
    int ncells2 = mesh2_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    auto &psurf = *(S_->GetFieldData(pv2_key)->ViewComponent("cell"));
    for (int c = 0; c < ncells2; ++c) {
      psurf[0][c] = surfp[c];
    }

    //
    Teuchos::RCP<Amanzi::FieldEvaluator> fm = S_->GetFieldEvaluator(pv2_key);
    Teuchos::RCP<Amanzi::PrimaryVariableFieldEvaluator> solution_evaluator =
      Teuchos::rcp_dynamic_cast<Amanzi::PrimaryVariableFieldEvaluator>(fm);
    if (solution_evaluator != Teuchos::null){
      solution_evaluator->SetFieldAsChanged(S_.ptr());
      std::cout<<" Changed? "<< pv2_key <<std::endl;
    }

  }

  //
  //StateField_refresh(*S_);

  *S_next_ = *S_;
  if (subcycled_ts_) {
    *S_inter_ = *S_;
  } else {
    S_inter_ = S_;
  }
  pk_->CommitStep(S_->time(), S_->time(), S_);
  pk_->set_states(S_, S_inter_, S_next_);

  // visualization at IC (after re-setting)
  double dt = dt_restart_ > 0 ? dt_restart_ : get_dt(false);
  visualize(false);
  checkpoint(dt);

}

// reset ats boundray conditions (BC)
void ats_elm_drv::bc_reset() {

  // TODO
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  std::string pk_name_ = "flow";

}

// reset ats source-sink terms (SS)
void ats_elm_drv::ss_reset() {

  //PKs
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  Teuchos::RCP<Teuchos::ParameterList> surfflow_plist_ = Teuchos::sublist(pk_plist_, "overland flow");

  // overland flow SS, e.g. rain/snow-melting/rain-throughfall etc.
  // NOTE: here all treated as potential
  std::string pk_name = "overland flow";
  if((surfflow_plist_->get<bool>("source term"))){
    std::string ss_key = surfflow_plist_->get<std::string>("source key");
    if (S_next_->HasField(ss_key)){
		auto mesh_ = S_next_->GetMesh("surface");
		int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
		auto &ss = *(S_next_->GetFieldData(ss_key)->ViewComponent("cell"));
  		for (int c = 0; c < ncells; ++c) {
  		  const auto& xyc = mesh_->cell_centroid(c);
          ss[0][c] = soil_infl[c];
  	    }

        //
   	    Teuchos::RCP<Amanzi::FieldEvaluator> fm = S_next_->GetFieldEvaluator(ss_key);
   	    Teuchos::RCP<Amanzi::PrimaryVariableFieldEvaluator> solution_evaluator =
   	      Teuchos::rcp_dynamic_cast<Amanzi::PrimaryVariableFieldEvaluator>(fm);
   	    if (solution_evaluator != Teuchos::null){
   	      solution_evaluator->SetFieldAsChanged(S_next_.ptr());
   	      std::cout<<" Changed? "<< ss_key <<std::endl;
  	    }
   	    S_next_->GetFieldEvaluator(ss_key)->HasFieldChanged(S_next_.ptr(), ss_key);

	  }
  }

  //flow SS, e.g. root water extraction (i.e. transpiration) + top-soil evaporation
  pk_name = "flow";
  if((flow_plist_->get<bool>("source term"))){
    std::string ss2_key = flow_plist_->get<std::string>("source key");

    if (S_next_->HasField(ss2_key)){
  		auto mesh2_ = S_next_->GetMesh("domain");
  		int ncells = mesh2_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  		std::string ss1_name = "soil_potential_evaporation";
  		auto &ss1 = *(S_next_->GetFieldData(ss1_name)->ViewComponent("cell"));
  		std::string ss2_name = "potential_transpiration";
  		auto &ss2 = *(S_next_->GetFieldData(ss2_name)->ViewComponent("cell"));
   		for (int c = 0; c < ncells; ++c) {
  		  const auto& xyzc = mesh2_->cell_centroid(c);
  	      ss2[0][c] = root_waterextract[c];
  	      if (c==0) {
  	    	  ss1[0][c] = soil_evap[0];
  	      }
  	    }

        //
   	    Teuchos::RCP<Amanzi::FieldEvaluator> fm = S_next_->GetFieldEvaluator(ss1_name);
   	    Teuchos::RCP<Amanzi::PrimaryVariableFieldEvaluator> solution_evaluator =
   	      Teuchos::rcp_dynamic_cast<Amanzi::PrimaryVariableFieldEvaluator>(fm);
   	    if (solution_evaluator != Teuchos::null){
   	      solution_evaluator->SetFieldAsChanged(S_next_.ptr());
   	      std::cout<<" Changed? "<< ss1_name <<std::endl;
   	    }
   	    S_next_->GetFieldEvaluator(ss1_name)->HasFieldChanged(S_next_.ptr(), ss1_name);

    }

  }

  //
  *S_ = *S_next_;
  *S_inter_ = *S_;
  //pk_->CommitStep(S_->time(), S_->time(), S_);
  //pk_->set_states(S_, S_inter_, S_next_);

}

// read data
void ats_elm_drv::get_data(Teuchos::RCP<Amanzi::State> SS_) {

  // read primary variable in PKs
  Teuchos::RCP<Teuchos::ParameterList> pk_plist_ = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::RCP<Teuchos::ParameterList> flow_plist_ = Teuchos::sublist(pk_plist_, "flow");
  std::string pk_name = "flow";
  std::string domain_name = "domain";

  auto mesh_ = SS_->GetMesh(domain_name);
  int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  std::string pv_key = flow_plist_->get<std::string>("primary variable key");
  std::cout<<"flow pv_key: "<< pv_key <<std::endl;

  if (SS_->HasField(pv_key)){
	Teuchos::RCP<Amanzi::Field> field = SS_->GetField(pv_key, pk_name);
	std::cout <<"state field: "<< field->fieldname()<<" - type: "<< field->type()<<std::endl;
    auto &pv = *(SS_->GetFieldData(pv_key)->ViewComponent("cell"));
    auto &cp = *(SS_->GetFieldData("capillary_pressure_gas_liq")->ViewComponent("cell"));
    auto &ss = *(SS_->GetFieldData("water_source")->ViewComponent("cell"));
    auto &sev = *(SS_->GetFieldData("soil_potential_evaporation")->ViewComponent("cell"));
    auto &vt = *(SS_->GetFieldData("potential_transpiration")->ViewComponent("cell"));
    for (int c = 0; c < ncells; ++c) {
      const auto& xyzc = mesh_->cell_centroid(c);
      std::cout <<"coords: "<<xyzc<<"  pressure:  "<<pv[0][c]
								  <<"  capillary-pressure:  "<<cp[0][c]
								  <<"  water-source/sink:  "<<ss[0][c]
								  <<"  soil evap:  "<<sev[0][c]
								  <<"  transpiration:  "<<vt[0][c]
								  <<std::endl;
    }

  }

}

// refresh State after resetting some field
void StateField_refresh(Amanzi::State& SS){

	for (Amanzi::State::field_iterator field=SS.field_begin(); field!=SS.field_end(); ++field) {
      if (SS.HasFieldEvaluator(field->first)) {
        Teuchos::RCP<Amanzi::FieldEvaluator> fm = SS.GetFieldEvaluator(field->first);
        Teuchos::RCP<Amanzi::PrimaryVariableFieldEvaluator> solution_evaluator =
          Teuchos::rcp_dynamic_cast<Amanzi::PrimaryVariableFieldEvaluator>(fm);
        if (solution_evaluator != Teuchos::null){
          Teuchos::Ptr<Amanzi::State> S_ptr(&SS);
          solution_evaluator->SetFieldAsChanged(S_ptr);
        }
      }
      field->second->set_initialized();
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

