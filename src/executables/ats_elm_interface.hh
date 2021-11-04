/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* Simulation controller and top-level driver For ONE E3SM Land Model (ELM) Timestep

 The interface imitates ats 'main.cc', 'simulation_driver.cc', 'coordinator.cc' very much, WITHOUT heavy print-out info
 So the spec. are exactly same.

 The building system will compile and link into a package of dynamic libraries,
 with libats_elm_interface.dylib as hooking class to ELM's EMI-EM_ATS interface.

 ALL inputs, inc. *.xml, required by ATS would be in E3SM's input data directory under 'lnd/clm2/ats'

 Authors: F.-M. Yuan (yuanf@ornl.gov), Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

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

#ifndef ATS_ELM_INTERFACE_HH_
#define ATS_ELM_INTERFACE_HH_

#include "ats_elm_drv.hh"

ATS::ats_elm_drv ats_elm();

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

//int ats_elm_init(std::string input_filename,  // input .xml file passing from ELM
//		  const Teuchos::RCP<const Amanzi::Comm_type>& comm); // communicator passing from ELM
int ats_elm_input(std::string input_filenam); // input .xml file passing from ELM
void ats_elm_setmpicomm(const int comm);   //TODO - comm is a fake type now
void ats_elm_setIC();      // initial conditions
void ats_elm_setBC();      // boundary conditions
void ats_elm_setSS();      // source/sink terms
void ats_elm_onestep();    // advance one ELM-timestep
void ats_elm_getdata();    // extract data and return to ELM


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
