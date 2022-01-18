/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* extern "C" interface of ATS-ELM coupling

   (1) Its fortran-code interface is named as 'ELM_ATS_InterfaceMod.F90'.

   (2) The building system will compile and link into a package of dynamic libraries,
       with libats_elm_interface.dylib as hooking portal, including this interface and 'ats_elm_drv' class.

   (3) ALL inputs, inc. *.xml, required by ATS would be in E3SM's input data directory under 'lnd/clm2/ats'

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
 ELM-EMI/em-ats calling subroutines include (NOT YET finalized):
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


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

ATS::ats_elm_drv ats_elm;

// input .xml file passing from ELM as c_char array
int ats_elm_init(const char* input_filename, const int comm, const double start_time);  // input .xml file passing from ELM

void ats_elm_setmesh(const double* surf_gridsX, const double* surf_gridsY,
		const double* surf_gridsZ, const double *col_nodes,
		const int len_gridsX, const int len_gridsY, const int len_nodes);
void ats_elm_setIC(const double* patm,
		const double* soilpressure,
		const double* wtd);      // initial conditions
void ats_elm_setBC();      // boundary conditions
void ats_elm_setSS();      // source/sink terms
void ats_elm_onestep(const double start_time, const double end_time, // advance one ELM-timestep
		const int resetIC_from_elm=0,                                // 0 or 1 pass from ELM
        const int restart_from_ats=0);                               // 0 or 1 pass from ELM, must have 'checkingpoint_final.h5' properly
void ats_elm_getdata();    // extract data and return to ELM


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
