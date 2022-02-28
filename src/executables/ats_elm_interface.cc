/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS-ELM interface, direct portal to ELM-ATS fortran interface.

License: see $ATS_DIR/COPYRIGHT
Authors: Ethan Coon, F.-M. Yuan @ ORNL

The interface is a series of extern "C" functions of c-interface,
which are directly called by ELM-ATS fortran interface in ELM-emi/em-ats subroutines.

NOTE: this interface appears not be in a class, otherwise hard to be called by fortran interface.
      (not sure why currently)
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

// -----------------------------------------------------------------------------
//
// functions called by ELM-emi/ATS interface (in *.F90)
//
// -----------------------------------------------------------------------------
#include "ats_elm_interface.hh"

// input file read from ELM calling
int ats_elm_init(const char* c_input_file, const int comm) {

  ats_elm = ATS::ats_elm_drv();

  // TODO comm
  int rank = ats_elm.comm_rank;

  if (rank == 0) {
    std::cout << "INFO: inputfile from ELM " << std:: endl
		    << c_input_file << std::endl;
  }

  //initializing ats_elm_driver (1) input file and communicator (TODO)
  int iret = 0;
  try {
	std::string input_filename(c_input_file);  // converting from c_char* to std:string
	iret = ats_elm.drv_init(input_filename);
  } catch (std::string& s) {
	if (rank == 0) {
      std::cerr << "ERROR:" << std::endl
                << s << std::endl;
    }
    return 1;
  } catch (int& ierr) {
	if (rank == 0) {
      std::cerr << "ERROR: unknown error code " << ierr << std::endl;
	}
	return ierr;
  }

  return iret;
}

//
void ats_elm_setmesh(const double* surf_gridsX, const double* surf_gridsY,
		const double* surf_gridsZ, const double *col_nodes,
		const int len_gridsX, const int len_gridsY, const int len_nodes) {

  // need to check the sizes of pointer data
  ats_elm.length_gridsX = len_gridsX;
  ats_elm.length_gridsY = len_gridsY;
  ats_elm.length_nodes = len_nodes;

  // coordinates pass from ELM
  ats_elm.elm_surf_gridsX = new double[len_gridsX];              // elm surface-grid X coord in m
  for (int i=0; i<len_gridsX; i++) {ats_elm.elm_surf_gridsX[i] = surf_gridsX[i];}
  ats_elm.elm_surf_gridsY = new double[len_gridsY];              // elm surface-grid Y coord in m
  for (int i=0; i<len_gridsY; i++) {ats_elm.elm_surf_gridsY[i] = surf_gridsY[i];}
  ats_elm.elm_surf_gridsZ = new double[len_gridsX*len_gridsY];;  // elm surface-grid elevation in m
  // (TODO) surface node elevation passing

  ats_elm.elm_col_nodes = new double[len_nodes];  // elm soil column nodes in m (elevation)
  for (int i=0; i<len_nodes; i++) {ats_elm.elm_col_nodes[i] = col_nodes[i];}

}

//
void ats_elm_setmat(const double* porosity, const double* hksat,
		const double* CH_bsw, const double* CH_smpsat, const double* CH_sr,
		const double* eff_porosity){

  int n = ats_elm.length_nodes-1;
  ats_elm.porosity = new double[n];
  ats_elm.hksat    = new double[n];
  ats_elm.CH_bsw   = new double[n];
  ats_elm.CH_smpsat= new double[n];
  ats_elm.CH_sr    = new double[n];
  ats_elm.eff_porosity = new double[n];
  for (int i=0; i<ats_elm.length_nodes-1; i++) {
    ats_elm.porosity[i] = porosity[i];
    ats_elm.hksat[i]    = hksat[i];
    ats_elm.CH_bsw[i]   = CH_bsw[i];
    ats_elm.CH_smpsat[i]= CH_smpsat[i];
    ats_elm.CH_sr[i]    = CH_sr[i];
    ats_elm.eff_porosity[i] = eff_porosity[i];

  }

  std::cout<<"INFO:: material data READY! "<<std::endl;

}

// setup ats and initialize fully
void ats_elm_setup(const double start_ts){
  try {

    ats_elm.plist_reset(start_ts);

    int iret = ats_elm.drv_setup(start_ts, true);

  }catch (int& ierr) {
	if (ats_elm.comm_rank == 0) {
      std::cerr << "ERROR in ats_elm driver setup and initialization with error code: " << ierr << std::endl;
	}
  }

}

// initial conditions
void ats_elm_setIC(const double* patm,
		const double* soilpressure,
		const double* wtd){
  //

  int g = (ats_elm.length_gridsX-1)*(ats_elm.length_gridsY-1); //TODO - need to check how c++ 2D-array arranged
  ats_elm.patm = new double[g];
  ats_elm.wtd = new double[g];
  ats_elm.surfp = new double[g];
  for (int i=0; i<g; i++) {
	    ats_elm.patm[i] = patm[i];
	    ats_elm.wtd[i] = wtd[i];
	    ats_elm.surfp[i] = soilpressure[0]; // temporarilly set (TODO: better get surface water depth from ELM)
  }

  int n = ats_elm.length_nodes-1;
  ats_elm.soilp = new double[n];
  for (int j=0; j<ats_elm.length_nodes-1; j++) {ats_elm.soilp[j] = soilpressure[j];}

}

// boundary conditions
void ats_elm_setBC(){
	//TODO
}

// source/sink terms
void ats_elm_setSS(const double* ss_soilinfl, const double* ss_soilevap, const double* ss_soilbottom,
		const double* ss_roottran, const double* ss_other){

  int c = (ats_elm.length_gridsX-1)*(ats_elm.length_gridsY-1); //TODO - need to check how c++ 2D-array arranged
  ats_elm.soil_infl = new double[c];
  ats_elm.soil_evap = new double[c];
  for (int i=0; i<c; i++) {
    ats_elm.soil_infl[i] = ss_soilinfl[i];
    ats_elm.soil_evap[i] = ss_soilevap[i];
    if (ats_elm.soil_infl[i]>0.0) {
    	std::cout<<"soil gross infiltration from ELM: "<<i<<" - "<<ats_elm.soil_infl[i]<<std::endl;
    }
    if (ats_elm.soil_evap[i]!=0.0) {
        std::cout<<"soil potential evaporation from ELM: "<<i<<" - "<<ats_elm.soil_evap[i]<<std::endl;
    }
  }

  int n = ats_elm.length_nodes-1;
  ats_elm.root_waterextract = new double[n];  //TODO - need to check how c++ 2D-array arranged
  for (int j=0; j<n; j++) {
	  ats_elm.root_waterextract[j] = ss_roottran[j];
	  if (ats_elm.root_waterextract[j]>0.0) {
	     std::cout<<"root-transpiration from ELM: "<<j<<" - "<<ats_elm.root_waterextract[j]<<std::endl;
	  }
  }
}

// run one timestep
void ats_elm_onestep(const double start_ts, const double end_ts,
		const int resetIC_elm, const int restart){

  int rank = ats_elm.comm_rank;

  if (rank == 0) {
    std::cout << std::endl;
    std::cout << "INFO: cycle-driver activiated from ELM for time period of "
			<< start_ts << " -- " << end_ts << " second" << std::endl;
  }

  //over-ride ATS initial conditions, if need
  if (resetIC_elm==1){
      ats_elm.ic_reset();
  }

  // always put into ATS BCs/SS, each ELM time-step
  ats_elm.bc_reset();  // NOT YET
  ats_elm.ss_reset();

  //
  try {

    bool reset = true;
    if (restart == 0) {
    	reset = false;
    }

    ats_elm.cycle_driver(start_ts, end_ts, reset);
	if (rank == 0) {
      std::cout << "INFO: ATS runs sucessfully! " << std::endl;
      std::cout << std::endl;
    }

  } catch (std::string& s) {
	if (rank == 0) {
      std::cerr << "ERROR:" << std::endl
                << s << std::endl;
    }
  } catch (int& ierr) {
	if (rank == 0) {
      std::cerr << "ERROR: unknown error code " << ierr << std::endl;
	}
  }
}

// return data to ELM's ATS interface
void ats_elm_getdata(){
	int rank = ats_elm.comm_rank;

	//
	try {

	    ats_elm.get_data(ats_elm.S_);
		if (rank == 0) {
	      std::cout << "---------------------------------------------------------"
	    		  "------ "<<std::endl<<std::endl;
		}
	} catch (std::string& s) {
		if (rank == 0) {
	      std::cerr << "ERROR:" << std::endl
	                << s << std::endl;
	    }
	} catch (int& ierr) {
		if (rank == 0) {
	      std::cerr << "ERROR: unknown error code " << ierr << std::endl;
		}
	}


	//TODO: data-passing to ELM
}

