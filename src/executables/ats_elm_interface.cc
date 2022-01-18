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
int ats_elm_init(const char* c_input_file, const int comm, const double start_ts) {

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

  // setup ats and initialize fully
  try {
	int iret = ats_elm.drv_setup(0.0, true);
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

  //ats_elm.patm = patm;
  //ats_elm.wtd = wtd;
  int n = ats_elm.length_nodes-1;
  ats_elm.soilp = new double[n];
  for (int i=0; i<ats_elm.length_nodes-1; i++) {ats_elm.soilp[i] = soilpressure[i];}

}

// boundary conditions
void ats_elm_setBC(){
	//TODO
}

// source/sink terms
void ats_elm_setSS(){
	//TODO
}

// run one timestep
void ats_elm_onestep(const double start_ts, const double end_ts,
		const int resetIC_elm, const int restart){

  // TODO comm
  int rank = ats_elm.comm_rank;

  if (rank == 0) {
    std::cout << std::endl;
    std::cout << "INFO: cycle-driver activiated from ELM for time period of "
			<< start_ts << " -- " << end_ts << " second" << std::endl;
    std::cout << "reset ATS IC: " << resetIC_elm  <<std::endl;
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
	//TODO
}

