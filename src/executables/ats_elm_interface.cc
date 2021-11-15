/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS interface for ELM

License: see $ATS_DIR/COPYRIGHT
Authors: Ethan Coon, F.-M. Yuan @ ORNL

Implementation for the ats_elm_interface.  The interface is a class to:
(1) read-in inputs by filename passing from ELM, and communicator synchorizing;
(2) pass data from ELM to ATS;
(3) in ONE single ELM-timestep, run cycle driver, which runs the overall, top level timestep loop.
As in 'coordinator.cc', this portion of codes instantiates states, ensures they are initialized,
and runs the timestep loop, including Vis and restart/checkpoint dumps. It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.
(4) return data to ELM.
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

  std::cout << "INFO: inputfile from ELM " << std:: endl
		    << c_input_file << std::endl;

  // TODO comm
  int rank = 0;

  //initializing ats_elm_driver
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

// initial conditions
void ats_elm_setIC(){
	//TODO
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
void ats_elm_onestep(){
	//TODO
}

// return data to ELM's ATS interface
void ats_elm_getdata(){
	//TODO
}

