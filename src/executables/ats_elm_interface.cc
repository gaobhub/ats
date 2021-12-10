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

  // setup ats and initialize fully
  try {
	iret = ats_elm.drv_setup(start_ts, true);
  }catch (int& ierr) {
	if (rank == 0) {
    std::cerr << "ERROR in ats_elm driver setup and initialization with error code: " << ierr << std::endl;
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
void ats_elm_onestep(const double start_ts, const double end_ts, const int resetIC){
  //TODO
  // TODO comm
  int rank = ats_elm.comm_rank;

  if (rank == 0) {
    std::cout << std::endl;
    std::cout << "INFO: cycle-driver activiated from ELM for time period of "
			<< start_ts << " -- " << end_ts << " second" << std::endl;
  }

  //
  try {

    bool reset = true;
    if (resetIC == 0) {
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

