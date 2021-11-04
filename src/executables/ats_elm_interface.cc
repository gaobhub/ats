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
int ats_elm_input(std::string input_filename) {
    //TODO
	return 0;
}

void ats_elm_setmpicomm(const int comm){
	//TODO
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

