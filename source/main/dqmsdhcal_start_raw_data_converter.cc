  /// \file start_sdhcal_shm_reader.cc
/*
 *
 * start_sdhcal_shm_reader.cc main source file template automatically generated
 * Creation date : mer. nov. 5 2014
 *
 * This file is part of DQMSDHCAL libraries.
 *
 * DQMSDHCAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * DQMSDHCAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DQMSDHCAL.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @author Remi Ete
 * @copyright CNRS , IPNL
 */

#include <iostream>

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMLogging.h"

// -- tclap headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

// -- dim headers
#include "dis.hxx"

// -- sdhcal dmq headers
#include "dqmsdhcal/services/DQMRawDataConverterApplication.h"

using namespace std;
using namespace dqm4hep;
using namespace dqm_sdhcal;

DQMRawDataConverterApplication *pApplication = NULL;

// simple function to exit the program
void exit_application(int returnCode)
{
	streamlog_out(MESSAGE) << "Exiting application !" << std::endl;

	if(NULL != pApplication)
	{
		pApplication->exit( returnCode );
		delete pApplication;
		pApplication = NULL;
	}

	exit(returnCode);
}

// key interrupt signal handling
void int_key_signal_handler(int signal)
{
	if(NULL == pApplication)
		exit(0);

	streamlog_out(WARNING) << "*** SIGN INT ***" << std::endl;

	streamlog_out(MESSAGE) << "Caught signal " << signal << ". Stopping the application." << std::endl;
	exit_application( static_cast<int>(STATUS_CODE_SUCCESS) );
}

// segmentation violation signal handling
void seg_viol_signal_handling(int signal)
{
	if(NULL == pApplication)
		exit(1);

	streamlog_out(MESSAGE) << "*** SEG VIOL ***" << std::endl;
	streamlog_out(MESSAGE) << "Caught signal : " << signal << std::endl;
	exit_application( static_cast<int>(STATUS_CODE_INVALID_PTR) );
}

int main(int argc, char* argv[])
{
	DQM4HEP::screenSplash();

	std::string cmdLineFooter = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine *pCommandLine = new TCLAP::CmdLine(cmdLineFooter, ' ', DQM4HEP_VERSION_STR);

	TCLAP::ValueArg<std::string> steeringFileNameArg(
				  "f"
				 , "steering-file"
				 , "The steering file used for this application (JSON format)"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(steeringFileNameArg);

	TCLAP::ValueArg<std::string> verbosityArg(
				  "v"
				 , "verbosity"
				 , "The verbosity used for this application"
				 , false
				 , "WARNING"
				 , "string");
	pCommandLine->add(verbosityArg);

	// parse command line
	pCommandLine->parse(argc, argv);

	std::string verbosity = verbosityArg.getValue();
	streamlog_init( "DQM4HEP MODULE APPLICATION" , verbosity );

	// install signal handlers
	streamlog_out(MESSAGE) << "Installing signal handlers ... " << std::endl;
	signal(SIGINT,  int_key_signal_handler);
	signal(SIGSEGV, seg_viol_signal_handling);

	streamlog_out(MESSAGE) << "Creating raw data converter application ... " << std::endl;

	try
	{
		pApplication = new DQMRawDataConverterApplication();
	}
	catch(StatusCodeException &exception)
	{
		streamlog_out(ERROR) << "StatusCodeException caught : " << exception.toString() << std::endl;
		exit_application( exception.getStatusCode() );
	}

	streamlog_out(MESSAGE) << "Creating raw data converter application ... OK" << std::endl;

	streamlog_out(MESSAGE) << "Raw data converter application read settings ..." << std::endl;

	try
	{
		THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, pApplication->readSettings(steeringFileNameArg.getValue()));
	}
	catch(StatusCodeException &exception)
	{
		streamlog_out(ERROR) << "StatusCodeException caught : " << exception.toString() << std::endl;
		exit_application( exception.getStatusCode() );
	}

	streamlog_out(MESSAGE) << "Raw data converter application read settings ... OK" << std::endl;


	streamlog_out(MESSAGE) << "Running raw data converter application ... " << std::endl;
	try
	{
		THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, pApplication->run());
	}
	catch(StatusCodeException &exception)
	{
		streamlog_out(ERROR) << "StatusCodeException caught : " << exception.toString() << std::endl;
		exit_application( exception.getStatusCode() );
	}

	delete pCommandLine;
	exit_application( static_cast<int>(STATUS_CODE_SUCCESS) );
}
