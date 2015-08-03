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

//// -- streamlog headers
//#include "streamlog/streamlog.h"

// -- dqm4hep headers
#include "dqm4hep/core/DQM4HEP.h"
#include "dqm4hep/core/DQMLogging.h"

// -- tclap headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

// -- dim headers
#include "dis.hxx"

// -- sdhcal dmq headers
#include "dqmsdhcal/daq/ShmReader.h"

using namespace std;
using namespace dqm4hep;
using namespace dqm_sdhcal;

int main(int argc, char* argv[])
{
	DQM4HEP::screenSplash();

	std::string cmdLineFooter = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine *pCommandLine = new TCLAP::CmdLine(cmdLineFooter, ' ', DQM4HEP_VERSION_STR);

	TCLAP::ValueArg<std::string> dataCollectorNameArg(
				  "c"
				 , "collector-name"
				 , "The data collector name"
				 , true
				 , "DEFAULT"
				 , "string");
	pCommandLine->add(dataCollectorNameArg);

	TCLAP::ValueArg<unsigned int> numberOfDifsArg(
				  "n"
				 , "number-of-difs"
				 , "The number of difs in the setup to consider"
				 , true
				 , 144
				 , "unsigned int");
	pCommandLine->add(numberOfDifsArg);

	TCLAP::ValueArg<std::string> detectorNameArg(
				  "d"
				 , "detector-name"
				 , "The detector name related to the setup"
				 , false
				 , "SDHCAL"
				 , "string");
	pCommandLine->add(detectorNameArg);

	TCLAP::ValueArg<std::string> collectionNameArg(
				  "e"
				 , "collection-name"
				 , "The name of the LCGenericObject collection send to the dqm system"
				 , false
				 , "RU_XDAQ"
				 , "string");
	pCommandLine->add(collectionNameArg);

	TCLAP::ValueArg<std::string> monitoringDirectoryArg(
				  "m"
				 , "monitoring-directory"
				 , "The name of directory in which to read data files"
				 , true
				 , "/dev/shm/monitoring"
				 , "string");
	pCommandLine->add(monitoringDirectoryArg);

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
	std::string logHeader = detectorNameArg.getValue() + " SHM READER";
	streamlog_init( logHeader , verbosity );

	ShmReader *pShmReader = NULL;

	try
	{
		pShmReader = new ShmReader(dataCollectorNameArg.getValue());

		pShmReader->setNumberOfDifs(numberOfDifsArg.getValue());
		pShmReader->setDetectorName(detectorNameArg.getValue());
		pShmReader->setCollectionName(collectionNameArg.getValue());
		pShmReader->setMonitoringDirectory(monitoringDirectoryArg.getValue());
	}
	catch(dqm4hep::StatusCodeException &exception)
	{
		std::cout << "Couldn't instantiate shm reader : " << exception.getStatusCode() << std::endl;
		return 1;
	}

	try
	{
		pShmReader->start();
	}
	catch(dqm4hep::StatusCodeException &exception)
	{
		std::cout << "Exception caught while running shm reader : " << exception.getStatusCode() << std::endl;
		return 1;
	}

	delete pCommandLine;
	delete pShmReader;

	return 0;
}
