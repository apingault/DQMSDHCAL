  /// \file dqmsdhcal_dump_db_geometry.cc
/*
 *
 * dqmsdhcal_dump_db_geometry.cc main source file template automatically generated
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

// -- sdhcal dmq headers
#include "Geometry.h"

using namespace std;
using namespace dqm4hep;
using namespace dqm_sdhcal;

int main(int argc, char* argv[])
{
	std::string cmdLineFooter = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine *pCommandLine = new TCLAP::CmdLine(cmdLineFooter, ' ', DQM4HEP_VERSION_STR);

	TCLAP::ValueArg<std::string> passwordArg(
				  "p"
				 , "password"
				 , "The database password"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(passwordArg);

	TCLAP::ValueArg<std::string> userArg(
				  "u"
				 , "user"
				 , "The database user"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(userArg);

	TCLAP::ValueArg<std::string> hostArg(
				  "k"
				 , "host"
				 , "The database host"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(hostArg);

	TCLAP::ValueArg<std::string> databaseArg(
				  "d"
				 , "database"
				 , "The database to use"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(databaseArg);

	TCLAP::ValueArg<std::string> fileArg(
				  "f"
				 , "file"
				 , "The xml geometry output file"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(fileArg);

	TCLAP::ValueArg<std::string> beamTestArg(
				  "b"
				 , "beam-test"
				 , "The beam test version (string) to dump"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(beamTestArg);

	// parse command line
	pCommandLine->parse(argc, argv);

	std::string log4cxx_file = std::string(DQMCore_DIR) + "/conf/defaultLoggerConfig.xml";
	log4cxx::xml::DOMConfigurator::configure(log4cxx_file);

	try
	{
		GeometryDBInterface interface;

		THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, interface.connect(
				hostArg.getValue() ,
				userArg.getValue() ,
				passwordArg.getValue() ,
				databaseArg.getValue() ) );

		THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, interface.dumpGeometry( fileArg.getValue(), beamTestArg.getValue() ) );
	}
	catch(dqm4hep::StatusCodeException &exception)
	{
		std::cout << "Couldn't dump geometry to xml file : " << exception.toString() << std::endl;
		return exception.getStatusCode();
	}

	std::cout << "OK" << std::endl;

	return 0;
}
