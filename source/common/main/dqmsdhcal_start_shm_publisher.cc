  /// \file start_sdhcal_shm_publisher.cc
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
#include "dqm4hep/DQMCoreTool.h"

// -- tclap headers
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

// -- levbdim headers
#include "buffer.hh"
#include "shmdriver.hh"

// -- lcio headers
#include "EVENT/LCIO.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"

namespace local_namespace
{

class LMGeneric: public IMPL::LCGenericObjectImpl
{
public:
	inline std::vector<int>& getIntVector() { return _intVec; }
};

}

void copyRuContent( local_namespace::LMGeneric *pLMGeneric , levbdim::buffer *pBuffer );
void publishBufferContent( levbdim::buffer *pBuffer , const std::string &directory );


int main(int argc, char* argv[])
{
	dqm4hep::DQM4HEP::screenSplash();

	std::string cmdLineFooter = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine *pCommandLine = new TCLAP::CmdLine(cmdLineFooter, ' ', DQM4HEP_VERSION_STR);
	std::string log4cxx_file = std::string(DQMCore_DIR) + "/conf/defaultLoggerConfig.xml";

	TCLAP::ValueArg<std::string> collectionNameArg(
				  "c"
				 , "collection-name"
				 , "The SDHCAL raw data collection name"
				 , false
				 , "RU_XDAQ"
				 , "string");
	pCommandLine->add(collectionNameArg);

	TCLAP::ValueArg<unsigned int> sleepTimeArg(
				  "t"
				 , "sleep-time"
				 , "The time sleep between two consecutive published events (unit ms)"
				 , false
				 , 1000
				 , "unsigned int");
	pCommandLine->add(sleepTimeArg);

	TCLAP::ValueArg<std::string> monitoringDirectoryArg(
				  "m"
				 , "monitoring-directory"
				 , "The name of directory in which to write data"
				 , false
				 , "/dev/shm/levbdim"
				 , "string");
	pCommandLine->add(monitoringDirectoryArg);

	TCLAP::ValueArg<unsigned int> skipNEventsArg(
				  "n"
				 , "skip-event"
				 , "The number of events to skip form file beginning"
				 , false
				 , 0
				 , "unsigned int");
	pCommandLine->add(skipNEventsArg);

	std::vector<std::string> allowedLevels;
	allowedLevels.push_back("INFO");
	allowedLevels.push_back("WARN");
	allowedLevels.push_back("DEBUG");
	allowedLevels.push_back("TRACE");
	allowedLevels.push_back("ERROR");
	allowedLevels.push_back("FATAL");
	allowedLevels.push_back("OFF");
	allowedLevels.push_back("ALL");
	TCLAP::ValuesConstraint<std::string> allowedLevelsContraint( allowedLevels );

	TCLAP::ValueArg<std::string> verbosityArg(
				  "v"
				 , "verbosity"
				 , "The verbosity level used for this application"
				 , false
				 , "INFO"
				 , &allowedLevelsContraint);
	pCommandLine->add(verbosityArg);

	TCLAP::ValueArg<std::string> loggerConfigArg(
				  "l"
				 , "logger-config"
				 , "The xml logger file to configure log4cxx"
				 , false
				 , log4cxx_file
				 , "string");
	pCommandLine->add(loggerConfigArg);

	TCLAP::ValueArg<std::string> lcioFileNamesArg(
				  "f"
				 , "lcio-files"
				 , "The list of lcio files to process (separated by a ':' character)"
				 , true
				 , ""
				 , "string");
	pCommandLine->add(lcioFileNamesArg);

	// parse command line
	pCommandLine->parse(argc, argv);

	log4cxx_file = loggerConfigArg.getValue();
	log4cxx::xml::DOMConfigurator::configure(log4cxx_file);

	if( verbosityArg.isSet() )
		dqm4hep::dqmMainLogger->setLevel( log4cxx::Level::toLevel( verbosityArg.getValue() ) );

	std::vector<std::string> lcioInputFiles;
	dqm4hep::DQM4HEP::tokenize(lcioFileNamesArg.getValue(), lcioInputFiles, ":");

	// lcio file reader
	IO::LCReader *pLCReader = IOIMPL::LCFactory::getInstance()->createLCReader(0);

	// open files
	pLCReader->open( lcioInputFiles );

	// skip first events if asked
	if( skipNEventsArg.getValue() )
	   pLCReader->skipNEvents( skipNEventsArg.getValue() );

	// global buffer to publish data
	levbdim::buffer buffer(0x20000);

	EVENT::LCEvent *pLCEvent = NULL;

	int ret(0);

	while( ( pLCEvent = pLCReader->readNextEvent() ) )
	{
		if( ret != 0 )
			break;

		EVENT::LCCollection *pLCCollection = pLCEvent->getCollection( collectionNameArg.getValue() );

		if( pLCCollection->getTypeName() != EVENT::LCIO::LCGENERICOBJECT )
		{
			LOG4CXX_FATAL( dqm4hep::dqmMainLogger , "Wrong collection type. Expected LCGenericObject collection !" );
			ret = 1;
			continue;
		}

		for( unsigned int e=0 ; e<pLCCollection->getNumberOfElements() ; e++ )
		{
			local_namespace::LMGeneric *pLMGeneric = (local_namespace::LMGeneric *)( pLCCollection->getElementAt(e) );

			if( ! pLMGeneric )
			{
				LOG4CXX_FATAL( dqm4hep::dqmMainLogger , "Wrong object type. Expected LCGenericObject collection !" );
				ret = 1;
				break;
			}


			copyRuContent( pLMGeneric , &buffer );
			publishBufferContent( &buffer , monitoringDirectoryArg.getValue() );
		}

		dqm4hep::DQMCoreTool::sleep( std::chrono::milliseconds( sleepTimeArg.getValue() ) );
	}

	delete pCommandLine;
	delete pLCReader;

	return ret;
}




void copyRuContent( local_namespace::LMGeneric *pLMGeneric , levbdim::buffer *pBuffer )
{
	// get raw buffer and size
	int *pRawBuffer = & pLMGeneric->getIntVector()[0];
	unsigned char *RawBufferUChar = (unsigned char *) pRawBuffer;
	uint32_t bufferSize = pLMGeneric->getNInt()*sizeof(int32_t);

	// copy ru to buffer
	memcpy( pBuffer->ptr(), RawBufferUChar, bufferSize );
	pBuffer->setPayloadSize( bufferSize - 3*sizeof(uint32_t) - sizeof(uint64_t) );
}




void publishBufferContent( levbdim::buffer *pBuffer , const std::string &directory )
{
	levbdim::shmdriver::store(
			pBuffer->detectorId(),
			pBuffer->dataSourceId(),
			pBuffer->eventId(),
			pBuffer->bxId(),
			pBuffer->ptr(),
			pBuffer->size(),
			directory );
}



