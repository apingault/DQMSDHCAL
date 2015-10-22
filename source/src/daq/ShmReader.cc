  /// \file ShmReader.cc
/*
 *
 * ShmReader.cc source template automatically generated by a class generator
 * Creation date : mar. mai 5 2015
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

// -- sdhcal dqm headers
#include "dqmsdhcal/daq/ShmReader.h"
#include "dqmsdhcal/daq/Decoder.h"
#include "dqmsdhcal/daq/DIFReadoutConstant.h"

// -- std headers
#include <dirent.h>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// -- lcio headers
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCGenericObjectImpl.h"

// -- dqm4hep headers
#include "dqm4hep/core/DQMPluginManager.h"
#include "dqm4hep/network/DQMDataSender.h"
#include "dqm4hep/network/DQMNetworkTool.h"
#include "dqm4hep/lcio/DQMLCEventStreamer.h"
#include "dqm4hep/lcio/DQMLCEvent.h"

namespace dqm_sdhcal
{

int file_select(const struct dirent *entry)
{
	if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))
		return (0);
	else
		return (1);
}

//-------------------------------------------------------------------------------------------------

ShmReader::ShmReader(const std::string &eventCollectorName) :
		m_currentEventNumber(1),
		m_detectorName("SDHCAL"),
		m_collectionName("RU_XDAQ"),
		m_numberOfDifs(0)
{
	try
	{
		m_pEventSender = new dqm4hep::DQMDataSender();
		dqm4hep::DQMLCEventStreamer *pLCEventStreamer = NULL;
		THROW_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMPluginManager::instance()->getCastedPluginClone("LCIOStreamer", pLCEventStreamer));
		m_pEventSender->setEventStreamer(pLCEventStreamer);
		m_pEventSender->setCollectorName(eventCollectorName);

		if(!dqm4hep::DQMNetworkTool::dataCollectorExists(eventCollectorName))
			streamlog_out(WARNING) << "SDHCAL ShmReader : collector '" << eventCollectorName << "' not yet registered over the network" << std::endl;
	}
	catch(dqm4hep::StatusCodeException &statusCodeException)
	{
		if(m_pEventSender)
		{
			delete m_pEventSender;
			m_pEventSender = NULL;
		}

		throw statusCodeException;
	}
}

//-------------------------------------------------------------------------------------------------

ShmReader::~ShmReader() 
{
	if(m_pEventSender)
		delete m_pEventSender;
}

//-------------------------------------------------------------------------------------------------

void ShmReader::start()
{
	// Requires n dif > 0
	if(0 == m_numberOfDifs)
	{
		std::cout << "Invalid number of difs (0). Excepted > 0 !" << std::endl;
		throw dqm4hep::StatusCodeException(dqm4hep::STATUS_CODE_INVALID_PARAMETER);
	}

	// the map of raw sdhcal buffers for each dif (see class definition)
	BufferMap bufferMap;

	// infinite loop over incoming events
	while (1)
	{
		streamlog_out(DEBUG) << "On liste" << std::endl;

		struct dirent **files;
		unsigned char cbuf[0x20000];
		std::string monitoringClosedDir = m_monitoringDirectory + "/closed";

		// read the directory contents. Returns the number of entries found in the directory
		int count = scandir(monitoringClosedDir.c_str(), &files, file_select, alphasort);

		if(count <= 0)
		{
			usleep(100);
			continue;
		}

		bool shouldContinue = false;

		// loop over found files
		for (int i=1; i<count+1; ++i)
		{
			uint32_t dif,dtc,gtc;
			unsigned long long abcid;

			// get the previous variables from the file name
			sscanf(files[i-1]->d_name,"%lld_%d_%d_%d",&abcid,&dtc,&gtc,&dif);

			// open the event file (read only)
			std::stringstream eventFileName;
			eventFileName << m_monitoringDirectory << "/Event_" << abcid << "_" << dtc << "_" << gtc << "_" << dif;

			int fd = open(eventFileName.str().c_str(), O_RDONLY);

			// check for error
			if (fd < 0)
			{
				printf("%d rc\n",fd);
				usleep(100);
				shouldContinue = true;
				break;
			}

			// copy the content in cbuf buffer
			int size_buf = read(fd, cbuf, 0x20000);

		    close(fd);
		    unlink(eventFileName.str().c_str());

		    eventFileName.str("");
		    eventFileName << m_monitoringDirectory << "/closed/" << abcid << "_" << dtc << "_" << gtc << "_" << dif;
		    unlink(eventFileName.str().c_str());

		    uint32_t isizeb = size_buf/4+1;
		    unsigned char* cdata = new unsigned char[isizeb*4];

		    memcpy(cdata, cbuf, size_buf);

		    uint32_t ib0 = 24;
		    gtc = Decoder::getBufferGTC(cbuf, ib0);
		    dtc = Decoder::getBufferDTC(cbuf, ib0);
		    uint32_t difid = Decoder::getBufferDIF(cbuf, ib0);
		    uint32_t idx_abcid = gtc;

		    streamlog_out(DEBUG) << "FIFO read " << difid << " : "
		    		<< dtc << " "
		    		<< gtc << " "
		    		<< abcid << std::endl;

		    // find the entry in the buffer map
		    BufferMap::iterator it_gtc = bufferMap.find(idx_abcid);

		    // add buffer in the entry if found
		    if (it_gtc != bufferMap.end())
		    {
		      it_gtc->second.push_back(cdata);
		    }
		    // else create the entry
		    else
		    {
				std::vector<unsigned char*> v;
				v.push_back(cdata);
				bufferMap.insert( BufferMap::value_type(idx_abcid, v) );
				it_gtc = bufferMap.find(gtc);
			}
		}

		if(shouldContinue)
		{
			usleep(100);
			continue;
		}

		// Now loop on the Buffer Map and analyze completed events
		for (BufferMap::iterator it = bufferMap.begin() ; it != bufferMap.end() ; )
		{
			if (it->second.size() != m_numberOfDifs)
			{
				it++;
				continue;
			}

			// create the event and collection
			EVENT::LCEvent *pLCEvent = createEvent();
			EVENT::LCCollection *pCollection = createRawDataCollection();
			pLCEvent->addCollection(pCollection, m_collectionName);

			// create the LCGenericObject collection from the vector of difs
			for (std::vector<unsigned char*>::iterator iv = it->second.begin() ; iv != it->second.end() ; iv++)
			{
				int32_t* idata = (int32_t*) (*iv);
				uint32_t ib0 = 24;
				unsigned char* cdif = (unsigned char*) idata;

				uint32_t dif_id = Decoder::getBufferDIF(&cdif[ib0], 0);

				addDifEntryInCollection(pCollection, &idata[0], idata[SHM_BUFFER_SIZE]/sizeof(uint32_t)+1);
			}

			// create a dqm event
			dqm4hep::DQMLCEvent dqmEvent;
			dqmEvent.setEvent(pLCEvent);

			// send it to the system
			try
			{
				THROW_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventSender->sendEvent(&dqmEvent));
			}
			catch(dqm4hep::StatusCodeException &exception)
			{
				streamlog_out(WARNING) << "Couldn't send event : " << exception.getStatusCode() << std::endl;
			}

			// purge the map
			for (std::vector<unsigned char*>::iterator iv=it->second.begin();iv!=it->second.end();iv++)
				delete (*iv);

			it->second.clear();
			bufferMap.erase(it++);
		}
	}
}

//-------------------------------------------------------------------------------------------------

EVENT::LCEvent *ShmReader::createEvent()
{
	IMPL::LCEventImpl *pLCEventImpl = new IMPL::LCEventImpl();
	pLCEventImpl->setRunNumber(0);  // see in dqm run control service
	pLCEventImpl->setTimeStamp(time(0));
	pLCEventImpl->setDetectorName(m_detectorName);
	pLCEventImpl->setEventNumber(m_currentEventNumber);

	m_currentEventNumber++;

	return pLCEventImpl;
}

//-------------------------------------------------------------------------------------------------

EVENT::LCCollection *ShmReader::createRawDataCollection() const
{
	return new IMPL::LCCollectionVec(EVENT::LCIO::LCGENERICOBJECT);
}

//-------------------------------------------------------------------------------------------------

void ShmReader::addDifEntryInCollection(EVENT::LCCollection *pCollection, int *pBuffer, int bufferSize) const
{
	IMPL::LCGenericObjectImpl *pGenericObjectImpl = new IMPL::LCGenericObjectImpl();

	for(int i=0 ; i<bufferSize ; i++)
		pGenericObjectImpl->setIntVal(i, pBuffer[i]);

	pCollection->addElement(pGenericObjectImpl);
}

} 

