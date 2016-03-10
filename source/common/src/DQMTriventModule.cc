  /// \file DQMTriventModule.cc
/*
 *
 * DQMTriventModule.cc source template automatically generated by a class generator
 * Creation date : mar. mars 8 2016
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


#include "DQMTriventModule.h"
#include "Streamout.h"
#include "Trivent.h"

#include "dqm4hep/DQMLogging.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMEvent.h"

#include "Exceptions.h"

namespace dqm_sdhcal
{

DQMTriventModule::DQMTriventModule() :
		DQMAnalysisModule()
{
	m_pStreamout = new Streamout();
	m_pTrivent = new Trivent();

	m_pTrivent->addListener(this);
}

DQMTriventModule::~DQMTriventModule() 
{
	delete 	m_pStreamout;
	delete 	m_pTrivent;
}

dqm4hep::StatusCode DQMTriventModule::readSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	m_shouldProcessStreamout = true;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"ShouldProcessStreamout", m_shouldProcessStreamout));

	if(m_shouldProcessStreamout)
	{
		dqm4hep::TiXmlElement *pXmlElement = xmlHandle.FirstChild("Streamout").Element();

		if( ! pXmlElement )
			return dqm4hep::STATUS_CODE_NOT_FOUND;

		dqm4hep::TiXmlHandle streamoutHandle(pXmlElement);

		std::string inputCollectionName = "RU_XDAQ";
		RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(streamoutHandle,
				"InputCollectionName", inputCollectionName));

		std::string outputCollectionName = "DHCALRawHits";
		RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(streamoutHandle,
				"OutputCollectionName", outputCollectionName));

		unsigned int xdaqShift = 24;
		RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(streamoutHandle,
				"XDaqShift", xdaqShift));

		m_pStreamout->setInputCollectionName(inputCollectionName);
		m_pStreamout->setOutputCollectionName(outputCollectionName);
		m_pStreamout->setXDaqShift(xdaqShift);
	}

	dqm4hep::TiXmlElement *pXmlElement = xmlHandle.FirstChild("Trivent").Element();

	if( ! pXmlElement )
		return dqm4hep::STATUS_CODE_NOT_FOUND;

	dqm4hep::TiXmlHandle triventHandle(pXmlElement);

	// forward parsing to Trivent
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pTrivent->readSettings(triventHandle));
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->userReadSettings(xmlHandle));

	return dqm4hep::STATUS_CODE_SUCCESS;
}

dqm4hep::StatusCode DQMTriventModule::initModule()
{
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pTrivent->init());
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->userInitModule());

	return dqm4hep::STATUS_CODE_SUCCESS;
}

dqm4hep::StatusCode DQMTriventModule::processEvent(dqm4hep::DQMEvent *const pEvent)
{
	EVENT::LCEvent *pLCEvent = pEvent->getEvent<EVENT::LCEvent>();

	if(NULL == pLCEvent)
		return dqm4hep::STATUS_CODE_FAILURE;

	try
	{
		// process Streamout if needed
		if(m_shouldProcessStreamout)
			THROW_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pStreamout->processEvent(pLCEvent));

		THROW_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pTrivent->processEvent(pLCEvent));
	}
	catch(dqm4hep::StatusCodeException &exception)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Caught StatusCodeException : " << exception.toString() );
		return exception.getStatusCode();
	}
	catch(EVENT::DataNotAvailableException &exception)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Caught EVENT::DataNotAvailableException : " << exception.what() );
		return dqm4hep::STATUS_CODE_SUCCESS;
	}
	catch(...)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Caught unknown exception !" );
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

} 
