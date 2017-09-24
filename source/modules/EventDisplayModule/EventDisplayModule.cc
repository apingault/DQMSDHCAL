/*
 *
 * EventDisplayModule.cc source template automatically generated by a class generator
 * Creation date : jeu. mars 10 2016
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


#include "EventDisplayModule.h"
#include "AnalysisTools.h"

//-- lcio headers
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// -- dqm4hep headers
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMEvent.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMPluginManager.h"

namespace dqm_sdhcal
{

DQM_PLUGIN_DECL( EventDisplayModule , "EventDisplayModule" )

//-------------------------------------------------------------------------------------------------

EventDisplayModule::EventDisplayModule() :
		DQMTriventModule(),
		m_pEventClassifier(NULL)
{

}

//-------------------------------------------------------------------------------------------------

EventDisplayModule::~EventDisplayModule() 
{

}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::userInitModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::userReadSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"InputCaloHitCollections", m_inputCaloHitCollections, [] (const dqm4hep::StringVector &vec) { return ! vec.empty(); }));

	const unsigned int nCollections = m_inputCaloHitCollections.size();

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"ColorWeightList", m_colorWeightList, [&] (const dqm4hep::IntVector &vec) { return vec.size() == nCollections; }));

	int markerColor3D = kBlack;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"MarkerColor3D", markerColor3D));

	int markerSize3D = 1;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"MarkerSize3D", markerSize3D));

	int markerStyle3D = 21;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"MarkerStyle3D", markerStyle3D));

	/*------ Monitor element booking ------*/
	std::vector<std::string> directories =
	  { "All", "Noise", "CosmicMuons", "BeamMuons",
		  "ChargedHadrons", "NeutralHadrons", "EmShowers" ,
		  "Others"};

	const std::string baseDirectory("/");
//	dqm4hep::DQMModuleApi::mkdir(this, baseDirectory);

	for(auto iter = directories.begin(), endIter = directories.end();
			endIter != iter ; ++iter)
	{
		// enter global directory
		dqm4hep::DQMModuleApi::cd(this, baseDirectory);

		// create and enter per particle directory
		dqm4hep::DQMModuleApi::mkdir(this, *iter);
		dqm4hep::DQMModuleApi::cd(this, *iter);

		dqm4hep::DQMMonitorElementPtrList list;
		DisplayElements elements;

		elements.m_pEventDisplay3D = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"EventDisplay3D", elements.m_pEventDisplay3D));
		list.push_back(elements.m_pEventDisplay3D);

		elements.m_pEventDisplay3D->get<TH3>()->SetMarkerColor(markerColor3D);
		elements.m_pEventDisplay3D->get<TH3>()->SetMarkerSize(markerSize3D);
		elements.m_pEventDisplay3D->get<TH3>()->SetMarkerStyle(markerStyle3D);

		// create and enter profiles directory
		dqm4hep::DQMModuleApi::mkdir(this, "Profiles");
		dqm4hep::DQMModuleApi::cd(this, "Profiles");

		elements.m_pLastProfileZX = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"LastProfileZX", elements.m_pLastProfileZX));
		list.push_back(elements.m_pLastProfileZX);

		elements.m_pLastProfileZY = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"LastProfileZY", elements.m_pLastProfileZY));
		list.push_back(elements.m_pLastProfileZY);

		elements.m_pLastProfileXY = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"LastProfileXY", elements.m_pLastProfileXY));
		list.push_back(elements.m_pLastProfileXY);

		elements.m_pCycleStackedProfileZX = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"CycleStackedProfileZX", elements.m_pCycleStackedProfileZX));
		list.push_back(elements.m_pCycleStackedProfileZX);

		elements.m_pCycleStackedProfileZY = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"CycleStackedProfileZY", elements.m_pCycleStackedProfileZY));
		list.push_back(elements.m_pCycleStackedProfileZY);

		elements.m_pCycleStackedProfileXY = NULL;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"CycleStackedProfileXY", elements.m_pCycleStackedProfileXY));
		list.push_back(elements.m_pCycleStackedProfileXY);

		for(auto lIter = list.begin(), endLIter = list.end() ;
				endLIter != lIter ; ++lIter)
			(*lIter)->setTitle( (*lIter)->getTitle() + " [" + *iter + "]" );

		m_displayElementsMap[ *iter ] = elements;
	}

	/*-----------------------------------------------------*/

	dqm4hep::TiXmlElement *pEventClassifierElement = xmlHandle.FirstChild("eventClassifier").Element();

	if( ! pEventClassifierElement )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Couldn't found xml element eventClassifier !" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	std::string plugin;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::getAttribute(pEventClassifierElement, "plugin", plugin));

	m_pEventClassifier = dqm4hep::DQMPluginManager::instance()->createPluginClass<EventClassifier>(plugin);

	if( ! m_pEventClassifier )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Couldn't found eventClassifier plugin called : " << plugin );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventClassifier->readSettings(dqm4hep::TiXmlHandle(pEventClassifierElement)));

	/*-----------------------------------------------------*/

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::processEvent(EVENT::LCEvent *pLCEvent)
{
	LOG4CXX_INFO( dqm4hep::dqmMainLogger , "Processing physics event no " << pLCEvent->getEventNumber() );

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventClassifier->processEvent(pLCEvent));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["All"]));
	
	if( m_pEventClassifier->isNoisyEvent() )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["Noise"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	if( m_pEventClassifier->getEventType() == EventClassifier::COSMIC_MUON_EVENT )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["CosmicMuons"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	if( m_pEventClassifier->getEventType() == EventClassifier::BEAM_MUON_EVENT )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["BeamMuons"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	if( m_pEventClassifier->getEventType() == EventClassifier::CHARGED_HAD_SHOWER_EVENT )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["ChargedHadrons"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	if( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_HAD_SHOWER_EVENT )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["NeutralHadrons"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	// both neutral and charged EM showers
	if( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_EM_SHOWER_EVENT
	 || m_pEventClassifier->getEventType() == EventClassifier::CHARGED_EM_SHOWER_EVENT )
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["EmShowers"]));
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	// fill what remains ...
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillElements(pLCEvent, m_displayElementsMap["Others"]));
	return dqm4hep::STATUS_CODE_SUCCESS;


	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::fillElements(EVENT::LCEvent *pLCEvent, DisplayElements &elements)
{


	const unsigned int nCollections = m_inputCaloHitCollections.size();

	try
	{
		for(unsigned int c=0 ; c<nCollections ; c++)
		{
			const std::string collectionName(m_inputCaloHitCollections.at(c));
			const int colorWeight(m_colorWeightList.at(c));

			try
			{
				EVENT::LCCollection *pCalorimeterHitCollection = pLCEvent->getCollection(collectionName);
				UTIL::CellIDDecoder<EVENT::CalorimeterHit> decoder(pCalorimeterHitCollection);

				LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , "Processing collection " << collectionName << ", n elts = " << pCalorimeterHitCollection->getNumberOfElements() );
				elements.m_pEventDisplay3D->reset();
				elements.m_pLastProfileZX->reset();
				elements.m_pLastProfileZY->reset();
				elements.m_pLastProfileXY->reset();

				// loop over hits in this event
				for(unsigned int h=0 ; h<pCalorimeterHitCollection->getNumberOfElements() ; h++)
				{
					EVENT::CalorimeterHit *pCaloHit = dynamic_cast<EVENT::CalorimeterHit*>(pCalorimeterHitCollection->getElementAt(h));

					if(NULL == pCaloHit)
						continue;

					const float x(pCaloHit->getPosition()[0]);
					const float y(pCaloHit->getPosition()[1]);
					const float z(pCaloHit->getPosition()[2]);

					elements.m_pEventDisplay3D->get<TH3>()->Fill(z, x, y, colorWeight);

					elements.m_pLastProfileZX->get<TH2>()->Fill(z, x);
					elements.m_pLastProfileZY->get<TH2>()->Fill(z, y);
					elements.m_pLastProfileXY->get<TH2>()->Fill(x, y);

					elements.m_pCycleStackedProfileZX->get<TH2>()->Fill(z, x);
					elements.m_pCycleStackedProfileZY->get<TH2>()->Fill(z, y);
					elements.m_pCycleStackedProfileXY->get<TH2>()->Fill(x, y);
				}

			}
			catch(EVENT::DataNotAvailableException &exception)
			{
				LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Caught EVENT::DataNotAvailableException : " << exception.what() );
				continue;
			}
		}
	}
	catch(...)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Caught unknown exception !" );
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::startOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::endOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::startOfRun(dqm4hep::DQMRun *const pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::endOfRun(dqm4hep::DQMRun *const pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode EventDisplayModule::endModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

} 

