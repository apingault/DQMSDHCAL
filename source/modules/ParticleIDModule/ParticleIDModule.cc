/*
 *
 * ParticleIDModule.cc source template automatically generated by a class generator
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


#include "ParticleIDModule.h"
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

DQM_PLUGIN_DECL( ParticleIDModule , "ParticleIDModule" )

//-------------------------------------------------------------------------------------------------

ParticleIDModule::ParticleIDModule() :
		DQMTriventModule(),
		m_pEventClassifier(NULL),
		m_nNoiseWithinRun(0),
		m_nCosmicMuonsWithinRun(0),
		m_nBeamMuonWithinRun(0),
		m_nChargedHadronsWithinRun(0),
		m_nNeutralHadronsWithinRun(0),
		m_nElectronsWithinRun(0),
		m_nPhotonsWithinRun(0),
		m_nOthersWithinRun(0)
{
}

//-------------------------------------------------------------------------------------------------

ParticleIDModule::~ParticleIDModule()
{

}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::userInitModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::userReadSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"CaloHitCollectionNames", m_caloHitCollectionNames));

	/*------ Monitor element booking ------*/
	m_pParticleIDSummary = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"ParticleIDSummary", m_pParticleIDSummary));

	m_pNHitNoise = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitNormal", "Noise", m_pNHitNoise));
	m_pNHitNoise->setTitle( m_pNHitNoise->getTitle() + " [Noise]" );

	m_pNHitCosmicMuons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitTight", "CosmicMuons", m_pNHitCosmicMuons));
	m_pNHitCosmicMuons->setTitle( m_pNHitCosmicMuons->getTitle() + " [CosmicMuons]" );

	m_pNHitBeamMuons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitTight", "BeamMuons", m_pNHitBeamMuons));
	m_pNHitBeamMuons->setTitle( m_pNHitBeamMuons->getTitle() + " [BeamMuons]" );

	m_pNHitNeutralHadrons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "NeutralHadrons", m_pNHitNeutralHadrons));
	m_pNHitNeutralHadrons->setTitle( m_pNHitNeutralHadrons->getTitle() + " [NeutralHadrons]" );

	m_pNHitChargedHadrons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "ChargedHadrons", m_pNHitChargedHadrons));
	m_pNHitChargedHadrons->setTitle( m_pNHitChargedHadrons->getTitle() + " [ChargedHadrons]" );

	m_pNHitElectrons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "Electrons", m_pNHitElectrons));
	m_pNHitElectrons->setTitle( m_pNHitElectrons->getTitle() + " [Electrons]" );

	m_pNHitPhotons = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "Photons", m_pNHitPhotons));
	m_pNHitPhotons->setTitle( m_pNHitPhotons->getTitle() + " [Photons]" );

	m_pNHitOthers = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "Others", m_pNHitOthers));
	m_pNHitOthers->setTitle( m_pNHitOthers->getTitle() + " [Others]" );

	m_pNHitAll = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NHitWide", "All", m_pNHitAll));
	m_pNHitAll->setTitle( m_pNHitAll->getTitle() + " [All]" );

	/*------ Event classifier settings ------*/
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

dqm4hep::StatusCode ParticleIDModule::processEvent(EVENT::LCEvent *pLCEvent)
{
	LOG4CXX_INFO( dqm4hep::dqmMainLogger , "Processing physics event no " << pLCEvent->getEventNumber() );

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventClassifier->processEvent(pLCEvent));

	unsigned int nHits(0);
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->getNHits(pLCEvent, nHits));

	m_pNHitAll->get<TH1>()->Fill( nHits );

	if( m_pEventClassifier->isNoisyEvent() )
	{
		m_nNoiseWithinRun++;
		m_pNHitNoise->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::COSMIC_MUON_EVENT )
	{
		m_nCosmicMuonsWithinRun++;
		m_pNHitCosmicMuons->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::BEAM_MUON_EVENT )
	{
		m_nBeamMuonWithinRun++;
		m_pNHitBeamMuons->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::CHARGED_HAD_SHOWER_EVENT )
	{
		m_nChargedHadronsWithinRun++;
		m_pNHitChargedHadrons->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_HAD_SHOWER_EVENT )
	{
		m_nNeutralHadronsWithinRun++;
		m_pNHitNeutralHadrons->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_EM_SHOWER_EVENT )
	{
		m_nPhotonsWithinRun++;
		m_pNHitPhotons->get<TH1>()->Fill( nHits );
	}
	else if( m_pEventClassifier->getEventType() == EventClassifier::CHARGED_EM_SHOWER_EVENT )
	{
		m_nElectronsWithinRun++;
		m_pNHitElectrons->get<TH1>()->Fill( nHits );
	}
	else
	{
		m_nOthersWithinRun++;
		m_pNHitOthers->get<TH1>()->Fill( nHits );
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::getNHits(EVENT::LCEvent *pLCEvent, unsigned int &nHits)
{
	try
	{
		for(auto iter = m_caloHitCollectionNames.begin(), endIter = m_caloHitCollectionNames.end() ;
				endIter != iter ; ++iter)
		{
		  if  (std::find(pLCEvent->getCollectionNames()->begin(), pLCEvent->getCollectionNames()->end(), *iter) == pLCEvent->getCollectionNames()->end())
		    continue;

			EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

			if( pLCCollection->getTypeName() != EVENT::LCIO::CALORIMETERHIT )
			{
				LOG4CXX_WARN( dqm4hep::dqmMainLogger , "Collection '" << *iter << "' is not a calorimeter hit collection ! Skipping collection ..." );
				continue;
			}

			nHits += pLCCollection->getNumberOfElements();
		}
	}
	catch(EVENT::DataNotAvailableException &exception)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "ParticleIDModule::getNHits : " << exception.what() );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::fillSummary()
{
	std::stringstream summary;

	summary << "*****  Particle ID summary  *****\n";
	summary << " - Noise           : " << m_nNoiseWithinRun << "\n";
	summary << " - Cosmic muons    : " << m_nCosmicMuonsWithinRun << "\n";
	summary << " - Beam muons      : " << m_nBeamMuonWithinRun << "\n";
	summary << " - Charged hadrons : " << m_nChargedHadronsWithinRun << "\n";
	summary << " - Neutral hadrons : " << m_nNeutralHadronsWithinRun << "\n";
	summary << " - Photons         : " << m_nPhotonsWithinRun << "\n";
	summary << " - Electrons       : " << m_nElectronsWithinRun << "\n";
	summary << " - Others          : " << m_nOthersWithinRun << "\n";
	summary << "*********************************";

	m_pParticleIDSummary->get< dqm4hep::TScalarString >()->Set( summary.str() );

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::startOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::endOfCycle()
{
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillSummary());
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::startOfRun(dqm4hep::DQMRun *const /*pRun*/)
{
	m_nNoiseWithinRun = 0;
	m_nCosmicMuonsWithinRun = 0;
	m_nBeamMuonWithinRun = 0;
	m_nChargedHadronsWithinRun = 0;
	m_nNeutralHadronsWithinRun = 0;
	m_nElectronsWithinRun = 0;
	m_nPhotonsWithinRun = 0;
	m_nOthersWithinRun = 0;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::endOfRun(dqm4hep::DQMRun *const /*pRun*/)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode ParticleIDModule::endModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

} 

