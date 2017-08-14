/// \file HitAnalysisModule.cc
/*
 *
 * HitAnalysisModule.cc source template automatically generated by a class generator
 * Creation date : ven. ao�t 28 2015
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
 * @author Remi Ete, Arnaud Steen
 * @copyright CNRS , IPNL
 */


#include "HitAnalysisModule.h"

// -- dqm4hep headers
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMEvent.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMPlugin.h"
#include "dqm4hep/DQMPluginManager.h"


// -- lcio headers
#include <EVENT/LCCollection.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTOOLS.h>

// -- dqmsdhcal headers
#include <AnalysisTools.h>

using namespace dqm4hep;

namespace dqm_sdhcal
{

// plugin declaration
DQM_PLUGIN_DECL( HitAnalysisModule, "HitAnalysisModule" )

HitAnalysisModule::HitAnalysisModule() :
	DQMTriventModule(),
	m_moduleLogStr("[HitAnalysisModule]")
{
}

//-------------------------------------------------------------------------------------------------

HitAnalysisModule::~HitAnalysisModule()
{
}

//-------------------------------------------------------------------------------------------------
dqm4hep::StatusCode HitAnalysisModule::userReadSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	m_nActiveLayers = 48;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "NActiveLayers", m_nActiveLayers));

	m_cellIDDecoderString = "M:3,S-1:3,I:9,J:9,K-1:6";
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "CellIDDecoderString", m_cellIDDecoderString));

	/*-----------------------------------------------------*/
	m_clusteringSettings.maxTransversal = 1;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "Clustering.MaxTranverseCellID", m_clusteringSettings.maxTransversal));

	m_clusteringSettings.maxLongitudinal = 0;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "Clustering.MaxLongitudinalCellID", m_clusteringSettings.maxLongitudinal));

	m_clusteringSettings.useDistanceInsteadCellID = false;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "Clustering.UseDistanceInsteadCellID", m_clusteringSettings.useDistanceInsteadCellID));

	m_clusteringSettings.maxTransversalDistance = 11.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "Clustering.MaxTranverseDistance", m_clusteringSettings.maxTransversalDistance));

	m_clusteringSettings.maxLongitudinalDistance = 11.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "Clustering.MaxLongitudinalDistance", m_clusteringSettings.maxLongitudinalDistance));

	/*-----------------------------------------------------*/
	m_clusteringHelperSettings.longitudinalDistance = 100.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "ClusteringHelper.LongitudinalDistanceForIsolation", m_clusteringHelperSettings.longitudinalDistance));

	m_clusteringHelperSettings.transversalDistance = 200.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "ClusteringHelper.TransverseDistanceForIsolation", m_clusteringHelperSettings.transversalDistance));

	/*-----------------------------------------------------*/
	// Number of hits in full detector
	//
	//
	m_pInstantRate = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "InstantRate", m_pInstantRate));

	m_pMeanRunRate = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "MeanRunRate", m_pMeanRunRate));

	m_pRateVsClusterProfile = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "RateVsClusterProfile", m_pRateVsClusterProfile));

	m_pRateVsClusterProfileNoClassification = NULL;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "RateVsClusterProfile", m_pRateVsClusterProfileNoClassification));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits0", m_pNHit0));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits1", m_pNHit1));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits2", m_pNHit2));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHitsTotal", m_pNHit));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits0PerLayer", m_pNHit0PerLayer));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits1PerLayer", m_pNHit1PerLayer));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHits2PerLayer", m_pNHit2PerLayer));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
	                 "NumberOfHitsTotalPerLayer", m_pNHitTotPerLayer));


	// Number of hits per layer
	for (unsigned int layerId = 0 ; layerId < m_nActiveLayers ; layerId++)
	{
		DQMParameters parameters;
		parameters["layerId"] = std::to_string(layerId);

		std::string folderPath;
		RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::getAttribute(xmlHandle.FirstChild("monitorElement").Element(), "path", folderPath));
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle, "ChamberHitsMap1", m_layerElementMap[layerId].m_pChamberHitsMap0, parameters));
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle, "ChamberHitsMap2", m_layerElementMap[layerId].m_pChamberHitsMap1, parameters));
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle, "ChamberHitsMap3", m_layerElementMap[layerId].m_pChamberHitsMap2, parameters));
	    
	}

	m_firstLayerCut = 1;// 2;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                 "FirstLayerCut", m_firstLayerCut));

	m_lastLayerCut = 48;//46;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                 "LastLayerCut", m_lastLayerCut));
	m_nMipInLayer = 25;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                 "NMipInLayer", m_nMipInLayer));
	m_nMipMinimum = 40;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                 "NMipMinimum", m_nMipMinimum));

	m_inputCollectionName = "SDHCAL_HIT";
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
	                        "InputCollectionName", m_inputCollectionName));

	/*------ Event classifier settings ------*/
	dqm4hep::TiXmlElement *pEventClassifierElement = xmlHandle.FirstChild("eventClassifier").Element();

	if ( ! pEventClassifierElement )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't found xml element eventClassifier !" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	std::string plugin;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::getAttribute(pEventClassifierElement, "plugin", plugin));

	m_pEventClassifier = dqm4hep::DQMPluginManager::instance()->createPluginClass<EventClassifier>(plugin);

	if ( ! m_pEventClassifier )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't found eventClassifier plugin called : " << plugin );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventClassifier->readSettings(dqm4hep::TiXmlHandle(pEventClassifierElement)));

	/*-----------------------------------------------------*/

	/*------ Event helper settings ------*/
	dqm4hep::TiXmlElement *pEventHelperElement = xmlHandle.FirstChild("eventHelper").Element();

	if ( ! pEventHelperElement )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't found xml element eventHelper !" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::getAttribute(pEventHelperElement, "plugin", plugin));

	m_pEventHelper = dqm4hep::DQMPluginManager::instance()->createPluginClass<EventHelper>(plugin);

	if ( ! m_pEventHelper )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't found eventHelper plugin called : " << plugin );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventHelper->readSettings(dqm4hep::TiXmlHandle(pEventHelperElement)));

	/*-----------------------------------------------------*/



	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::userInitModule()
{

	m_nTrigger = 0;
	m_nSpill = 0;

    m_nParticleWithinRun = 0;
	m_nParticleWithinSpill = 0;
    m_nBeamMuonWithinRun = 0;
	m_nBeamMuonWithinSpill = 0;
    m_nChargedHadronsWithinRun = 0;
	m_nChargedHadronsWithinSpill = 0;
    m_nNeutralHadronsWithinRun = 0;
	m_nNeutralHadronsWithinSpill = 0;
    m_nPhotonsWithinRun = 0;
	m_nPhotonsWithinSpill = 0;
    m_nElectronsWithinRun = 0;
	m_nElectronsWithinSpill = 0;
    m_nOthersWithinRun = 0;
    m_nOthersWithinSpill = 0;
    m_nCosmicMuonsWithinRun = 0;
    m_nCosmicMuonsWithinSpill = 0;
    m_nUndefinedWithinRun = 0;
    m_nUndefinedWithinSpill = 0;
    m_nNoiseWithinRun = 0;
    m_nNoiseWithinSpill = 0;
    // m_eventIntegratedTime = 0;
    // m_spillIntegratedTime = 0;
    // m_totalIntegratedTime = 0;
    // m_timeLastTrigger = 0;
    // m_timeLastSpill = 0;

	// initialize algorithms
	m_clusteringAlgorithm.SetClusterParameterSetting(m_clusteringSettings);
	m_clusteringHelper.SetClusteringHelperParameterSetting(m_clusteringHelperSettings);

	return STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------
dqm4hep::StatusCode HitAnalysisModule::processEvent(EVENT::LCEvent *pLCEvent)
{
	LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Processing physics event no " << pLCEvent->getEventNumber() );
	int nHitProcessedEvent = 0;

	// content management
	caloobject::CaloHitMap caloHitMap;
	caloobject::CaloHitList hits;
	caloobject::CaloClusterList clusters;

	CLHEP::Hep3Vector globalHitShift(0, 0, 0);

	try
	{
		EVENT::LCCollection *pCalorimeterHitCollection = pLCEvent->getCollection(m_inputCollectionName);

		if (NULL == pCalorimeterHitCollection)
			return dqm4hep::STATUS_CODE_SUCCESS;

		// Classify and count particle in event
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pEventClassifier->processEvent(pLCEvent));
		if (m_pEventClassifier->isUndefined())
		{
			m_nUndefinedWithinRun++;
			m_nUndefinedWithinSpill++;
		}
		else if (m_pEventClassifier->isNoisyEvent())
		{
			m_nNoiseWithinRun++;
			m_nNoiseWithinSpill++;
		}
		else if ( m_pEventClassifier->isPhysicsEvent() )
		{
			if ( m_pEventClassifier->getEventType() != EventClassifier::COSMIC_MUON_EVENT)
			{
				m_nParticleWithinRun++;
				m_nParticleWithinSpill++;

				if ( m_pEventClassifier->getEventType() == EventClassifier::BEAM_MUON_EVENT )
				{
					m_nBeamMuonWithinRun++;
					m_nBeamMuonWithinSpill++;
				}

				else if ( m_pEventClassifier->getEventType() == EventClassifier::CHARGED_HAD_SHOWER_EVENT )
				{
					m_nChargedHadronsWithinRun++;
					m_nChargedHadronsWithinSpill++;
				}
				else if ( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_HAD_SHOWER_EVENT )
				{
					m_nNeutralHadronsWithinRun++;
					m_nNeutralHadronsWithinSpill++;
				}
				else if ( m_pEventClassifier->getEventType() == EventClassifier::NEUTRAL_EM_SHOWER_EVENT )
				{
					m_nPhotonsWithinRun++;
					m_nPhotonsWithinSpill++;
				}
				else if ( m_pEventClassifier->getEventType() == EventClassifier::CHARGED_EM_SHOWER_EVENT )
				{
					m_nElectronsWithinRun++;
					m_nElectronsWithinSpill++;
				}
			}
			else
			{
				m_nCosmicMuonsWithinSpill++;
				m_nCosmicMuonsWithinRun++;
			}
		}

		// Find New Trigger/Spill and fill Rates per particle
		LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - findTrigger for Trigger event " << pLCEvent->getEventNumber() << "...");
		StatusCode status = m_pEventHelper->findTrigger(pCalorimeterHitCollection, m_eventParameters);

	// m_eventParameters.dumpParameters();
		if (dqm4hep::STATUS_CODE_SUCCESS != status)
		{
			LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - findTrigger for Trigger event " << pLCEvent->getEventNumber() << " failed with status : " << status);
			LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Going to next Trigger event ");
			return dqm4hep::STATUS_CODE_SUCCESS;
		}
		LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - findTrigger for Trigger event " << pLCEvent->getEventNumber() << "...OK");

		if (m_eventParameters.newTrigger)
		{
	    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr <<  " - New Trigger at time : " << m_eventParameters.timeTrigger*m_pEventHelper->getDAQ_BC_Period() << "s\t spillIntegratedTime : " << m_eventParameters.spillIntegratedTime*m_pEventHelper->getDAQ_BC_Period() << "s");

			if (m_eventParameters.newSpill)
			{
		LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr <<  " - New Spill -  time since last startOfSpill : " <<  (m_eventParameters.timeTrigger - m_eventParameters.timeLastSpill)*m_pEventHelper->getDAQ_BC_Period() << " s.  Last spill Stat: Length : " <<  m_eventParameters.spillIntegratedTime*m_pEventHelper->getDAQ_BC_Period() << "s\t triggers : " << m_nTrigger << "\t particles : " << m_nParticleWithinSpill << "\t totalIntegratedTime: " << m_eventParameters.totalIntegratedTime*m_pEventHelper->getDAQ_BC_Period());

				m_nSpill++;
				m_nTrigger = 0;

				LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr <<  " - Filling Rates..." );
				if (m_eventParameters.lastSpillIntegratedTime != 0)
					RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->fillRates());
				LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr <<  " - Filling Rates...OK" );

				// Reinitialize rates
				m_nUndefinedWithinSpill = 0;
				m_nNoiseWithinSpill = 0;
				m_nCosmicMuonsWithinSpill = 0;
				m_nParticleWithinSpill = 0;
				m_nBeamMuonWithinSpill = 0;
				m_nChargedHadronsWithinSpill = 0;
				m_nNeutralHadronsWithinSpill = 0;
				m_nPhotonsWithinSpill = 0;
				m_nElectronsWithinSpill = 0;
			}
		}

		UTIL::CellIDDecoder<EVENT::CalorimeterHit> cellIDDecoder(m_cellIDDecoderString);

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << "Creating wrapper hits");

		// loop over hits in this event
		for (unsigned int h = 0 ; h < pCalorimeterHitCollection->getNumberOfElements() ; h++)
		{
			EVENT::CalorimeterHit *pCaloHit = dynamic_cast<EVENT::CalorimeterHit*>(pCalorimeterHitCollection->getElementAt(h));

			if (NULL == pCaloHit)
				continue;

			int cellID[3];
			cellID[0] = cellIDDecoder(pCaloHit)["I"];
			cellID[1] = cellIDDecoder(pCaloHit)["J"];
			cellID[2] = cellIDDecoder(pCaloHit)["K-1"];

			if ( cellID[2] >= m_nActiveLayers )
			{
				LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Wrong number of layer in your configuration file!");
				LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Found a hit in layer " << cellID[2] << " - Last layer in xml configuration file is " << m_nActiveLayers);
				continue;
			}

			CLHEP::Hep3Vector positionVector(
			  pCaloHit->getPosition()[0],
			  pCaloHit->getPosition()[1],
			  pCaloHit->getPosition()[2] );

			caloobject::CaloHit *pWrapperHit = new caloobject::CaloHit(
			  cellID,
			  positionVector,
			  pCaloHit->getEnergy(),
			  pCaloHit->getTime(),
			  globalHitShift);

			caloHitMap[ cellID[2] ].push_back(pWrapperHit);
			hits.push_back(pWrapperHit);
		}

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << "Creating intra layer clusters");

		for (caloobject::CaloHitMap::iterator iter = caloHitMap.begin(), endIter = caloHitMap.end() ;
		     iter != endIter ; ++iter)
			m_clusteringAlgorithm.Run(iter->second, clusters);

		std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);

	Float_t rate =  m_nParticleWithinSpill/m_eventParameters.spillIntegratedTime*m_pEventHelper->getDAQ_BC_Period();
		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Number of clusters : " << clusters.size());
	LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Number of particleInSpill : " << m_nParticleWithinSpill);
	LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " SpillIntegratedTime : " << m_eventParameters.spillIntegratedTime*m_pEventHelper->getDAQ_BC_Period());
	LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Rate : " << rate );

		// Fill Number Of Clusters vs Rate
		// if (m_eventParameters.spillIntegratedTime > 0 && m_pEventClassifier->isPhysicsEvent() && m_pEventClassifier->getEventType() != EventClassifier::COSMIC_MUON_EVENT)
		if (m_eventParameters.spillIntegratedTime > 0 && m_pEventClassifier->isPhysicsEvent())
			m_pRateVsClusterProfile->get<TProfile>()->Fill(m_nParticleWithinSpill / m_eventParameters.lastSpillIntegratedTime * m_pEventHelper->getDAQ_BC_Period(), clusters.size());
		if (m_eventParameters.spillIntegratedTime > 0)
			m_pRateVsClusterProfileNoClassification->get<TProfile>()->Fill(m_nParticleWithinSpill / m_eventParameters.lastSpillIntegratedTime * m_pEventHelper->getDAQ_BC_Period(), clusters.size());

		int nHit0 = 0;
		int nHit1 = 0;
		int nHit2 = 0;
		int nMip = 0;
		int nHit0Layer = 0;
		int nHit1Layer = 0;
		int nHit2Layer = 0;
		int nMipLayer  = 0;

		for (caloobject::CaloClusterList::const_iterator clusterIter = clusters.begin(), clusterEndIter = clusters.end(); clusterEndIter != clusterIter ; ++clusterIter)
		{
			unsigned int nHitProcessed = 0; 										// Number of hit processed in the current cluster

			int layerId = (*clusterIter)->getLayerID();
			if (layerId >= m_nActiveLayers)
			{
				LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " Found a cluster after last layer... : layer=" << layerId << " maxLayer= " << m_nActiveLayers);
				continue;
			}

			for (caloobject::CaloHitList::const_iterator hitIter = (*clusterIter)->getHits().begin(), hitEndIter = (*clusterIter)->getHits().end(); hitEndIter != hitIter; ++hitIter)
			{
				uint32_t hitWeight = 0;
				int hitThreshold = (*hitIter)->getEnergy();

				if (hitThreshold == 3)
				{
					hitWeight = 15;
					nHit2Layer++;

		    // if ( m_pEventClassifier->isPhysicsEvent() && m_pEventClassifier->getEventType() != EventClassifier::COSMIC_MUON_EVENT)
		    m_layerElementMap.at(layerId).m_pChamberHitsMap2->get<TH2>()->Fill((*hitIter)->getCellID()[0], (*hitIter)->getCellID()[1]);
				}
		else if (hitThreshold == 2) // || hitThreshold == 3)
				{
					hitWeight = 3;
					nHit1Layer++;

		    // if ( m_pEventClassifier->isPhysicsEvent() && m_pEventClassifier->getEventType() != EventClassifier::COSMIC_MUON_EVENT)
		    m_layerElementMap.at(layerId).m_pChamberHitsMap1->get<TH2>()->Fill((*hitIter)->getCellID()[0], (*hitIter)->getCellID()[1]);

				}
		else if (hitThreshold == 1 ) //|| hitThreshold == 2 || hitThreshold == 3)
				{
					hitWeight = 1;
					nHit0Layer++;
		    // if ( m_pEventClassifier->isPhysicsEvent() && m_pEventClassifier->getEventType() != EventClassifier::COSMIC_MUON_EVENT)
		    m_layerElementMap.at(layerId).m_pChamberHitsMap0->get<TH2>()->Fill((*hitIter)->getCellID()[0], (*hitIter)->getCellID()[1]);
				}
				else
				{
					LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " Wait what? Hit with threshold : " << hitThreshold << " in layer " << layerId);
					continue;
				}

				nMipLayer += hitWeight;
				nHitProcessedEvent++; 	// Total number of hits processed in current event
				nHitProcessed++;			  // Total number of hits processed in current cluster

				// if (hitEndIter == hitIter + 1)
				// {
				// 	int numberOfHit = (*clusterIter)->getHits().size(); // Number of hits in the current cluster
				// 	LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Processed : " << nHitProcessed << "/" << numberOfHit << " hits in current cluster of layer " << layerId << "\t nHitTotLayer : " << nHit0Layer + nHit1Layer + nHit2Layer  << "\t nHitLayer0 : " << nHit0Layer << "\t nHitLayer1 : " << nHit1Layer << "\t nHitLayer2 : " << nHit2Layer);
				// }
			}

			// uint32_t firstLayer = m_nActiveLayers, lastLayer = 1;
			// for (uint32_t i = 1; i < m_nActiveLayers + 1; ++i)
			// {
			// 	if (nMipLayer > m_nMipInLayer && i < firstLayer) firstLayer = i;
			// 	if (nMipLayer > m_nMipInLayer && i > lastLayer) lastLayer = i;
			// }

			// if (firstLayer >= m_firstLayerCut && lastLayer <= m_lastLayerCut)
			// {
			// 	for (uint32_t i = 0; i < m_nActiveLayers; ++i)
			// 	{
			// 		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Adding hits in layer : " << i << "\t hits0 : " << nHit0Layer << "\t hits1 : " << nHit1Layer<< "\t hits2 : " << nHit2Layer);

			// Stack hits from multiple clusters in the same layer before filling histograms
			caloobject::CaloClusterList::const_iterator nextClusterIter;

			if (clusterEndIter != clusterIter + 1)
				nextClusterIter = std::next(clusterIter, 1);
			else
				nextClusterIter = clusterIter;

			// Fill hit histograms only if last cluster of current layer or last cluster of the event
			if ( (*nextClusterIter)->getLayerID() != layerId || nextClusterIter == clusterIter)
			{
				if (0 != nHit0Layer)
				{
					// m_layerElementMap[layerId].m_pNHit0Layer->get<TH1>()->Fill(nHit0Layer);
					m_pNHit0PerLayer->get<TH1>()->SetBinContent(layerId, m_pNHit0PerLayer->get<TH1>()->GetBinContent(layerId) + nHit0Layer);
				}
				if (0 != nHit1Layer)
				{
					// m_layerElementMap[layerId].m_pNHit1Layer->get<TH1>()->Fill(nHit1Layer);
					m_pNHit1PerLayer->get<TH1>()->SetBinContent(layerId, m_pNHit1PerLayer->get<TH1>()->GetBinContent(layerId) +  nHit1Layer);
				}
				if (0 != nHit2Layer)
				{
					// m_layerElementMap[layerId].m_pNHit2Layer->get<TH1>()->Fill(nHit2Layer);
					m_pNHit2PerLayer->get<TH1>()->SetBinContent(layerId, m_pNHit2PerLayer->get<TH1>()->GetBinContent(layerId) +  nHit2Layer);
				}
				// m_layerElementMap[layerId].m_pNHitTotLayer->get<TH1>()->Fill(nHit0Layer + nHit1Layer + nHit2Layer);
				m_pNHitTotPerLayer->get<TH1>()->SetBinContent(layerId, m_pNHitTotPerLayer->get<TH1>()->GetBinContent(layerId) +  nHit0Layer + nHit1Layer + nHit2Layer);

				nHit0 += nHit0Layer;
				nHit1 += nHit1Layer;
				nHit2 += nHit2Layer;
				nMip += nMipLayer;
				nHit0Layer = 0;
				nHit1Layer = 0;
				nHit2Layer = 0;
			}
			// 	}
			// }
		} // Loop on TrackingClusters

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " Total hits processed : " << nHitProcessedEvent << "\t hits(0+1+2) : " << nHit0 + nHit1 + nHit2 << "\t hits0 : " << nHit0 << "\t hits1 : " << nHit1 << "\t hits2 : " << nHit2 << "\t nMip : " << nMip);

	// if (nMip > m_nMipMinimum)
	// {
			m_pNHit0->get<TH1>()->Fill(nHit0);
			m_pNHit1->get<TH1>()->Fill(nHit1);
			m_pNHit2->get<TH1>()->Fill(nHit2);
			m_pNHit->get<TH1>()->Fill(nHit0 + nHit1 + nHit2);
	// }
	}
	catch (EVENT::DataNotAvailableException &exception)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << "Caught EVENT::DataNotAvailableException : " << exception.what());
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger ,  m_moduleLogStr << "Skipping event" );
		this->clearEventContents(hits, clusters);
		return STATUS_CODE_SUCCESS;
	}
	catch (...)
	{
		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,  m_moduleLogStr << "Caught unknown exception !");
		this->clearEventContents(hits, clusters);
		return STATUS_CODE_FAILURE;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::fillRates()
{
	std::stringstream instantRate;

    if (m_eventParameters.lastSpillIntegratedTime == 0){
      LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - spillIntegratedTime is null, not filling rates" );
      return dqm4hep::STATUS_CODE_SUCCESS;
    }
    float spillTimeInSeconds = m_eventParameters.lastSpillIntegratedTime*m_pEventHelper->getDAQ_BC_Period();

	instantRate << "*****  Instant Rate  *****\n";
	instantRate << "*****  Physics  *****\n";
    instantRate << "Current spill Length : " << spillTimeInSeconds << "s\n";
    instantRate << " - Particles    	 : " << m_nParticleWithinSpill << "/spill\t" << m_nParticleWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Beam muons      : " 	     << m_nBeamMuonWithinSpill << "/spill\t" << m_nBeamMuonWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Charged hadrons : " 	     << m_nChargedHadronsWithinSpill << "/spill\t" << m_nChargedHadronsWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Neutral hadrons : " 	     << m_nNeutralHadronsWithinSpill << "/spill\t" << m_nNeutralHadronsWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Photons         : " 	     << m_nPhotonsWithinSpill << "/spill\t" << m_nPhotonsWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Electrons       : " 	     << m_nElectronsWithinSpill	<< "/spill\t" << m_nElectronsWithinSpill / spillTimeInSeconds << "/s" << "\n";
	instantRate << "*****  Non Physics  *****\n";
    instantRate << " - Undefined       : " 	     << m_nUndefinedWithinSpill << "/spill\t" << m_nUndefinedWithinSpill / spillTimeInSeconds << "/s" << "\n";
    instantRate << " - Noise           : " 	     << m_nNoiseWithinSpill << "/spill\t" << m_nNoiseWithinSpill / spillTimeInSeconds	<< "/s" << "\n";
    instantRate << " - Cosmic muons    : " 	     << m_nCosmicMuonsWithinSpill << "/spill\t" << m_nCosmicMuonsWithinSpill / spillTimeInSeconds << "/s" << "\n";
	instantRate << "*********************************";

    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,"*****  Instant Rate  *****\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,"*****  Physics  *****\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,"Current spill Length : " << spillTimeInSeconds << "s\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Particles    	 : " << m_nParticleWithinSpill << "/spill\t" << m_nParticleWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Beam muons      : " 	     << m_nBeamMuonWithinSpill << "/spill\t" << m_nBeamMuonWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Charged hadrons : " 	     << m_nChargedHadronsWithinSpill << "/spill\t" << m_nChargedHadronsWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Neutral hadrons : " 	     << m_nNeutralHadronsWithinSpill << "/spill\t" << m_nNeutralHadronsWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Photons         : " 	     << m_nPhotonsWithinSpill << "/spill\t" << m_nPhotonsWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Electrons       : " 	     << m_nElectronsWithinSpill	<< "/spill\t" << m_nElectronsWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,"*****  Non Physics  *****\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Undefined       : " 	     << m_nUndefinedWithinSpill << "/spill\t" << m_nUndefinedWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Noise           : " 	     << m_nNoiseWithinSpill << "/spill\t" << m_nNoiseWithinSpill / spillTimeInSeconds	<< "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ," - Cosmic muons    : " 	     << m_nCosmicMuonsWithinSpill << "/spill\t" << m_nCosmicMuonsWithinSpill / spillTimeInSeconds << "/s" << "\n");
    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger ,"*********************************");

    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - \n " << instantRate.str() );
	m_pInstantRate->get< dqm4hep::TScalarString >()->Set( instantRate.str() );

	std::stringstream meanRunRate;

    if (m_eventParameters.totalIntegratedTime == 0 && m_nSpill == 0){
      LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - totalIntegratedTime is null, not filling rates" );
      return dqm4hep::STATUS_CODE_SUCCESS;
    }
    float totalTimeInSeconds = m_eventParameters.totalIntegratedTime*m_pEventHelper->getDAQ_BC_Period();

	meanRunRate << "*****  Mean Run Rate  *****\n";
	meanRunRate << "*****  Physics  *****\n";
    meanRunRate << "Current run length : " <<  totalTimeInSeconds << "s\n";
	meanRunRate << " - Particles    	 : " << m_nParticleWithinRun << "/Run\t"
	            << m_nParticleWithinRun / m_nSpill << "/Spill\t"
		<< m_nParticleWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Beam muons      : " << m_nBeamMuonWithinRun 	<< "/Run\t"
	            << m_nBeamMuonWithinRun / m_nSpill << "/Spill\t"
		<< m_nBeamMuonWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Charged hadrons : " << m_nChargedHadronsWithinRun << "/Run\t"
	            << m_nChargedHadronsWithinRun / m_nSpill << "/Spill\t"
		<< m_nChargedHadronsWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Neutral hadrons : " << m_nNeutralHadronsWithinRun << "/Run\t"
	            << m_nNeutralHadronsWithinRun / m_nSpill << "/Spill\t"
		<< m_nNeutralHadronsWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Photons         : " << m_nPhotonsWithinRun << "/Run\t"
	            << m_nPhotonsWithinRun / m_nSpill << "/Spill\t"
		<< m_nPhotonsWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Electrons       : " << m_nElectronsWithinRun	<< "/Run\t"
	            << m_nElectronsWithinRun / m_nSpill << "/Spill\t"
		<< m_nElectronsWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << "*****  Non Physics  *****\n";
	meanRunRate << " - Undefined       : " << m_nUndefinedWithinRun << "/Run\t"
	            << m_nUndefinedWithinRun / m_nSpill << "/Spill\t"
		<< m_nUndefinedWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Noise           : " << m_nNoiseWithinRun	<< "/Run\t"
	            << m_nNoiseWithinRun / m_nSpill << "/Spill\t"
		<< m_nNoiseWithinRun / totalTimeInSeconds << "/s\n";

	meanRunRate << " - Cosmic muons    : " << m_nCosmicMuonsWithinRun << "/Run\t"
	            << m_nCosmicMuonsWithinRun / m_nSpill << "/Spill\t"
		<< m_nCosmicMuonsWithinRun / totalTimeInSeconds << "/s\n";
	meanRunRate << "*********************************";

    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - \n "<< meanRunRate.str() );
	m_pMeanRunRate->get< dqm4hep::TScalarString >()->Set( meanRunRate.str() );

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::startOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::endOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::startOfRun(DQMRun * pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::endOfRun(DQMRun * pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode HitAnalysisModule::endModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

void HitAnalysisModule::clearEventContents(caloobject::CaloHitList & hits, caloobject::CaloClusterList & clusters)
{
	for_each(hits.begin(), hits.end(), [] (caloobject::CaloHit * pCaloHit) { delete pCaloHit; });
	for_each(clusters.begin(), clusters.end(), [] (caloobject::CaloCluster * pCluster) { delete pCluster; });

	hits.clear();
	clusters.clear();
}

//-------------------------------------------------------------------------------------------------

void HitAnalysisModule::resetElements()
{
	for (std::map<unsigned int, LayerElements>::iterator iter = m_layerElementMap.begin(), endIter = m_layerElementMap.end() ;
	     endIter != iter ; ++iter)
	{
		iter->second.m_pNHit0Layer->reset();
		iter->second.m_pNHit1Layer->reset();
		iter->second.m_pNHit2Layer->reset();
		iter->second.m_pNHitTotLayer->reset();
	}
	m_pNHit0->reset();
	m_pNHit1->reset();
	m_pNHit2->reset();
	m_pNHit->reset();
}

}

