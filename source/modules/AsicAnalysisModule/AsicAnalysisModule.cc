  /// \file AsicAnalysisModule.cc
/*
 *
 * AsicAnalysisModule.cc source template automatically generated by a class generator
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


#include "AsicAnalysisModule.h"

// -- dqm4hep headers
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMEvent.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMPlugin.h"
#include "dqm4hep/DQMPluginManager.h"

// -- std headers
#include <iostream>
#include <fstream>

//-- lcio headers
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

using namespace dqm4hep;

namespace dqm_sdhcal
{

// plugin declaration
DQM_PLUGIN_DECL( AsicAnalysisModule, "AsicAnalysisModule" )

AsicAnalysisModule::AsicAnalysisModule() :
	DQMTriventModule(),
	m_nAsicX(0),
	m_nAsicY(0),
	m_nStartLayerShift(0),
	m_pElectronicsMapping(NULL),
	m_nActiveLayers(0),
	m_moduleLogStr("[AsicAnalysisModule]")
{
}

//-------------------------------------------------------------------------------------------------

AsicAnalysisModule::~AsicAnalysisModule() 
{
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::userReadSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	m_nActiveLayers = 48;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"NActiveLayers", m_nActiveLayers));

	m_nStartLayerShift = 0;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"NStartLayerShift", m_nStartLayerShift));

	m_nAsicX = 12;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"NAsicX", m_nAsicX));

	m_nAsicY = 12;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"NAsicY", m_nAsicY));

	m_cellIDDecoderString = "M:3,S-1:3,I:9,J:9,K-1:6";
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellIDDecoderString", m_cellIDDecoderString));

	m_inputCollectionName = "SDHCAL_HIT";
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InputCollectionName", m_inputCollectionName));

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
	m_trackingSettings.chiSquareLimit = 100.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Tracking.ChiSquareLimit", m_trackingSettings.chiSquareLimit));

	m_trackingSettings.maxTransverseRatio = 0.05f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Tracking.MaxTransverseRatio", m_trackingSettings.maxTransverseRatio));

	m_trackingSettings.cosThetaLimit = 0.9f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Tracking.CosThetaLimit", m_trackingSettings.cosThetaLimit));

	m_trackingSettings.printDebug = false;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Tracking.PrintDebug", m_trackingSettings.printDebug));

	/*-----------------------------------------------------*/
	m_efficiencySettings.maxRadius = 25.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Efficiency.MaxRadius", m_efficiencySettings.maxRadius));

	m_efficiencySettings.semiDigitalReadout = true;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Efficiency.SDHCALReadout", m_efficiencySettings.semiDigitalReadout));

	m_efficiencySettings.printDebug = false;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Efficiency.PrintDebug", m_efficiencySettings.printDebug));

	/*-----------------------------------------------------*/
	m_interactionFinderSettings.minSize = 4;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InteractionFinder.MinSize", m_interactionFinderSettings.minSize));

	m_interactionFinderSettings.maxRadius = 50.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InteractionFinder.MaxRadius", m_interactionFinderSettings.maxRadius));

	m_interactionFinderSettings.maxDepth = 4;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InteractionFinder.MaxDepth", m_interactionFinderSettings.maxDepth));

	m_interactionFinderSettings.minNumberOfCluster = 3;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InteractionFinder.MinNumberOfCluster", m_interactionFinderSettings.minNumberOfCluster));

	/*-----------------------------------------------------*/
	m_layerSettings.edgeX_min = 0.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Layer.EdgeX_Min", m_layerSettings.edgeX_min));

	m_layerSettings.edgeX_max = 1000.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Layer.EdgeX_Max", m_layerSettings.edgeX_max));

	m_layerSettings.edgeY_min = 0.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Layer.EdgeY_Min", m_layerSettings.edgeY_min));

	m_layerSettings.edgeY_max = 1000.f;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"Layer.EdgeY_Max", m_layerSettings.edgeY_max));

	/*-----------------------------------------------------*/

	TiXmlElement *pElecMapElement = xmlHandle.FirstChild("electronicsMapping").Element();

	if( ! pElecMapElement )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't find xml element electronicsMapping !" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	std::string plugin;
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::getAttribute(pElecMapElement, "plugin", plugin));

	m_pElectronicsMapping = dqm4hep::DQMPluginManager::instance()->createPluginClass<dqm4hep::DQMElectronicsMapping>(plugin);

	if( ! m_pElectronicsMapping )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Couldn't find electronicsMapping plugin called : " << plugin );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, m_pElectronicsMapping->readSettings(dqm4hep::TiXmlHandle(pElecMapElement)));

	/*-----------------------------------------------------*/
	for(unsigned int l=0 ; l<m_nActiveLayers ; l++)
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"EfficiencyMap", l, m_layerElementMap[l].m_pEfficiencyMap));
		std::string meTitle = m_layerElementMap[l].m_pEfficiencyMap->getTitle();
		meTitle += " layer " + std::to_string(l);
		m_layerElementMap[l].m_pEfficiencyMap->setTitle(meTitle.c_str());

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"Efficiency2Map", l, m_layerElementMap[l].m_pEfficiency2Map));
		meTitle = m_layerElementMap[l].m_pEfficiency2Map->getTitle();
		meTitle += " layer " + std::to_string(l);
		m_layerElementMap[l].m_pEfficiency2Map->setTitle(meTitle.c_str());

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"Efficiency3Map", l, m_layerElementMap[l].m_pEfficiency3Map));
		meTitle = m_layerElementMap[l].m_pEfficiency3Map->getTitle();
		meTitle += " layer " + std::to_string(l);
		m_layerElementMap[l].m_pEfficiency3Map->setTitle(meTitle.c_str());

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"MultiplicityMap", l, m_layerElementMap[l].m_pMultiplicityMap));
		meTitle = m_layerElementMap[l].m_pMultiplicityMap->getTitle();
		meTitle += " layer " + std::to_string(l);
		m_layerElementMap[l].m_pMultiplicityMap->setTitle(meTitle.c_str());
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerEfficiency", m_pLayerEfficiency));

	m_pLayerEfficiency->get<TGraph>()->SetMarkerSize(1);
	m_pLayerEfficiency->get<TGraph>()->SetMarkerColor(kBlack);
	m_pLayerEfficiency->get<TGraph>()->SetMarkerStyle(23);

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerEfficiency2", m_pLayerEfficiency2));

	m_pLayerEfficiency2->get<TGraph>()->SetMarkerSize(1);
	m_pLayerEfficiency2->get<TGraph>()->SetMarkerColor(kBlack);
	m_pLayerEfficiency2->get<TGraph>()->SetMarkerStyle(23);

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerEfficiency3", m_pLayerEfficiency3));

	m_pLayerEfficiency3->get<TGraph>()->SetMarkerSize(1);
	m_pLayerEfficiency3->get<TGraph>()->SetMarkerColor(kBlack);
	m_pLayerEfficiency3->get<TGraph>()->SetMarkerStyle(23);

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerMultiplicity", m_pLayerMultiplicity));

	m_pLayerMultiplicity->get<TGraph>()->SetMarkerSize(1);
	m_pLayerMultiplicity->get<TGraph>()->SetMarkerColor(kBlack);
	m_pLayerMultiplicity->get<TGraph>()->SetMarkerStyle(23);

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicEfficiency", m_pAsicEfficiency));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicEfficiency2", m_pAsicEfficiency2));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicEfficiency3", m_pAsicEfficiency3));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicMultiplicity", m_pAsicMultiplicity));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedEfficiencyMap", m_pStackedEfficiencyMap));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedEfficiency2Map", m_pStackedEfficiency2Map));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedEfficiency3Map", m_pStackedEfficiency3Map));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedMultiplicityMap", m_pStackedMultiplicityMap));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalEfficiency", m_pGlobalEfficiency));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalEfficiency2", m_pGlobalEfficiency2));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalEfficiency3", m_pGlobalEfficiency3));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalMultiplicity", m_pGlobalMultiplicity));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NTracksPerAsic", m_pNTracksPerAsic));


	DQMModuleApi::cd(this);
	DQMModuleApi::ls(this, true); // true for recursive

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::userInitModule()
{
	// initialize algorithms
	m_clusteringAlgorithm.SetClusterParameterSetting(m_clusteringSettings);
	m_clusteringHelper.SetClusteringHelperParameterSetting(m_clusteringHelperSettings);
	m_trackingAlgorithm.SetTrackingParameterSetting(m_trackingSettings);
	m_interactionFinderAlgorithm.SetInteractionFinderParameterSetting(m_interactionFinderSettings);
	m_efficiencyAlgorithm.SetEfficiencyParameterSetting(m_efficiencySettings);

	return STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::processEvent(EVENT::LCEvent *pLCEvent)
{
	LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Processing physics event no " << pLCEvent->getEventNumber() );

	// check for event rejection : noise ?
	bool rejectEvent = this->shouldRejectEvent(pLCEvent);

	if( rejectEvent )
	  {
	    LOG4CXX_INFO( dqm4hep::dqmMainLogger , m_moduleLogStr << " - rejecting noise " );
	    return dqm4hep::STATUS_CODE_SUCCESS;
	  }
	
	// content management
	caloobject::CaloHitMap caloHitMap;
	std::vector< caloobject::CaloHit *> hits;
	std::vector< caloobject::CaloCluster *> clusters;
	std::vector< caloobject::CaloTrack *>   tracks;

	CLHEP::Hep3Vector globalHitShift(0, 0, 0);

	try
	{
		EVENT::LCCollection *pCalorimeterHitCollection = pLCEvent->getCollection(m_inputCollectionName);

		if(NULL == pCalorimeterHitCollection)
		  {
		    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - NULL caloHitCllection " );
		    return dqm4hep::STATUS_CODE_SUCCESS;
		  }
		UTIL::CellIDDecoder<EVENT::CalorimeterHit> cellIDDecoder(m_cellIDDecoderString);

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Creating wrapper hits");

		// loop over hits in this event
		for(unsigned int h=0 ; h<pCalorimeterHitCollection->getNumberOfElements() ; h++)
		{
			EVENT::CalorimeterHit *pCaloHit = dynamic_cast<EVENT::CalorimeterHit*>(pCalorimeterHitCollection->getElementAt(h));

			if(NULL == pCaloHit)
			  {
			    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - rejecting nullCaoHit " ); 
			    continue;
			  }
			int cellID[3];
			cellID[0] = cellIDDecoder(pCaloHit)["I"];
			cellID[1] = cellIDDecoder(pCaloHit)["J"];
			cellID[2] = cellIDDecoder(pCaloHit)["K-1"] + m_nStartLayerShift;

			if( cellID[2] >= m_nActiveLayers )
			  {
			    LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - layer > max numbr of layers" );
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

			std::cout << "Input hit : " << pWrapperHit->getPosition() << std::endl;
		}

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Creating intra layer clusters");

		for(caloobject::CaloHitMap::iterator iter = caloHitMap.begin(), endIter = caloHitMap.end() ;
				iter != endIter ; ++iter)
		  m_clusteringAlgorithm.Run(iter->second, clusters);

		std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Filter non - isolated clusters");

		caloobject::CaloClusterList trackingClusters;

		std::cout << "before isolation : clusters.size() = " << clusters.size() << std::endl;
		for(std::vector<caloobject::CaloCluster*>::iterator iter = clusters.begin(), endIter = clusters.end() ;
				endIter != iter ; ++iter)

			if( ! m_clusteringHelper.IsIsolatedCluster(*iter, clusters) )
			  {
			    // LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Found non isolated hit " );
			    std::cout << "Position: " << (*iter)->getPosition() << std::endl;
			    std::cout << "Layer: " << (*iter)->getLayerID() << std::endl;
				trackingClusters.push_back(*iter);
			  }
			else
			  {
			    std::cout << "Isol Position: " << (*iter)->getPosition() << std::endl;
			    std::cout << "Isol Layer: " << (*iter)->getLayerID() << std::endl;
			  }
		
		std::cout << "After isolation : trackingCluster.size() = " << trackingClusters.size() << std::endl;

		caloobject::CaloTrack *pTrack = NULL;

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Run tracking algorithm");

		m_trackingAlgorithm.Run(trackingClusters, pTrack);

		// stop processing if no reconstructed track
		if( NULL == pTrack )
		{
		  LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - No track found !");
		  throw dqm4hep::StatusCodeException(dqm4hep::STATUS_CODE_SUCCESS);
		}

		tracks.push_back(pTrack);

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Run interaction finder");

		bool isInteraction = m_interactionFinderAlgorithm.Run(trackingClusters, pTrack->getTrackParameters());

		// tracking on muons only
		// stop processing if event is an interaction
		if( isInteraction )
			throw dqm4hep::StatusCodeException(dqm4hep::STATUS_CODE_SUCCESS);

		// update analysis contents : asics and layers

		int trackBegin = (*clusters.begin())->getLayerID();
		int trackEnd = (*(clusters.rbegin()))->getLayerID();

		if(1 == trackBegin)
			trackBegin = 0;

		if(m_nActiveLayers-2 == trackEnd)
			trackEnd = m_nActiveLayers-1;

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Analyzing track, layer per layer");

		for(unsigned int l = trackBegin ; l<=trackEnd ; l++)
		{
			caloobject::CaloLayer *pLayer = this->getOrCreateLayer(l);

			// reset layer properties
			pLayer->Reset();

			// and re-evaluate efficiency of layer
			m_efficiencyAlgorithm.Run(pLayer, clusters);

			if( 0 == pLayer->getNTracks() )
				continue;

			dqm4hep::DQMCartesianVector position(
					m_efficiencyAlgorithm.getExpectedPosition().x(),
					m_efficiencyAlgorithm.getExpectedPosition().y(),
					m_efficiencyAlgorithm.getExpectedPosition().z());

			if( ! this->isValid( position ) )
				continue;

			dqm4hep::DQMElectronicsMapping::Electronics electronics;
			dqm4hep::DQMElectronicsMapping::Cell cell;

			if(dqm4hep::STATUS_CODE_SUCCESS != m_pElectronicsMapping->positionToCell(position, cell))
				continue;

			cell.m_layer = pLayer->getLayerID();

			if(dqm4hep::STATUS_CODE_SUCCESS != m_pElectronicsMapping->cellToElectronics(cell, electronics))
				continue;

			this->updateAsic(electronics.m_difId, electronics.m_asicId, pLayer);
		}

		LOG4CXX_DEBUG( dqm4hep::dqmMainLogger , m_moduleLogStr << " - End of event, clearing content");

		this->clearEventContents(hits, clusters, tracks);
	}
	catch(EVENT::DataNotAvailableException &exception)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Caught EVENT::DataNotAvailableException : " << exception.what() << ". Skipping event ..." );
		this->clearEventContents(hits, clusters, tracks);
		return dqm4hep::STATUS_CODE_SUCCESS;
	}
	catch(dqm4hep::StatusCodeException &exception)
	{
		this->clearEventContents(hits, clusters, tracks);

		if(dqm4hep::STATUS_CODE_SUCCESS == exception.getStatusCode() )
			return dqm4hep::STATUS_CODE_SUCCESS;

		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Caught StatusCodeException : " << exception.toString() << ". Skipping event ..." );
		return dqm4hep::STATUS_CODE_SUCCESS;
	}
	catch(...)
	{
		this->clearEventContents(hits, clusters, tracks);
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - Caught unknown exception !" );
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::startOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endOfCycle()
{
	std::map<unsigned int, LayerInfo> layerInfoMap;

	this->resetElements();

	for(AsicMap::iterator iter = m_asicMap.begin(), endIter = m_asicMap.end() ;
			endIter != iter ; ++iter)
	{
		if( ! (iter->second->m_nTracks > 0) )
			continue;

		const unsigned int nTracks = iter->second->m_nTracks;
		const unsigned int layerID = iter->second->m_layerId;

		const float efficiency1 = nTracks ? iter->second->m_efficiency/nTracks : 0.f;
		const float efficiency2 = nTracks ? iter->second->m_efficiency2/nTracks : 0.f;
		const float efficiency3 = nTracks ? iter->second->m_efficiency3/nTracks : 0.f;

		const float x = iter->second->m_x;
		const float y = iter->second->m_y;

		const bool isEfficient = efficiency1 > 0.f;

		// const float efficiencyError1 = std::sqrt(efficiency1 * (1 - efficiency1) / nTracks);
		// const float efficiencyError2 = std::sqrt(efficiency2 * (1 - efficiency2) / nTracks);
		// const float efficiencyError3 = std::sqrt(efficiency3 * (1 - efficiency3) / nTracks);

		std::map<unsigned int, LayerInfo>::iterator infoIter = layerInfoMap.find(layerID);

		if( layerInfoMap.end() == infoIter )
		{
			infoIter = layerInfoMap.insert(std::map<unsigned int, LayerInfo>::value_type(layerID, LayerInfo())).first;
			infoIter->second.m_efficiency = 0.f;
			infoIter->second.m_efficiency2 = 0.f;
			infoIter->second.m_efficiency3 = 0.f;
			infoIter->second.m_multiplicity = 0.f;
			infoIter->second.m_count = 0;
			infoIter->second.m_efficientCount = 0;
		}

		infoIter->second.m_efficiency += efficiency1;
		infoIter->second.m_efficiency2 += efficiency2;
		infoIter->second.m_efficiency3 += efficiency3;
		infoIter->second.m_count ++;

		std::map<unsigned int, LayerElements>::iterator layerIter = m_layerElementMap.find(layerID);

		// global fill
		m_pNTracksPerAsic->get<TH1>()->Fill( nTracks );

		// fill efficiency
		m_pAsicEfficiency->get<TH1F>()->Fill( efficiency1 );
		m_pAsicEfficiency2->get<TH1F>()->Fill( efficiency2 );
		m_pAsicEfficiency3->get<TH1F>()->Fill( efficiency3 );

//		std::cout << "Filling stacked eff at (" << x << " , " << y << ")" << std::endl;
		m_pStackedEfficiencyMap->get<TH2F>()->Fill(x, y, efficiency1 / m_nActiveLayers);

		if( m_layerElementMap.end() != layerIter )
		{
			layerIter->second.m_pEfficiencyMap->get<TH2>()->Fill(x, y, efficiency1);
			layerIter->second.m_pEfficiency2Map->get<TH2>()->Fill(x, y, efficiency2);
			layerIter->second.m_pEfficiency3Map->get<TH2>()->Fill(x, y, efficiency3);
		}

		// fill multiplicity
		if( isEfficient )
		{
			const float multiplicity = (efficiency1 > std::numeric_limits<float>::epsilon()) ? iter->second->m_multiplicity / iter->second->m_efficiency : 0.f;
			m_pAsicMultiplicity->get<TH1F>()->Fill( multiplicity );
			m_pStackedMultiplicityMap->get<TH2F>()->Fill(x, y, multiplicity / m_nActiveLayers);

			if( m_layerElementMap.end() != layerIter )
				layerIter->second.m_pMultiplicityMap->get<TH2>()->Fill(x, y, multiplicity);

			infoIter->second.m_multiplicity += multiplicity;
			infoIter->second.m_efficientCount ++;
		}
	}

	LayerInfo globalInfo;
	globalInfo.m_efficiency = 0.f;
	globalInfo.m_efficiency2 = 0.f;
	globalInfo.m_efficiency3 = 0.f;
	globalInfo.m_multiplicity = 0.f;
	globalInfo.m_count = 0;
	globalInfo.m_efficientCount = 0;

	// fill efficiency and multiplicity per layer
	// compute the global efficiency and multiplicity for all layers
	for(std::map<unsigned int, LayerInfo>::iterator iter = layerInfoMap.begin(), endIter = layerInfoMap.end() ;
			endIter != iter ; ++iter)
	{
		const unsigned int layerID = iter->first;

		if(iter->second.m_count > 0)
		{
			const float layerEfficiency = ( iter->second.m_efficiency / iter->second.m_count );
			const float layerEfficiency2 = ( iter->second.m_efficiency2 / iter->second.m_count );
			const float layerEfficiency3 = ( iter->second.m_efficiency3 / iter->second.m_count );

			Int_t pointID = m_pLayerEfficiency->get<TGraph>()->GetN();
			m_pLayerEfficiency->get<TGraph>()->SetPoint(pointID, layerID , layerEfficiency * 100 );
			m_pLayerEfficiency2->get<TGraph>()->SetPoint(pointID, layerID , layerEfficiency2 * 100 );
			m_pLayerEfficiency3->get<TGraph>()->SetPoint(pointID, layerID , layerEfficiency3 * 100 );

			globalInfo.m_efficiency += layerEfficiency;
			globalInfo.m_efficiency2 += layerEfficiency2;
			globalInfo.m_efficiency3 += layerEfficiency3;
			globalInfo.m_count ++;
		}

		if(iter->second.m_efficientCount > 0)
		{
			const float layerMultiplicity = ( iter->second.m_multiplicity / iter->second.m_efficientCount );

			Int_t pointID = m_pLayerMultiplicity->get<TGraph>()->GetN();
			m_pLayerMultiplicity->get<TGraph>()->SetPoint( pointID , layerID , layerMultiplicity );

			globalInfo.m_multiplicity += layerMultiplicity;
			globalInfo.m_efficientCount ++;
		}
	}

	// set global detector efficiency and multiplicity
	m_pGlobalEfficiency->get< TScalarObject<float> >()->Clear();
	m_pGlobalEfficiency2->get< TScalarObject<float> >()->Clear();
	m_pGlobalEfficiency3->get< TScalarObject<float> >()->Clear();
	m_pGlobalMultiplicity->get< TScalarObject<float> >()->Clear();

	if( globalInfo.m_count > 0 )
	{
		m_pGlobalEfficiency->get< TScalarObject<float> >()->Set( ( globalInfo.m_efficiency / globalInfo.m_count) * 100 );
		m_pGlobalEfficiency2->get< TScalarObject<float> >()->Set( ( globalInfo.m_efficiency2 / globalInfo.m_count) * 100 );
		m_pGlobalEfficiency3->get< TScalarObject<float> >()->Set( ( globalInfo.m_efficiency3 / globalInfo.m_count) * 100 );
	}

	if( globalInfo.m_efficientCount > 0 )
		m_pGlobalMultiplicity->get< TScalarObject<float> >()->Set( globalInfo.m_multiplicity / globalInfo.m_efficientCount );

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::startOfRun(DQMRun *pRun)
{
	this->clearContents();
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endOfRun(DQMRun *pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endModule()
{
	this->clearContents();
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::clearEventContents(caloobject::CaloHitList &hits, caloobject::CaloClusterList &clusters, caloobject::CaloTrackList &tracks)
{
	for_each(hits.begin(), hits.end(), [] (caloobject::CaloHit *pCaloHit) { delete pCaloHit; });
	for_each(clusters.begin(), clusters.end(), [] (caloobject::CaloCluster *pCluster) { delete pCluster; });
	for_each(tracks.begin(), tracks.end(), [] (caloobject::CaloTrack *pTrack) { delete pTrack; });

	hits.clear();
	clusters.clear();
	tracks.clear();
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::clearContents()
{
	for(caloobject::CaloLayerMap::iterator iter = m_caloLayerMap.begin(), endIter = m_caloLayerMap.end() ;
			endIter != iter ; ++iter)
		delete iter->second;

	for(AsicMap::iterator iter = m_asicMap.begin(), endIter = m_asicMap.end() ;
			endIter != iter ; ++iter)
		delete iter->second;

	m_caloLayerMap.clear();
	m_asicMap.clear();
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::createAsicKey(unsigned int difId, unsigned int asicId, unsigned int &key)
{
	key = asicId + difId*1000;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::decodeAsicKey(unsigned int key, unsigned int &difId, unsigned int &asicId)
{
	difId = key / 1000;
	asicId = key%difId;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::updateAsic(unsigned int difId, unsigned int asicId, caloobject::CaloLayer *pLayer)
{
	unsigned int key(0);
	this->createAsicKey(difId, asicId, key);

	AsicMap::iterator iter = m_asicMap.find(key);

	if( m_asicMap.end() == iter )
	{
		Asic *pAsic = new Asic();

		// fill electronics ids
		pAsic->m_difId = difId;
		pAsic->m_asicId = asicId;
		pAsic->m_layerId = pLayer->getLayerID();

		// find a generic asic position
		dqm4hep::DQMElectronicsMapping::Electronics electronics;
		electronics.m_difId = pAsic->m_difId;
		electronics.m_asicId = pAsic->m_asicId;
		electronics.m_channelId = 1;

		dqm4hep::DQMElectronicsMapping::Cell cell;

		if(dqm4hep::STATUS_CODE_SUCCESS != m_pElectronicsMapping->electronicstoCell(electronics, cell))
			return;

		pAsic->m_x = cell.m_iCell;
		pAsic->m_y = cell.m_jCell;

		// insert and set iterator to the inserted one
		iter = m_asicMap.insert( AsicMap::value_type(key, pAsic) ).first;
	}

	iter->second->m_nTracks++;
	iter->second->m_efficiency += pLayer->getEfficiency();
	iter->second->m_multiplicity += pLayer->getMultiplicity();

	if( pLayer->getEfficiencyEnergy() == 3 )
	{
		iter->second->m_efficiency3 += pLayer->getEfficiency();
		iter->second->m_efficiency2 += pLayer->getEfficiency();
	}
	else if( pLayer->getEfficiencyEnergy() == 2 )
	{
		iter->second->m_efficiency2 += pLayer->getEfficiency();
	}
}

//-------------------------------------------------------------------------------------------------

caloobject::CaloLayer *AsicAnalysisModule::getOrCreateLayer(unsigned int layerId)
{
	caloobject::CaloLayerMap::iterator iter = m_caloLayerMap.find(layerId);

	if( m_caloLayerMap.end() == iter )
	{
		caloobject::CaloLayer *pLayer = new caloobject::CaloLayer( layerId );

		pLayer->setLayerParameterSetting(m_layerSettings);

		iter = m_caloLayerMap.insert( caloobject::CaloLayerMap::value_type( layerId , pLayer ) ).first;
	}

	return iter->second;
}

//-------------------------------------------------------------------------------------------------

bool AsicAnalysisModule::isValid(const dqm4hep::DQMCartesianVector &vector)
{
	const float x = vector.getX();
	const float y = vector.getY();
	const float z = vector.getZ();

	if( std::isinf(x) || std::isnan(x)
	 || std::isinf(y) || std::isnan(y)
	 || std::isinf(z) || std::isnan(z) )
		return false;

	return true;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::resetElements()
{
	for(std::map<unsigned int, LayerElements>::iterator iter = m_layerElementMap.begin(), endIter = m_layerElementMap.end() ;
			endIter != iter ; ++iter)
	{
		iter->second.m_pEfficiencyMap->reset();
		iter->second.m_pEfficiency2Map->reset();
		iter->second.m_pEfficiency3Map->reset();
		iter->second.m_pMultiplicityMap->reset();
	}

	m_pLayerEfficiency->get<TGraph>()->Set(0);
	m_pLayerEfficiency2->get<TGraph>()->Set(0);
	m_pLayerEfficiency3->get<TGraph>()->Set(0);
	m_pLayerMultiplicity->get<TGraph>()->Set(0);

	m_pAsicEfficiency->reset();
	m_pAsicEfficiency2->reset();
	m_pAsicEfficiency3->reset();
	m_pAsicMultiplicity->reset();

	m_pStackedEfficiencyMap->reset();
	m_pStackedEfficiency2Map->reset();
	m_pStackedEfficiency3Map->reset();
	m_pStackedMultiplicityMap->reset();

	m_pGlobalEfficiency->reset();
	m_pGlobalEfficiency2->reset();
	m_pGlobalEfficiency3->reset();
	m_pGlobalMultiplicity->reset();

	m_pNTracksPerAsic->reset();
}


bool AsicAnalysisModule::shouldRejectEvent(EVENT::LCEvent *pLCEvent)
{
	return false;
}

}

