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
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMEvent.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMPlugin.h"

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

#include "streamlog/streamlog.h"

using namespace dqm4hep;

namespace dqm_sdhcal
{

// plugin declaration
DQM_PLUGIN_DECL( AsicAnalysisModule, "AsicAnalysisModule" )

AsicAnalysisModule::AsicAnalysisModule() :
	DQMTriventModule()
{
}

//-------------------------------------------------------------------------------------------------

AsicAnalysisModule::~AsicAnalysisModule() 
{
}
//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::userInitModule()
{
	return STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::userReadSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	m_nActiveLayers = 48;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"NActiveLayers", m_nActiveLayers));

	for(unsigned int l=0 ; l<m_nActiveLayers ; l++)
	{
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"EfficiencyMap", l, m_layerElementMap[l].m_pEfficiencyMap));

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
				"MultiplicityMap", l, m_layerElementMap[l].m_pMultiplicityMap));
	}

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerEfficiency", m_pLayerEfficiency));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"LayerMultiplicity", m_pLayerMultiplicity));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicEfficiency", m_pAsicEfficiency));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"AsicMultiplicity", m_pAsicMultiplicity));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedEfficiencyMap", m_pStackedEfficiencyMap));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"StackedMultiplicityMap", m_pStackedMultiplicityMap));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalEfficiency", m_pGlobalEfficiency));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"GlobalMultiplicity", m_pGlobalMultiplicity));

	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::bookMonitorElement(this, xmlHandle,
			"NTracksPerAsic", m_pNTracksPerAsic));

	m_inputCollectionName = "SDHCAL_HIT";
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"InputCollectionName", m_inputCollectionName));

	return dqm4hep::STATUS_CODE_SUCCESS;
}

dqm4hep::StatusCode AsicAnalysisModule::processNoisyEvent(EVENT::LCEvent *pLCEvent)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

dqm4hep::StatusCode AsicAnalysisModule::processPhysicalEvent(EVENT::LCEvent *pLCEvent)
{
	LOG4CXX_INFO( dqm4hep::dqmMainLogger , "Processing physics event no " << pLCEvent->getEventNumber() );

	try
	{
		EVENT::LCCollection *pCalorimeterHitCollection = pLCEvent->getCollection(m_inputCollectionName);

		if(NULL == pCalorimeterHitCollection)
			return dqm4hep::STATUS_CODE_SUCCESS;

		UTIL::CellIDDecoder<EVENT::CalorimeterHit> cellIDDecoder(pCalorimeterHitCollection);

		// loop over hits in this event
		for(unsigned int h=0 ; h<pCalorimeterHitCollection->getNumberOfElements() ; h++)
		{
			EVENT::CalorimeterHit *pCaloHit = dynamic_cast<EVENT::CalorimeterHit*>(pCalorimeterHitCollection->getElementAt(h));

			if(NULL == pCaloHit)
				continue;

			m_calorimeterHitCollection.push_back(pCaloHit);
		}

		this->doTrackStudy();
		m_calorimeterHitCollection.clear();
	}
	catch(EVENT::DataNotAvailableException &exception)
	{
		streamlog_out(ERROR) << "Caught EVENT::DataNotAvailableException : " << exception.what() << std::endl;
		streamlog_out(ERROR) << "Skipping event" << std::endl;
		this->clearContents();
		return STATUS_CODE_SUCCESS;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::doTrackStudy()
{
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
	clusters.clear();

	std::vector<EVENT::CalorimeterHit*> _temp;
	int ID = 0;
	int nclusters = 0;
	Cluster* cluster = NULL;

	for(std::vector<EVENT::CalorimeterHit*>::iterator it=m_calorimeterHitCollection.begin(); it!=m_calorimeterHitCollection.end(); ++it)
	{
		if(std::find(_temp.begin(),_temp.end(), (*it) )!=_temp.end())
			continue;

		cluster = new Cluster(IDdecoder(*it)["K-1"]);
		cluster->AddHits(*it);
		nclusters++;
		ID++;
		_temp.push_back(*it);

		cluster->BuildCluster(_temp, m_calorimeterHitCollection, (*it));
		cluster->buildClusterPosition();
		cluster->setClusterID(ID);

		clusters.push_back(cluster);
	}

	for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
		(*it)->IsolatedCluster(clusters);

	for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
	{
		if( (*it)->isIsolated() )
		{
			streamlog_out( DEBUG ) << "cluster at " << (*it)->getClusterPosition().x() << " " << (*it)->getClusterPosition().y() << " " << (*it)->getClusterPosition().z()
					<< " is isolated and rejected" << std::endl;
			delete *it;
			clusters.erase(it);
			it--;
		}
	}

	std::sort(clusters.begin(), clusters.end(), ClusterClassFunction::sortDigitalClusterByLayer);

	if(clusters.size() > 5)
		if(TrackSelection(clusters))
			LayerProperties(clusters);

	for(std::vector<Cluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
		delete *it;
}

//-------------------------------------------------------------------------------------------------

bool AsicAnalysisModule::TrackSelection(std::vector<Cluster*> &clVec)
{
	TrackingAlgo* aTrackingAlgo = new TrackingAlgo();

	aTrackingAlgo->Init(clVec);
	aTrackingAlgo->DoTracking();
	bool success = aTrackingAlgo->TrackFinderSuccess();

	delete aTrackingAlgo;

	return success;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::LayerProperties(std::vector<Cluster*> &clVec)
{
	int trackBegin= (*clVec.begin())->getLayerID();
	int trackEnd=(*(clVec.end()-1))->getLayerID();

	if(trackBegin == 1)
		trackBegin = 0;

	if(trackEnd == 46)
		trackEnd = 47;

	for(int K=trackBegin; K<=trackEnd; K++)
	{
		Layer* aLayer = new Layer(K);
		aLayer->Init(clVec);
		aLayer->ComputeLayerProperties();

		int asicKey = findAsicKey(K,aLayer->getxExpected(),aLayer->getyExpected());

		if(asicKey<0)
		{
			delete aLayer;
			continue;
		}

		if(asicMap.find(asicKey)==asicMap.end())
		{
			Asic* asic = new Asic(asicKey);
			asicMap[asicKey] = asic;
		}

		if(aLayer->getLayerTag()==fUnefficientLayer)
			asicMap[asicKey]->Update(0);

		if(aLayer->getLayerTag()==fEfficientLayer)
			asicMap[asicKey]->Update(aLayer->getMultiplicity());

		delete aLayer;
	}
}

//-------------------------------------------------------------------------------------------------

int AsicAnalysisModule::findAsicKey(int layer,float x, float y)
{
	float I=round( x/10.408 );
	float J=round( y/10.408 );

	if(I>96||I<0||J>96||J<0)
		return -1;

	int jnum=(J-1)/8;
	int inum=(I-1)/8;
	int num=jnum*12+inum;

	return layer*1000+num;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::startOfCycle()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endOfCycle()
{
	DQMMonitorElement *pMonitorElement = NULL;

	float totalEfficiency = 0.f;
	unsigned int nEfficientAsics = 0;

	float totalMultiplicity = 0.f;
	unsigned int nMultiplicityAsics = 0;

	this->resetElements();

	for(unsigned int i=0; i<m_nActiveLayers; i++)
	{
		std::map<unsigned int, LayerElements>::iterator iter = m_layerElementMap.find(i);

		if(iter == m_layerElementMap.end())
			continue;

		float layerEfficiency = 0.f;
		unsigned int nLayerEfficientAsics = 0;

		float layerMultiplicity = 0.f;
		unsigned int nLayerMultiplicityAsics = 0;

		for( std::map<int,Asic*>::iterator it = asicMap.begin() ; it != asicMap.end() ; it++ )
		{
			if( it->first/1000 != i )
				continue;

			if( it->second->getAsicCounter() > 0 )
			{
				nLayerEfficientAsics++;
				nEfficientAsics++;

				float asicEfficiency = it->second->getAsicEfficiency()*1.f/it->second->getAsicCounter();
				std::cout << "asicEfficiency : " << asicEfficiency << std::endl;
				layerEfficiency += asicEfficiency;
				totalEfficiency += asicEfficiency;

				// fill elements
				iter->second.m_pEfficiencyMap->get<TH2F>()->Fill( it->second->getAsicPosition()[0], it->second->getAsicPosition()[1], asicEfficiency );
				m_pAsicEfficiency->get<TH1F>()->Fill( asicEfficiency );
				m_pStackedEfficiencyMap->get<TH2F>()->Fill( it->second->getAsicPosition()[0], it->second->getAsicPosition()[1], asicEfficiency / m_nActiveLayers );
			}

			if( it->second->getAsicEfficiency() > 0 )
			{
				nMultiplicityAsics++;
				nLayerMultiplicityAsics++;

				float asicMultiplicity = it->second->getAsicMultiplicity()*1.f/it->second->getAsicEfficiency();
				layerMultiplicity += asicMultiplicity;
				totalMultiplicity += asicMultiplicity;

				// fill elements with multiplicity
				iter->second.m_pMultiplicityMap->get<TH2F>()->Fill( it->second->getAsicPosition()[0], it->second->getAsicPosition()[1], asicMultiplicity );
				m_pAsicMultiplicity->get<TH1F>()->Fill( asicMultiplicity );
				m_pStackedEfficiencyMap->get<TH2F>()->Fill( it->second->getAsicPosition()[0], it->second->getAsicPosition()[1] , asicMultiplicity / m_nActiveLayers );
			}

			m_pNTracksPerAsic->get<TH1F>()->Fill( it->second->getAsicCounter() );
		}

		if( nLayerEfficientAsics != 0 )
		{
			std::cout << "m_pLayerEfficiency = " << m_pLayerEfficiency << std::endl;
			std::cout << "m_pLayerEfficiency->get<TH1I>() = " << m_pLayerEfficiency->get<TH1I>() << std::endl;
			m_pLayerEfficiency->get<TH1I>()->Fill( i , layerEfficiency/nLayerEfficientAsics );
		}

		if( nLayerMultiplicityAsics != 0 )
			m_pLayerMultiplicity->get<TH1I>()->Fill( i , layerMultiplicity/nLayerMultiplicityAsics );
	}

	if( nEfficientAsics != 0 )
		m_pGlobalEfficiency->get< TScalarObject<float> >()->Set( totalEfficiency / nEfficientAsics );

	if( nMultiplicityAsics != 0 )
		m_pGlobalMultiplicity->get< TScalarObject<float> >()->Set( totalMultiplicity / nMultiplicityAsics );

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::startOfRun(DQMRun *pRun)
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endOfRun(DQMRun *pRun)
{
	this->clearContents();
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode AsicAnalysisModule::endModule()
{
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::clearContents()
{
	m_calorimeterHitCollection.clear();

	for(std::vector<Cluster*>::iterator iter = clusters.begin(), endIter = clusters.end() ;
			endIter != iter ; ++iter)
	{
		delete *iter;
	}

	clusters.clear();

	for(std::map<int,Asic*>::iterator iter = asicMap.begin(), endIter = asicMap.end() ;
			endIter != iter ; ++iter)
	{
		delete iter->second;
	}

	asicMap.clear();
}

//-------------------------------------------------------------------------------------------------

void AsicAnalysisModule::resetElements()
{
	for(std::map<unsigned int, LayerElements>::iterator iter = m_layerElementMap.begin(), endIter = m_layerElementMap.end() ;
			endIter != iter ; ++iter)
	{
		iter->second.m_pEfficiencyMap->reset();
		iter->second.m_pMultiplicityMap->reset();
	}

	m_pLayerEfficiency->reset();
	m_pLayerMultiplicity->reset();
	m_pAsicEfficiency->reset();
	m_pAsicMultiplicity->reset();
	m_pStackedEfficiencyMap->reset();
	m_pStackedMultiplicityMap->reset();
	m_pGlobalEfficiency->reset();
	m_pGlobalMultiplicity->reset();
	m_pNTracksPerAsic->reset();
}


}
