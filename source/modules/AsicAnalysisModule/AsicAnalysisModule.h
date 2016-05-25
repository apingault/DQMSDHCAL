  /// \file AsicAnalysisModule.h
/*
 *
 * AsicAnalysisModule.h header template automatically generated by a class generator
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


#ifndef DQMSDHCAL_ASICANALYSISMODULE_H
#define DQMSDHCAL_ASICANALYSISMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMElectronicsMapping.h"

// -- lcio headers
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"

// -- calo software headers
#include "CaloObject/Asic.h"
#include "Algorithm/Cluster.h"
#include "CaloObject/CaloHit.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Efficiency.h"
#include "Algorithm/AsicKeyFinder.h"

// -- std headers
#include <string>
#include <cstring>
#include <vector>

// -- dqmsdhcal headers
#include "DQMTriventModule.h"

namespace caloobject
{
	typedef std::map<unsigned int, std::vector<CaloHit *> > CaloHitMap;
	typedef std::map<unsigned int, CaloLayer *> CaloLayerMap;
	typedef std::vector<CaloHit *> CaloHitList;
	typedef std::vector<CaloCluster *> CaloClusterList;
	typedef std::vector<CaloTrack *> CaloTrackList;
}

namespace dqm_sdhcal
{

/** AsicAnalysisModule class
 */ 
class AsicAnalysisModule : public DQMTriventModule
{
public:
	/** Asic helper structure to compute
	 *  efficiency and multiplicity
	 */
	struct Asic
	{
		Asic() :
			m_difId(0),
			m_asicId(0),
			m_layerId(0),
			m_nTracks(0),
			m_efficiency(0.f),
			m_efficiency2(0.f),
			m_efficiency3(0.f),
			m_multiplicity(0.f),
			m_x(0.f),
			m_y(0.f)
		{
			/* nop */
		}

		unsigned int      m_difId;
		unsigned int      m_asicId;
		unsigned int      m_layerId;
		unsigned int      m_nTracks;
		float             m_efficiency;
		float             m_efficiency2;
		float             m_efficiency3;
		float             m_multiplicity;
		float             m_x;
		float             m_y;
	};

	typedef std::map<unsigned int, Asic *>      AsicMap;

	/** LayerInfo helper structure to compute
	 *  efficiency and multiplicity
	 */
	struct LayerInfo
	{
		float           m_efficiency;
		float           m_efficiency2;
		float           m_efficiency3;
		float           m_multiplicity;
		unsigned int    m_count;
		unsigned int    m_efficientCount;
	};

public:
	/** Constructor
	 */
	AsicAnalysisModule();

	/** Destructor
	 */
	virtual ~AsicAnalysisModule();

	dqm4hep::StatusCode userInitModule();
	dqm4hep::StatusCode userReadSettings(const dqm4hep::TiXmlHandle xmlHandle);
	dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

	dqm4hep::StatusCode startOfCycle();
	dqm4hep::StatusCode endOfCycle();
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endModule();

private:

	//
	// UTILITY FUNCTIONS
	//

	/**
	 */
	void createAsicKey(unsigned int difId, unsigned int asicId, unsigned int &key);

	/**
	 */
	void decodeAsicKey(unsigned int key, unsigned int &difId, unsigned int &asicId);

	/**
	 */
	void updateAsic(unsigned int difId, unsigned int asicId, caloobject::CaloLayer *pLayer);

	/**
	 */
	caloobject::CaloLayer *getOrCreateLayer(unsigned int layerId);

	/**
	 */
	bool isValid(const dqm4hep::DQMCartesianVector &vector);


	//
	// Clear content function
	//

	/** Clear the contents related to one event
	 */
	void clearEventContents(caloobject::CaloHitList &hits,
			caloobject::CaloClusterList &clusters, caloobject::CaloTrackList &tracks);

	/** Clear the contents related
	 */
	void clearContents();

	/** Create the asic and layer contents needed
	 *  to extract the asic efficiency and multiplicity
	 */
	void createAsicAndLayerContents();

	/** Reset the monitor element tuned on cycle
	 *  Called at end of cycle before filling again
	 *  efficiencies and multiplicities
	 */
	void resetElements();

	//
	// ADDITIONAL ANALYSIS FUNCTIONS
	//

	/** Whether the vent has to be reject before any further processing
	 */
	bool shouldRejectEvent(EVENT::LCEvent *pLCEvent);

private:

	std::string 														m_moduleLogStr;

	// algorithm contents
	caloobject::CaloLayerMap                 m_caloLayerMap;
	AsicMap                                  m_asicMap;

	// algorithms
	algorithm::Cluster                       m_clusteringAlgorithm;
	algorithm::ClusteringHelper              m_clusteringHelper;
	algorithm::Tracking                      m_trackingAlgorithm;
	algorithm::InteractionFinder             m_interactionFinderAlgorithm;
	algorithm::Efficiency                    m_efficiencyAlgorithm;

	// algorithm parameters
	algorithm::clusterParameterSetting           m_clusteringSettings;
	algorithm::ClusteringHelperParameterSetting  m_clusteringHelperSettings;
	algorithm::TrackingParameterSetting          m_trackingSettings;
	algorithm::InteractionFinderParameterSetting m_interactionFinderSettings;
	algorithm::EfficiencyParameterSetting        m_efficiencySettings;
	caloobject::LayerParameterSetting            m_layerSettings;

	// module parameters
	std::string                              m_inputCollectionName;
	std::string                              m_cellIDDecoderString;
	dqm4hep::IntVector                       m_asicTable;
	dqm4hep::IntVector                       m_difList;
	unsigned int                             m_nAsicX;
	unsigned int                             m_nAsicY;
	int                                      m_nStartLayerShift;
	unsigned int                             m_nActiveLayers;
	dqm4hep::DQMElectronicsMapping          *m_pElectronicsMapping;

private:
	// monitor elements

	struct LayerElements
	{
		dqm4hep::DQMMonitorElementPtr           m_pEfficiencyMap;
		dqm4hep::DQMMonitorElementPtr           m_pEfficiency2Map;
		dqm4hep::DQMMonitorElementPtr           m_pEfficiency3Map;
		dqm4hep::DQMMonitorElementPtr           m_pMultiplicityMap;
	};

	dqm4hep::DQMMonitorElementPtr               m_pLayerEfficiency;
	dqm4hep::DQMMonitorElementPtr               m_pLayerEfficiency2;
	dqm4hep::DQMMonitorElementPtr               m_pLayerEfficiency3;
	dqm4hep::DQMMonitorElementPtr               m_pLayerMultiplicity;
	dqm4hep::DQMMonitorElementPtr               m_pAsicEfficiency;
	dqm4hep::DQMMonitorElementPtr               m_pAsicEfficiency2;
	dqm4hep::DQMMonitorElementPtr               m_pAsicEfficiency3;
	dqm4hep::DQMMonitorElementPtr               m_pAsicMultiplicity;
	dqm4hep::DQMMonitorElementPtr               m_pStackedEfficiencyMap;
	dqm4hep::DQMMonitorElementPtr               m_pStackedEfficiency2Map;
	dqm4hep::DQMMonitorElementPtr               m_pStackedEfficiency3Map;
	dqm4hep::DQMMonitorElementPtr               m_pStackedMultiplicityMap;
	dqm4hep::DQMMonitorElementPtr               m_pGlobalEfficiency;
	dqm4hep::DQMMonitorElementPtr               m_pGlobalEfficiency2;
	dqm4hep::DQMMonitorElementPtr               m_pGlobalEfficiency3;
	dqm4hep::DQMMonitorElementPtr               m_pGlobalMultiplicity;
	dqm4hep::DQMMonitorElementPtr               m_pNTracksPerAsic;

	std::map<unsigned int, LayerElements>       m_layerElementMap;
}; 

} 

#endif  //  DQMSDHCAL_ASICANALYSISMODULE_H
