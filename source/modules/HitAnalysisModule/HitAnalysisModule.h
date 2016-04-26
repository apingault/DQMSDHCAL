/// \file HitAnalysisModule.h
/*
 *
 * HitAnalysisModule.h header template automatically generated by a class generator
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


#ifndef DQMSDHCAL_HITANALYSISMODULE_H
#define DQMSDHCAL_HITANALYSISMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- lcio headers
#include "lcio.h"
#include "EVENT/CalorimeterHit.h"

// -- calo software headers
#include "Algorithm/Cluster.h"
#include "CaloObject/CaloHit.h"
#include "Algorithm/ClusteringHelper.h"

// -- std headers
#include <string>
#include <cstring>
#include <vector>

// -- dqmsdhcal headers
#include "DQMTriventModule.h"

namespace caloobject
{
typedef std::map<unsigned int, std::vector<CaloHit *> > CaloHitMap;
typedef std::vector<CaloHit *> CaloHitList;
typedef std::vector<CaloCluster *> CaloClusterList;
}

namespace dqm_sdhcal
{

class EventClassifier;
class EventHelper;
/** HitAnalysisModule class
 */
class HitAnalysisModule : public DQMTriventModule
{
public:
	/** Constructor
	 */
	HitAnalysisModule();

	/** Destructor
	 */
	virtual ~HitAnalysisModule();

	dqm4hep::StatusCode userInitModule();
	dqm4hep::StatusCode userReadSettings(const dqm4hep::TiXmlHandle xmlHandle);
	dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

	dqm4hep::StatusCode startOfCycle();
	dqm4hep::StatusCode endOfCycle();
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endModule();

private:
	/** Clear the contents related to one event
	 */
	void clearEventContents(caloobject::CaloHitList &hits, caloobject::CaloClusterList &clusters);

	/**
	 *
	 */
dqm4hep::StatusCode fillRates();


	/** Reset the monitor element tuned on cycle
	 *  Called at end of cycle before filling again
	 *  efficiencies and multiplicities
	 */
	void resetElements();

private:
	// algorithms
	algorithm::Cluster                       m_clusteringAlgorithm;
	algorithm::ClusteringHelper              m_clusteringHelper;

	// algorithm parameters
	algorithm::clusterParameterSetting           m_clusteringSettings;
	algorithm::ClusteringHelperParameterSetting  m_clusteringHelperSettings;

	// module parameters
	std::string 														 m_moduleLogStr;
	std::string                              m_inputCollectionName;
	std::string                              m_cellIDDecoderString;
	unsigned int                             m_nActiveLayers;
	unsigned int                             m_firstLayerCut;
	unsigned int                             m_lastLayerCut;
	unsigned int 													   m_nMipMinimum;
	unsigned int 													   m_nMipInLayer;

  EventHelper															*m_pEventHelper;
  double 																	 m_eventIntegratedTime;
  double 																	 m_timeLastTrigger;
  double 																	 m_timeLastSpill;
  double 																	 m_spillIntegratedTime;
  double 																	 m_totalIntegratedTime;
	int 																		 m_nTrigger;
	int 																		 m_nSpill;

  EventClassifier 											  *m_pEventClassifier;
  int 																	   m_nParticleWithinRun;
  int 																	   m_nParticleWithinSpill;
  int 																	   m_nBeamMuonWithinRun;
  int 																	   m_nBeamMuonWithinSpill;
  int 																	   m_nChargedHadronsWithinRun;
  int 																	   m_nChargedHadronsWithinSpill;
  int 																	   m_nNeutralHadronsWithinRun;
  int 																	   m_nNeutralHadronsWithinSpill;
  int 																	   m_nPhotonsWithinRun;
  int 																	   m_nPhotonsWithinSpill;
  int 																	   m_nElectronsWithinRun;
  int 																	   m_nElectronsWithinSpill;
  int 																	   m_nOthersWithinRun;
  int 																	   m_nOthersWithinSpill;
	int 																	   m_nCosmicMuonsWithinRun;
	int 																	   m_nCosmicMuonsWithinSpill;
  int 																	   m_nUndefinedWithinRun;
  int 																	   m_nUndefinedWithinSpill;
  int 																	   m_nNoiseWithinRun;
  int 																	   m_nNoiseWithinSpill;

	// monitor elements
	//
	dqm4hep::DQMMonitorElementPtr			m_pInstantRate;
	dqm4hep::DQMMonitorElementPtr			m_pMeanRunRate;
	dqm4hep::DQMMonitorElementPtr			m_pNHit0;
	dqm4hep::DQMMonitorElementPtr			m_pNHit1;
	dqm4hep::DQMMonitorElementPtr			m_pNHit2;
	dqm4hep::DQMMonitorElementPtr			m_pNHit;

	dqm4hep::DQMMonitorElementPtr			m_pNHit0PerLayer;
	dqm4hep::DQMMonitorElementPtr			m_pNHit1PerLayer;
	dqm4hep::DQMMonitorElementPtr			m_pNHit2PerLayer;
	dqm4hep::DQMMonitorElementPtr			m_pNHitTotPerLayer;
	dqm4hep::DQMMonitorElementPtr			m_pRateVsClusterProfile;


	struct LayerElements
	{
		dqm4hep::DQMMonitorElementPtr			m_pNHit0Layer;
		dqm4hep::DQMMonitorElementPtr			m_pNHit1Layer;
		dqm4hep::DQMMonitorElementPtr			m_pNHit2Layer;
		dqm4hep::DQMMonitorElementPtr			m_pNHitTotLayer;
		dqm4hep::DQMMonitorElementPtr			m_pChamberHitsMap0;
		dqm4hep::DQMMonitorElementPtr			m_pChamberHitsMap1;
		dqm4hep::DQMMonitorElementPtr			m_pChamberHitsMap2;

	};

	std::map<unsigned int, LayerElements>       m_layerElementMap;
};

}

#endif  //  DQMSDHCAL_HITANALYSISMODULE_H
