/// \file NoiseAnalysisModule.h
/*
 *
 * AsicAnalysisModule.h header template automatically generated by a class generator
 * Creation date : mon. march 21 2016
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
 * @author Remi Ete, Laurent Mirabito, Antoine Pingault
 * @copyright CNRS , IPNL
 */


#ifndef DQMSDHCAL_NOISEANALYSISMODULE_H
#define DQMSDHCAL_NOISEANALYSISMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMAnalysisModule.h"
#include "dqm4hep/DQMElectronicsMapping.h"
#include "dqm4hep/DQMDataConverter.h"

// -- dqm sdhcal headers
#include "Geometry.h"
#include "AnalysisTools.h"

// -- lcio headers
#include "lcio.h"
#include "EVENT/RawCalorimeterHit.h"

// -- std headers
#include <string>
#include <vector>


namespace dqm4hep { class TiXmlElement; class TiXmlHandle; }

namespace dqm_sdhcal
{

class RawCaloHitObject
{
public:
  RawCaloHitObject(dqm4hep::DQMCartesianVector vec, int chanId, int asicId, int difId, int layerId, int threshold, int time, dqm4hep::DQMCartesianVector posShift);
  // algorithms assume that the zero is located at the middle of the first layer.

  ~RawCaloHitObject() {;}

  inline const dqm4hep::DQMCartesianVector getPosition() {return m_rawHitPosition;}
  inline const int getThreshold() {return m_threshold;}
  inline const float getTime() {return m_time;}
  inline const int getChannelId() {return m_chanId;}
  inline const int getAsicId() {return m_asicId;}
  inline const int getDifId() {return m_difId;}
  inline const int getLayerId() {return m_layerId;}

  typedef std::vector<RawCaloHitObject *> RawCaloHitList;

private:
  int m_chanId;
  int m_asicId;
  int m_difId;
  int m_layerId;
  dqm4hep::DQMCartesianVector m_rawHitPosition;
  float m_threshold;
  int m_time;
};

class Streamout;
class EventHelper;

class NoiseAnalysisModule : public dqm4hep::DQMAnalysisModule
{
public:
  NoiseAnalysisModule();
  virtual ~NoiseAnalysisModule();

private:
  // from analysis module
  dqm4hep::StatusCode initModule();
  dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);
  dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *const pEvent);
  dqm4hep::StatusCode performOutputDataConversion(EVENT::LCEvent *pLCEvent);

  dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
  dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
  dqm4hep::StatusCode startOfCycle();
  dqm4hep::StatusCode endOfCycle();
  dqm4hep::StatusCode endModule();

  dqm4hep::StatusCode findTrigger(EVENT::LCCollection* const pLCCollection);
  dqm4hep::StatusCode doDIFStudy(RawCaloHitObject * const pRawCaloHitObject);
  dqm4hep::StatusCode fillAsicOccupancyMap( RawCaloHitObject * const pRawCaloHitObject);
  dqm4hep::StatusCode doAsicStudy();
  int createAsicKey(int chanId, int difId, int asicId);
  void resetElements();

private:
  // converter parameters
  typedef dqm4hep::DQMDataConverter<EVENT::LCCollection, EVENT::LCCollection> CaloHitCollectionConverter;
  dqm4hep::StringVector                      m_rawCollectionNames;
  dqm4hep::StringVector                      m_recCollectionNames;
  dqm4hep::StringVector                      m_rawDataConverterNames;
  std::vector< CaloHitCollectionConverter *> m_dataConverters;
  unsigned short                             m_amplitudeBitRotation;
  // electronicsMapping parameters
  dqm4hep::DQMElectronicsMapping            *m_pElectronicsMapping;
  dqm4hep::DQMCartesianVector                m_cellReferencePosition;
  float m_cellSize0;
  float m_cellSize1;
  // Geometry parameters
  Geometry                                   m_geometry;
  DifMapping                                 m_difMapping;

private:
  EventHelper                                  *m_pEventHelper;
  EventClassifier                              *m_pEventClassifier;


private:
  // module parameters
  bool                                     m_shouldProcessStreamout;
  Streamout                               *m_pStreamout;

  std::string                              m_inputCollectionName;
  std::string                              m_cellIDDecoderString;
  std::string                              m_moduleLogStr;
  std::map<int, int>                       m_asicMap;


  unsigned int                             m_nActiveLayers;
  unsigned int                             m_nAsicPerDif;
  unsigned int                             m_nChanPerAsic;
  int                                      m_nStartLayerShift;

  int                                      m_nEventProcessed;
  double                                   m_eventIntegratedTime;
  double                                   m_spillIntegratedTime;
  double                                   m_totalIntegratedTime;
  unsigned long long                       m_hitTimeMin;
  unsigned long long                       m_hitTimeMax;

  double                                   m_timeLastTrigger;
  double                                   m_timeLastSpill;
  float                                    m_DAQ_BC_Period;
  unsigned int                             m_nParticleLastSpill;
  unsigned int                             m_nTrigger;
  // Cuts
  int                                      m_skipEvent;
  double                                   m_newSpillTimeCut;


  // Monitor Elements
  dqm4hep::DQMMonitorElementPtr               m_pTimeDiffSpill;
  dqm4hep::DQMMonitorElementPtr               m_pTimeDiffTrigger;
  dqm4hep::DQMMonitorElementPtr               m_pTimeDiffTriggerToSpill;
  dqm4hep::DQMMonitorElementPtr               m_pTriggerPerSpill;
  dqm4hep::DQMMonitorElementPtr               m_pTriggerLastSpill;
  dqm4hep::DQMMonitorElementPtr               m_pSpillLength;
  dqm4hep::DQMMonitorElementPtr               m_pAsicOccupancyAll;
  dqm4hep::DQMMonitorElementPtr               m_pAsicOccupancyChamber;
  dqm4hep::DQMMonitorElementPtr               m_pAsicOccupancyDIF;
  dqm4hep::DQMMonitorElementPtr               m_pAcquisitionTime;
  dqm4hep::DQMMonitorElementPtr               m_pHitFrequencyMap;

  struct DifElements
  {
    dqm4hep::DQMMonitorElementPtr             m_pAsicHits1;
    dqm4hep::DQMMonitorElementPtr             m_pAsicHits2;
    dqm4hep::DQMMonitorElementPtr             m_pAsicHits3;
    dqm4hep::DQMMonitorElementPtr             m_pAsicFreq1;
    dqm4hep::DQMMonitorElementPtr             m_pAsicFreq2;
    dqm4hep::DQMMonitorElementPtr             m_pAsicFreq3;

    dqm4hep::DQMMonitorElementPtr             m_pAsicOccupancy;
    dqm4hep::DQMMonitorElementPtr             m_pAsicOccupancyNumber;
    dqm4hep::DQMMonitorElementPtr             m_pAsicEventTime;
    dqm4hep::DQMMonitorElementPtr             m_pAsicEventTimeZoom;
  };
  struct LayerElements
  {
    dqm4hep::DQMMonitorElementPtr             m_pChamberHitsMap1;
    dqm4hep::DQMMonitorElementPtr             m_pChamberHitsMap2;
    dqm4hep::DQMMonitorElementPtr             m_pChamberHitsMap3;
    std::map<int, DifElements>                m_difElementMap;
  };
  std::map<unsigned int, LayerElements>       m_layerElementMap;
};
}
#endif // DQMSDHCAL_NOISEANALYSISMODULE_H
