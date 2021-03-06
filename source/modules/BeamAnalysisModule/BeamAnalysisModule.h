/// \file BeamAnalysisModule.h
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


#ifndef DQMSDHCAL_BEAMANALYSISMODULE_H
#define DQMSDHCAL_BEAMANALYSISMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMAnalysisModule.h"

// -- dqm sdhcal headers
#include "AnalysisTools.h"

// -- lcio headers
#include "lcio.h"

// -- std headers
#include <string>
#include <vector>


namespace dqm4hep { class TiXmlElement; class TiXmlHandle; }

namespace dqm_sdhcal
{

class EventHelper;

class BeamAnalysisModule : public dqm4hep::DQMAnalysisModule
{
public:
  BeamAnalysisModule();
  virtual ~BeamAnalysisModule();

private:
  // from analysis module
  dqm4hep::StatusCode initModule();
  dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);
  dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *const pEvent);

  dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
  dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
  dqm4hep::StatusCode startOfCycle();
  dqm4hep::StatusCode endOfCycle();
  dqm4hep::StatusCode endModule();

  void resetElements();

private:
  EventHelper                             *m_pEventHelper;
  EventHelper::EventParameters             m_eventParameters={};

  std::string                              m_inputCollectionName;
  std::string                              m_moduleLogStr;

  int                                      m_nEventProcessed;
  float                                    m_DAQ_BC_Period;
  unsigned int                             m_nParticleLastSpill;

  // Cuts
  int                                      m_skipEvent;
  double                                   m_newSpillTimeCut;


  // Monitor Elements
  dqm4hep::DQMMonitorElementPtr            m_pTimeDiffSpill;
  dqm4hep::DQMMonitorElementPtr            m_pTimeDiffTrigger;
  dqm4hep::DQMMonitorElementPtr            m_pTimeDiffTriggerToSpill;
  dqm4hep::DQMMonitorElementPtr            m_pTriggerPerSpill;
  dqm4hep::DQMMonitorElementPtr            m_pTriggerLastSpill;
  dqm4hep::DQMMonitorElementPtr            m_pSpillLength;
  dqm4hep::DQMMonitorElementPtr            m_pAcquisitionTime;

};
}
#endif // DQMSDHCAL_BEAMANALYSISMODULE_H
