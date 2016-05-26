/*
 *
 * ParticleIDModule.h header template automatically generated by a class generator
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


#ifndef DQMSDHCAL_PARTICLEIDMODULE_H
#define DQMSDHCAL_PARTICLEIDMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- std headers
#include <string>
#include <cstring>
#include <vector>

#include "DQMTriventModule.h"

namespace dqm_sdhcal
{

class EventClassifier;

/** ParticleIDModule class
 */ 
class ParticleIDModule : public DQMTriventModule
{
public:
	ParticleIDModule();
	~ParticleIDModule();

	dqm4hep::StatusCode userInitModule();
	dqm4hep::StatusCode userReadSettings(const dqm4hep::TiXmlHandle xmlHandle);
	dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

	dqm4hep::StatusCode startOfCycle();
	dqm4hep::StatusCode endOfCycle();
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endModule();

private:
	/**
	 */
	dqm4hep::StatusCode fillSummary();

	/**
	 */
	dqm4hep::StatusCode getNHits(EVENT::LCEvent *pLCEvent, unsigned int &nHits);

private:
	EventClassifier                      *m_pEventClassifier;
	dqm4hep::StringVector                 m_caloHitCollectionNames;

	dqm4hep::DQMMonitorElementPtr         m_pParticleIDSummary;
	dqm4hep::DQMMonitorElementPtr         m_pNHitNoise;
	dqm4hep::DQMMonitorElementPtr         m_pNHitCosmicMuons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitBeamMuons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitNeutralHadrons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitChargedHadrons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitElectrons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitPhotons;
	dqm4hep::DQMMonitorElementPtr         m_pNHitOthers;
	dqm4hep::DQMMonitorElementPtr         m_pNHitAll;

	unsigned int                          m_nNoiseWithinRun;
	unsigned int                          m_nCosmicMuonsWithinRun;
	unsigned int                          m_nBeamMuonWithinRun;
	unsigned int                          m_nChargedHadronsWithinRun;
	unsigned int                          m_nNeutralHadronsWithinRun;
	unsigned int                          m_nElectronsWithinRun;
	unsigned int                          m_nPhotonsWithinRun;
	unsigned int                          m_nOthersWithinRun;
}; 

} 

#endif  //  DQMSDHCAL_PARTICLEIDMODULE_H
