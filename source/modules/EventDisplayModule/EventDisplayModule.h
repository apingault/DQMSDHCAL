/*
 *
 * EventDisplayModule.h header template automatically generated by a class generator
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


#ifndef EVENTDISPLAYMODULE_H
#define EVENTDISPLAYMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- std headers
#include <string>
#include <cstring>
#include <vector>

#include "DQMTriventModule.h"


namespace dqm_sdhcal
{

/** EventDisplayModule class
 */ 
class EventDisplayModule : public DQMTriventModule
{
public:
	/** Constructor
	 */
	EventDisplayModule();

	/** Destructor
	 */
	virtual ~EventDisplayModule();

	dqm4hep::StatusCode userInitModule();
	dqm4hep::StatusCode userReadSettings(const dqm4hep::TiXmlHandle xmlHandle);
	dqm4hep::StatusCode processNoisyEvent(EVENT::LCEvent *pLCEvent);
	dqm4hep::StatusCode processPhysicalEvent(EVENT::LCEvent *pLCEvent);

	dqm4hep::StatusCode startOfCycle();
	dqm4hep::StatusCode endOfCycle();
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *const pRun);
	dqm4hep::StatusCode endModule();

private:
	std::string                           m_inputCollectionName;

	dqm4hep::DQMMonitorElementPtr         m_pEventDisplay3D;

	dqm4hep::DQMMonitorElementPtr         m_pLastProfileZX;
	dqm4hep::DQMMonitorElementPtr         m_pLastProfileZY;
	dqm4hep::DQMMonitorElementPtr         m_pLastProfileXY;

	dqm4hep::DQMMonitorElementPtr         m_pCycleStackedProfileZX;
	dqm4hep::DQMMonitorElementPtr         m_pCycleStackedProfileZY;
	dqm4hep::DQMMonitorElementPtr         m_pCycleStackedProfileXY;
}; 

} 

#endif  //  EVENTDISPLAYMODULE_H
