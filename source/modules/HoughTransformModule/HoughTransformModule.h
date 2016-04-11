/*
 *
 * HoughTransformModule.h header template automatically generated by a class generator
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


#ifndef HOUGHTRANSFORMMODULE_H
#define HOUGHTRANSFORMMODULE_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- std headers
#include <string>
#include <cstring>
#include <vector>

#include "DQMTriventModule.h"

namespace caloobject { class CaloTrack; }

namespace dqm_sdhcal
{

/** EventDisplayModule class
 */ 
class HoughTransformModule : public DQMTriventModule
{
public:
	/** Constructor
	 */
	HoughTransformModule();

	/** Destructor
	 */
	virtual ~HoughTransformModule();

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
	/** Analyze tracks and fill monitor elements
	 */
	dqm4hep::StatusCode analyzeTracks( const std::vector<caloobject::CaloTrack*> &tracks );

private:
	std::string                           m_inputCollectionName;

	dqm4hep::DQMMonitorElementPtr         m_pNRecTracks;
	dqm4hep::DQMMonitorElementPtr         m_pTrackChi2;
	dqm4hep::DQMMonitorElementPtr         m_pTrackLength;
	dqm4hep::DQMMonitorElementPtr         m_pClusterSize;
	dqm4hep::DQMMonitorElementPtr         m_pThetaTrack;
	dqm4hep::DQMMonitorElementPtr         m_pPhiTrack;
}; 

} 

#endif  //  HOUGHTRANSFORMMODULE_H