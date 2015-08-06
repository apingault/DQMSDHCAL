  /// \file TestRawDataModule.h
/*
 *
 * TestRawDataModule.h header template automatically generated by a class generator
 * Creation date : mar. ao�t 4 2015
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


#ifndef TESTRAWDATAMODULE_H
#define TESTRAWDATAMODULE_H

// -- dqm4hep headers
#include "dqm4hep/core/DQM4HEP.h"
#include "dqm4hep/module/DQMModule.h"

namespace dqm4hep { class DQMMonitorElement; }

namespace dqm_sdhcal
{

class Streamout;

/** TestRawDataModule class
 */
class TestRawDataModule : public dqm4hep::DQMModule
{
public:
	/** Constructor
	 */
	TestRawDataModule();

	/** Destructor
	 */
	~TestRawDataModule();

	dqm4hep::StatusCode initModule();
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle &xmlHandle);
	dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *pEvent);
	dqm4hep::StatusCode startOfCycle();
	dqm4hep::StatusCode endOfCycle();
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *pRun);
	dqm4hep::StatusCode endOfRun(dqm4hep::DQMRun *pRun);
	dqm4hep::StatusCode resetModule();
	dqm4hep::StatusCode endModule();

protected:

	bool                   m_shouldProcessStreamout;
	bool                   m_shouldProcessTrivent;
	std::string             m_streamoutInputCollectionName;
	std::string             m_streamoutOutputCollectionName;
	std::string             m_triventInputCollectionName;
	std::string             m_triventOutputCollectionName;
//	std::string             m_

	dqm4hep::DQMMonitorElement       *m_pRawNumberOfElement;
	dqm4hep::DQMMonitorElement       *m_pNRawCalorimeterHits;
	dqm4hep::DQMMonitorElement       *m_pAmplitudeDistribution;
	dqm4hep::DQMMonitorElement       *m_pTimeDistribution;
	dqm4hep::DQMMonitorElement       *m_pICellDistribution;
	dqm4hep::DQMMonitorElement       *m_pJCellDistribution;
	dqm4hep::DQMMonitorElement       *m_pKCellDistribution;

	Streamout               *m_pStreamout;
}; 

} 

#endif  //  TESTRAWDATAMODULE_H
