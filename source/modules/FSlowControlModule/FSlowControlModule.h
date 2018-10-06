/*
 *
 * FSlowControlModule.h header template automatically generated by a class generator
 * Creation date : mar. mars 8 2016
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


#ifndef DQMSDHCAL_FSLOWCONTROLMODULE_H
#define DQMSDHCAL_FSLOWCONTROLMODULE_H

#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMStandaloneModule.h"
#include "dqm4hep/DQMQualityTest.h"

class TGraph;

namespace dqm_sdhcal
{
/** LVInfo structure
 */
struct LVInfo
{
  float       m_vSet;   ///< The voltage as supplied by shifters
  float       m_vRead;  ///< The read voltage by the device
  float       m_iRead;  ///< The read current by the device
  std::string m_status; ///< The read current by the device
// curl -g http://lyosdhcal12:36000/FSLOW/CMD?name=LVSTATUS
// {"answer":{"STATUS":{"iout":0.0,"name":"ZUP","status":"OFF","vout":0.001000000047497451,"vset":5.943999767303467}},"status":"OK"}
};

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

/** HVInfo structure
 */
struct HVInfo
{
  int   m_chamberID; ///< The chamber ID
  float m_vSet;      ///< The voltage as supplied by shifters
  float m_iSet;      ///< The current as supplied by shifters
  float m_vRead;     ///< The read voltage by the device
  float m_iRead;     ///< The read current by the device

  // curl - g http:     //lyosdhcal12:36000/FSLOW/CMD?name=HVSTATUS\&first=0\&last=5
  // {"answer":
  //  {"STATUS":
  //     {"channels":[
  //     {"id":0,"iout":-25.99999618530273,"iset":0.0002500000118743628,"rampup":-25.99999618530273,"status":"WIENER-CRATE-MIB::outputStatus.u0 = BITS: 00 \n", "vout":-0.8530640006065369, "vset":1234.0},
  //     {"id":1,"..."}
  //        ]
  //  }
  // }
  // }
};

typedef std::map<unsigned int, HVInfo> HVInfoMap;

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

class CurrentQualityTest : public dqm4hep::DQMQualityTest
{
public:
  class Factory : public dqm4hep::DQMQualityTestFactory
  {
public:
    dqm4hep::DQMQualityTest *createQualityTest(const std::string& name) const;
  };

public:
  CurrentQualityTest(const std::string& name);
  dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);
  dqm4hep::StatusCode init();
  dqm4hep::StatusCode run(dqm4hep::DQMMonitorElement *pMonitorElement);
  bool canRun(dqm4hep::DQMMonitorElement *pMonitorElement) const;

private:
  int m_maxAllowedCurrent;
  int m_maxDangerousCurrent;
};

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

dqm4hep::DQMQualityTest *CurrentQualityTest::Factory::createQualityTest(const std::string& name) const
{
  return new CurrentQualityTest(name);
}


//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

/** FSlowControlModule class
 */
class FSlowControlModule : public dqm4hep::DQMStandaloneModule
{
public:

  /** Constructor
   */
  FSlowControlModule();

  /** Destructor
   */
  ~FSlowControlModule();

  dqm4hep::StatusCode initModule();
  dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle handle);
  dqm4hep::StatusCode startOfCycle();
  dqm4hep::StatusCode endOfCycle();
  dqm4hep::StatusCode endModule();
  dqm4hep::StatusCode process();

  std::string sendCmd(const char *cmd);
  float getOutputVoltage(uint32_t module, uint32_t voie);
  void getLVInfo(LVInfo& lvInfo);


private:

  /** Get the temperature as seen by the global device
   */
  float getGlobalTemperature();

  /** Get the pressure as seen by the global device
   */
  float getGlobalPressure();

  // curl -g http://lyosdhcal12:36000/FSLOW/CMD?name=PTSTATUS
  // {"answer":{"STATUS":{"name":"BMP","pressure":973.8900146484375,"status":"READ","temperature":22.50}},"status":"OK"}


	/** Get the high voltage infos for single chamber
	*/
  dqm4hep::StatusCode getHighVoltageInfos(HVInfo& hvInfo, int chamberID);

	/** Get the high voltage infos for all chambers
	*/
  dqm4hep::StatusCode getHighVoltageInfos(HVInfoMap& hvInfoMap);

  /** Get the low voltage info
   */
  dqm4hep::StatusCode getLowVoltageInfo(LVInfo& lvInfo);

  /** Configure graph
   */
  void configureGraph(TGraph *pGraph);

private:
  typedef std::map<unsigned int, dqm4hep::DQMMonitorElementPtr> DQMMonitorElementIDMap;

  // parameters
  std::string m_slowControlName;
  std::string m_hostName;
  std::string m_hostPort;
	unsigned int m_nChamber;
  std::string m_lvInfoName;
  std::string m_temperatureInfoName;
  std::string m_pressureInfoName;
  std::string m_hvInfoName;
  dqm4hep::StringVector m_hvInfoServiceNames;
  unsigned int m_globalDynamicGraphRange;

  // monitor elements
  dqm4hep::DQMMonitorElementPtr m_pGlobalTemperatureElement;
  dqm4hep::DQMMonitorElementPtr m_pGlobalPressureElement;

  dqm4hep::DQMMonitorElementPtr m_pHighVoltageVSetElement;
  dqm4hep::DQMMonitorElementPtr m_pHighVoltageVReadElement;
  dqm4hep::DQMMonitorElementPtr m_pHighVoltageVSetReadDiffElement;
  dqm4hep::DQMMonitorElementPtr m_pHighVoltageISetElement;
  dqm4hep::DQMMonitorElementPtr m_pHighVoltageIReadElement;
  dqm4hep::DQMMonitorElementPtr m_pHighVoltageISetReadDiffElement;

  dqm4hep::DQMMonitorElementPtr m_pLowVoltageElement;
  dqm4hep::DQMMonitorElementPtr m_pLowVoltageStatus;

  DQMMonitorElementIDMap m_chamberHVElementMap;
  DQMMonitorElementIDMap m_chamberCurrentElementMap;

  time_t m_startTime;
};
}

#endif  //  DQMSDHCAL_FSLOWCONTROLMODULE_H