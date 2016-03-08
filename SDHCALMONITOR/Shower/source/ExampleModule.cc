// our header include
#include "ExampleModule.h"
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMRun.h"
#include "dqm4hep/DQMEvent.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMCoreTool.h"

#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/LCTOOLS.h"
#include "IMPL/LCRunHeaderImpl.h" 
#include "IMPL/LCEventImpl.h" 

#include "EVENT/RawCalorimeterHit.h" 
#include "IMPL/LCGenericObjectImpl.h"
#include <IMPL/LCCollectionVec.h>
#include "IMPL/RawCalorimeterHitImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
// #include "EVENT/SimTrackerHit.h" 
#include <UTIL/LCSplitWriter.h>
#include "UTIL/CellIDDecoder.h"
#include <IO/LCRunListener.h>
#include <IO/LCEventListener.h>


ExampleModule anExampleModule;


ExampleModule::ExampleModule() : DQMAnalysisModule("LaurentShowerModule")
{
setDetectorName("MySweetCalorimeter");
setVersion(1, 0, 0);
 m_configFile="/home/acqilc/SDHCAL/SDHCAL_Analysis/chinese_avril2015.xml";

}

StatusCode ExampleModule::readSettings(const TiXmlHandle xmlHandle)
{
  //  m_configFile = value.get("ConfigFile", m_configFile).asString();

  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readParameterValue(xmlHandle,"ConfigFile", m_configFile));

  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::initModule()
{
  theReader_= new DHCalEventReader();
  //theRootHandler_= new DCHistogramHandler();

  //  printf("Le file est %s \n",m_configFile.c_str());getchar();
  theReader_->ParseSteering(m_configFile);
  theReader_->setXdaqShift(24);
  theReader_->setDropFirstRU(false);
  theTrackAnalyzer_= new DQMTrackAnalyzer(theReader_,this);
  theTrackAnalyzer_->presetParameters();
  theTrackAnalyzer_->setrebuild(false);
  theReader_->registerAnalysis(theTrackAnalyzer_);

  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::processEvent(DQMEvent *pEvent)
{
  EVENT::LCEvent *pLCEvent = pEvent->getEvent<EVENT::LCEvent>();
  
  if(NULL == pLCEvent)
    return STATUS_CODE_FAILURE;
  
  theReader_->processEvent(pLCEvent);
  
  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::startOfCycle()
{
// no operation
  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::endOfCycle()
{
  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::startOfRun(DQMRun *pRun)
{
  std::cout << "Module : " << getName() << " -- startOfRun()" << std::endl;
  std::cout << "Run no " << pRun->getRunNumber() << std::endl;
  std::string timeStr;
  DQMCoreTool::timeToHMS(pRun->getStartTime(), timeStr);
  std::cout << "Starting time : " << timeStr << std::endl;
  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::endOfRun(DQMRun *pRun)
{
  std::cout << "Module : " << getName() << " -- startOfRun()" << std::endl;
  std::cout << "Run no " << pRun->getRunNumber() << std::endl;
  std::string timeStr;
  DQMCoreTool::timeToHMS(pRun->getStartTime(), timeStr);
  std::cout << "Ending time : " << timeStr << std::endl;
  return STATUS_CODE_SUCCESS;
}

StatusCode ExampleModule::endModule()
{
// no operation
  return STATUS_CODE_SUCCESS;
}
