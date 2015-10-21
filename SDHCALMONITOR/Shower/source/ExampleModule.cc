// our header include
#include "ExampleModule.h"
#include "dqm4hep/core/DQMMonitorElement.h"
#include "dqm4hep/core/DQMRun.h"
#include "dqm4hep/core/DQMEvent.h"
#include "dqm4hep/core/DQMXmlHelper.h"
#include "dqm4hep/module/DQMModuleApi.h"
#include "dqm4hep/core/DQMCoreTool.h"

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


ExampleModule::ExampleModule() : DQMAnalysisModule("ExampleModule")
{
setDetectorName("MySweetCalorimeter");
setVersion(1, 0, 0);
 m_configFile="/home/acqilc/SDHCAL/SDHCAL_Analysis/chinese_avril2015.xml";

}

StatusCode ExampleModule::readSettings(const Json::Value &value)
{
  m_configFile = value.get("ConfigFile", m_configFile).asString();

  //  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readValue(xmlHandle,"ConfigFile", m_configFile));
  /*
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readValue(xmlHandle,"CollectionName", m_collectionName));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readValue(xmlHandle,"MinNHitToPublish", m_minNHitToPublish));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readValue(xmlHandle,"MinHitPosition", m_minHitPosition));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMXmlHelper::readValue(xmlHandle,"MaxHitPosition", m_maxHitPosition));
  */
  return STATUS_CODE_SUCCESS;
}
StatusCode ExampleModule::initModule()
{
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::mkdir(this,"/Hits"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::mkdir(this,"/Energy"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::mkdir(this,"/Time"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::cd(this, "/Hits"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,DQMModuleApi::bookIntHistogram1D(this,m_pNumberOfHitsHistogram, "NumberOfHits", "Number of hits", 1501, 0, 1500));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::cd(this,"/Energy"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,DQMModuleApi::bookRealHistogram1D(this,m_pEnergyHistogram, "HitEnergy", "Hits energy", 101, 0, 100));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DQMModuleApi::cd(this, "/Time"));
  RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,DQMModuleApi::bookRealHistogram1D(this,m_pHitTimeWithinSpill, "HitTimeWithinSpill", "Hit time within a spill",101, 0, 100));
  m_pHitTimeWithinSpill->setResetPolicy(END_OF_CYCLE_RESET_POLICY);
  DQMModuleApi::cd(this);
  DQMModuleApi::ls(this, true);



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

  EVENT::LCCollection *pCaloHitCollection = pLCEvent->getCollection(m_collectionName);
  for(unsigned int e=0 ; e<pCaloHitCollection->getNumberOfElements() ; e++)
    {
      EVENT::CalorimeterHit *pCaloHit = (EVENT::CalorimeterHit*) pCaloHitCollection->getElementAt(e);
      if(NULL == pCaloHit)
	continue;
      m_pEnergyHistogram->get<TH1F>()->Fill(pCaloHit->getEnergy());
      m_pHitTimeWithinSpill->get<TH1F>()->Fill(pCaloHit->getTime());
    }
  m_pNumberOfHitsHistogram->get<TH1I>()->Fill(pCaloHitCollection->getNumberOfElements());
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

  double meanNHit = m_pNumberOfHitsHistogram->get<TH1I>()->GetMean();
  if(meanNHit > 160 && meanNHit < 180)
    m_pNumberOfHitsHistogram->setQuality(GOOD_QUALITY);
  else
    m_pNumberOfHitsHistogram->setQuality(BAD_QUALITY);
  if(m_pNumberOfHitsHistogram->get<TH1I>()->GetEntries() < 500)
    m_pNumberOfHitsHistogram->setToPublish(false);
  else
    m_pNumberOfHitsHistogram->setToPublish(true);
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
