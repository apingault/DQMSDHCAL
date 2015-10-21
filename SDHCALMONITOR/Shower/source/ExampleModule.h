// ExampleModule.h file
#include "dqm4hep/module/DQMAnalysisModule.h"
#include "DQMTrackAnalyzer.h"
using namespace dqm4hep;
class ExampleModule : public DQMAnalysisModule
{
public:
  ExampleModule();
 
  StatusCode initModule();
  StatusCode readSettings(const Json::Value &value);
  StatusCode processEvent(DQMEvent *pEvent);
  StatusCode startOfCycle();
  StatusCode endOfCycle();
  StatusCode startOfRun(DQMRun *pRun);
  StatusCode endOfRun(DQMRun *pRun);
  StatusCode endModule();
  DQMPlugin *clone() const { return new ExampleModule(); }

private:
  // elements
 
  DQMMonitorElement *m_pNumberOfHitsHistogram;
  DQMMonitorElement *m_pEnergyHistogram;
  DQMMonitorElement *m_pHitTimeWithinSpill;
  DQMMonitorElement *m_pXYHitPositionsHistogram;
  // additional parameters
  std::string  m_collectionName;
  std::string  m_configFile;
  unsigned int  m_minNHitToPublish;
  float  m_minHitPosition;
  float  m_maxHitPosition;
  DHCalEventReader* theReader_;
  DCHistogramHandler* theRootHandler_;
  DQMTrackAnalyzer* theTrackAnalyzer_;
};
