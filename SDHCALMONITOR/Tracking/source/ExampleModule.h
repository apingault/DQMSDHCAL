// ExampleModule.h file
#include "dqm4hep/DQMAnalysisModule.h"
#include "DQMTrackAnalyzer.h"

using namespace dqm4hep;

class ExampleModule : public DQMAnalysisModule
{
public:
  ExampleModule();
 
  StatusCode initModule();
  StatusCode readSettings(const TiXmlHandle xmlHandle);
  StatusCode processEvent(DQMEvent *const pEvent);
  StatusCode startOfCycle();
  StatusCode endOfCycle();
  StatusCode startOfRun(DQMRun *const pRun);
  StatusCode endOfRun(DQMRun *const pRun);
  StatusCode endModule();
  DQMPlugin *clone() const { return new ExampleModule(); }

private:
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
