#ifndef _DQMTRACKANALYZER_H

#define _DQMTRACKANALYZER_H
#include <limits.h>
#include <vector>
#include "DHCALAnalyzer.h"
#include "DHCalEventReader.h"
#include "dqm4hep/module/DQMModule.h"
using namespace dqm4hep;

#include <iostream>
#include <sys/timeb.h>
#include "IMPL/LCTOOLS.h"
#include "EVENT/RawCalorimeterHit.h" 
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"
#include "UTIL/CellIDDecoder.h"
#include "DifGeom.h"
#include "ChamberGeom.h"

#include "TFile.h"

#include <vector>
#include <map>
#include "UtilDefs.h"
#include "Array3D.h"
#include "RecoHit.h"
#include "RECOCluster.h"
#include "RecoPoint.h"
#include "RecoCandTk.h"
#include "Shower.h"
#include "Amas.h"
#include "DQMSDHCALMonitor.h"
#include "libhoughStruct.h"
#include <ComputerTrack.h>

class DQMTrackAnalyzer : public DHCALAnalyzer
{
public:

  DQMTrackAnalyzer(DHCalEventReader* r,DQMModule* m);
  virtual ~DQMTrackAnalyzer(){;}
  void initialise();
  virtual void processEvent();
  virtual void initHistograms();
  virtual void processRunHeader()
  {
    if (writing_)
      reader_->writeRunHeader();
  }
  void presetParameters();
  void setWriting(bool t){writing_=t;}
  virtual void setReader(DHCalEventReader* r){reader_=r;}

  virtual void initJob();
  virtual void endJob();
  virtual void initRun(){;}
  virtual void endRun(){;}

  bool decodeTrigger(LCCollection* rhcol, double tcut);
  void setCollectionName(std::string s){ collectionName_=s;}
  unsigned long long getExternalTriggerTime() { return (unsigned long long) long(externalTriggerTime_);}
  unsigned long long getAbsoluteFrameTime(int bc) { return (unsigned long long) long(externalTriggerTime_-bc);}

  
  void setDropFirstSpillEvent(bool t){dropFirstSpillEvent_=t;}
  
  void setClockSynchCut(unsigned int t){clockSynchCut_=t;}
  void setSpillSize(double t){spillSize_=t;}
  
  uint32_t buildTracks(std::vector<RecoHit*> &vreco,std::string vdir="/Track");
  double checkTime();
  inline void setuseSynchronised(bool t){useSynchronised_=t; }
  inline void setoldAlgo(bool t){oldAlgo_=t;}
  inline bool getoldAlgo(){return oldAlgo_;}
  inline void settkMinPoint(int t){tkMinPoint_=t;}
  inline void settkExtMinPoint(int t){tkExtMinPoint_=t;}
  inline void settkBigClusterSize(int t){tkBigClusterSize_=t;}
  inline void setspillSize(float t){spillSize_=t;}
  inline void settkChi2Cut(float t){tkChi2Cut_=t;}
  inline void settkDistCut(float t){tkDistCut_=t;}
  inline void settkExtChi2Cut(float t){tkExtChi2Cut_=t;}
  inline void settkExtDistCut(float t){tkExtDistCut_=t;}
  inline void settkAngularCut(float t){tkAngularCut_=t;}
  inline void setclockSynchCut(int t){clockSynchCut_=t;}
  inline void setminChambersInTime(int t){minChambersInTime_=t;}
  inline void setmaxHitCount(int t){maxHitCount_=t;}
  inline void setchamberEdge(float t){chamberEdge_=t;}
  inline void setrebuild(bool t) {rebuild_=t;}
  inline void setcollectionName(std::string t){collectionName_=t;}

  
  void findTimeSeeds( IMPL::LCCollectionVec* rhcol, int32_t nhit_min,std::vector<uint32_t>& candidate);
  void processSeed(IMPL::LCCollectionVec* rhcol,uint32_t seed);

  uint32_t CerenkovTagger(uint32_t difid,uint32_t seed);
  uint32_t PMAnalysis(uint32_t difid);


private:


  int nAnalyzed_;
  int nInSynch_,run_;
  double externalTriggerTime_,lastSpill_,lastPowerPulsedTime_;


  std::map<uint32_t,uint32_t> asicCount_;
  std::map<uint32_t,double*> theCorreff_;
  double integratedTime_;
  int hrtype_;
  int lastrunnb_;
  unsigned int currentTime_,currentEvent_,lastSpyEvent_;
  // Control
  bool rebuild_;
  bool dropFirstSpillEvent_;
  bool findTracks_;
  bool useSynchronised_;
  bool oldAlgo_;
  int tkMinPoint_;
  int tkExtMinPoint_;
  int tkBigClusterSize_;
  float spillSize_;
  float tkChi2Cut_;
  float tkDistCut_;
  float tkExtChi2Cut_;
  float tkExtDistCut_;
  float tkAngularCut_;
  int clockSynchCut_;
  int minChambersInTime_;
  int maxHitCount_;
  int tkFirstChamber_;
  int tkLastChamber_;
  float chamberEdge_;
  uint32_t trackIndex_;
  uint32_t houghIndex_;
  bool useTk4_;
  std::string collectionName_;
  uint32_t offTimePrescale_;
  uint32_t theBIFId_;
  uint32_t theNplans_;
  
  // Reader
  DHCalEventReader* reader_;
  DQMModule* theModule_;
  bool writing_;
  bool headerWritten_;
  bool draw_;

  IMPL::LCEventImpl* evt_;

  struct timeb theTime_;
  struct timeb theCurrentTime_;

  double theRhcolTime_;
  double theTimeSortTime_;
  double theTrackingTime_;
  double theHistoTime_;
  double zLastAmas_;
  int theSeuil_;
  int32_t theSkip_,npi_;
  unsigned long long theBCID_;
  unsigned long long theBCIDSpill_;
  unsigned long long theAbsoluteTime_;

  uint32_t theDTC_,theGTC_;
  std::vector<RecoHit*> theHitVector_;
	  


  DQMSDHCALMonitor* theMonitoring_;
  uint32_t theMonitoringPeriod_;
  std::string theMonitoringPath_;
	
  uint32_t theNbShowers_,theNbTracks_;
  ComputerTrack* theComputerTrack_;
  unsigned long long theLastBCID_,theIdxSpill_;
  float theTimeInSpill_[20],theCountSpill_[20],theLastRate_;
  float coreRatio_;
  bool isNewSpill_,isPion_,isElectron_,isMuon_,isShower_,isProton_;
  uint32_t theCerenkovTag_,theNplansUsed_;
};
#endif
