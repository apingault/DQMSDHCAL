#define NPLANS_USED 6
#include "dqm4hep/core/DQMMonitorElement.h"
#include "dqm4hep/core/DQMRun.h"
#include "dqm4hep/core/DQMXmlHelper.h"
#include "dqm4hep/module/DQMModuleApi.h"
#include "dqm4hep/core/DQMCoreTool.h"
#include "DQMShowerAnalyzer.h"
#include "DIFUnpacker.h"
#include <TLine.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitter.h>
#include <TF1.h>
#include <TPluginManager.h>
#include <stdint.h>
#include <math.h>
#include <ext/hash_map>
using namespace __gnu_cxx;
#include <boost/pool/poolfwd.hpp>
#include <boost/pool/singleton_pool.hpp>
#include "TPolyLine3D.h"
#include "TVirtualPad.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
//#include <lapacke.h>
#include "DHShower.h"
using namespace boost; 
typedef boost::singleton_pool<RecoHit, sizeof(RecoHit)> RecoHitPool;
typedef std::vector<RecoHit*>::iterator recit;

#include <HoughLocal.h>

//#include "ChamberAnalyzer.h"

#define posError 0.5


uint32_t DQMShowerAnalyzer::CerenkovTagger(uint32_t difid,uint32_t seed)
{
  uint32_t tag=0;
  //printf(" The Seed %d \n",seed);
  float fs=seed*1.;
  for (std::vector<DIFPtr*>::iterator it = reader_->getDIFList().begin();it!=reader_->getDIFList().end();it++)
    {
      DIFPtr* d = (*it);
      if (d->getID()!=difid) continue;
      // Loop on frames
      
      for (uint32_t i=0;i<d->getNumberOfFrames();i++)
  	{
  	  //if (abs(d->getFrameTimeToTrigger(i)-seed)<200)
	  //printf("\t Cerenkov %d  %f\n",d->getFrameTimeToTrigger(i),fabs(seed-d->getFrameTimeToTrigger(i)));
	  float tf=d->getFrameTimeToTrigger(i);
  	  if (fabs(fs-tf)<10)
  	    {

	      
	      for (uint32_t j=0;j<64;j++)
		{
		  if (d->getFrameLevel(i,j,0)) tag +=1;
		  if (d->getFrameLevel(i,j,1)) tag +=2;
		}
	      return tag;

  	    }
  	}

	    
    }
  //  getchar();
  return 0;
}

void DQMShowerAnalyzer::initHistograms()
{
  //  rootHandler_->BookTH1("/Clusters/EST1",100,0.,300.);
}


DQMShowerAnalyzer::DQMShowerAnalyzer(DHCalEventReader* r,DQMModule* m) :trackIndex_(0),nAnalyzed_(0),clockSynchCut_(8), spillSize_(90000),maxHitCount_(500000),
									     tkMinPoint_(3),tkExtMinPoint_(3),tkBigClusterSize_(32),tkChi2Cut_(0.01),tkDistCut_(5.),tkExtChi2Cut_(0.01),tkExtDistCut_(10.),tkAngularCut_(20.),zLastAmas_(134.),
									     findTracks_(true),dropFirstSpillEvent_(false),useSynchronised_(true),chamberEdge_(5.),rebuild_(false),oldAlgo_(true),collectionName_("DHCALRawHits"),
									     tkFirstChamber_(1),tkLastChamber_(61),useTk4_(false),offTimePrescale_(1),houghIndex_(0),theRhcolTime_(0.),theTimeSortTime_(0.),theTrackingTime_(0),
								      theHistoTime_(0),theSeuil_(0),draw_(false),theSkip_(0),theMonitoring_(NULL),theMonitoringPeriod_(0),theMonitoringPath_("/dev/shm/Monitoring"),theBCIDSpill_(0),theLastBCID_(0),theBIFId_(333),theNplansUsed_(6)
{
  reader_=r;
  theModule_ =m;
  this->initialise();

  HoughCut cuts;
  theComputerTrack_=new ComputerTrack(&cuts);

  theComputerTrack_->DefaultCuts();
  theMonitoring_=new DQMSDHCALMonitor(theModule_);

}

void DQMShowerAnalyzer::initialise()
{
  headerWritten_=false;
  //  TVirtualFitter::SetDefaultFitter("Minuit"); 
  //  gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit", "TMinuitMinimizer", "Minuit", "TMinuitMinimizer(const char *)"); 
  integratedTime_=0;
  asicCount_.clear();
  DCBufferReader::setDAQ_BC_Period(0.2);
  
  theTime_.time=0;
  theTime_.millitm=0;


}
void DQMShowerAnalyzer::initJob(){presetParameters();}
void DQMShowerAnalyzer::endJob(){
 
	
}
void DQMShowerAnalyzer::presetParameters()
{
  //  printf("On rentre dans preset \n");getchar();
  std::map<std::string,MarlinParameter> m=reader_->getMarlinParameterMap();
  std::map<std::string,MarlinParameter>::iterator it;
  try
    {
      if ((it=m.find("ClockSynchCut"))!=m.end()) clockSynchCut_=it->second.getIntValue();
      if ((it=m.find("SpillSize"))!=m.end()) spillSize_=it->second.getIntValue();
      if ((it=m.find("MaxHitCount"))!=m.end()) maxHitCount_=it->second.getIntValue();
      if ((it=m.find("MinChambersInTime"))!=m.end()) minChambersInTime_=it->second.getIntValue();
      if ((it=m.find("TkMinPoint"))!=m.end()) tkMinPoint_=it->second.getIntValue();
      if ((it=m.find("TkExtMinPoint"))!=m.end()) tkExtMinPoint_=it->second.getIntValue();
      if ((it=m.find("TkBigClusterSize"))!=m.end()) tkBigClusterSize_=it->second.getIntValue();
      if ((it=m.find("TkChi2Cut"))!=m.end()) tkChi2Cut_=it->second.getDoubleValue();
      if ((it=m.find("TkDistCut"))!=m.end()) tkDistCut_=it->second.getDoubleValue();
      if ((it=m.find("TkExtChi2Cut"))!=m.end()) tkExtChi2Cut_=it->second.getDoubleValue();
      if ((it=m.find("TkExtDistCut"))!=m.end()) tkExtDistCut_=it->second.getDoubleValue();
      if ((it=m.find("TkAngularCut"))!=m.end()) tkAngularCut_=it->second.getDoubleValue();
      if ((it=m.find("ChamberEdge"))!=m.end()) chamberEdge_=it->second.getDoubleValue();
      if ((it=m.find("FindTracks"))!=m.end()) findTracks_=it->second.getBoolValue();
      if ((it=m.find("DropFirstSpillEvent"))!=m.end()) dropFirstSpillEvent_=it->second.getBoolValue();
      if ((it=m.find("UseSynchronised"))!=m.end()) useSynchronised_=it->second.getBoolValue();
      if ((it=m.find("UseTk4"))!=m.end()) useTk4_=it->second.getBoolValue();
      if ((it=m.find("Rebuild"))!=m.end()) rebuild_=it->second.getBoolValue();
      if ((it=m.find("OldAlgo"))!=m.end()) oldAlgo_=it->second.getBoolValue();
      if ((it=m.find("CollectionName"))!=m.end()) collectionName_=it->second.getStringValue();
      if ((it=m.find("TkFirstChamber"))!=m.end()) tkFirstChamber_=it->second.getIntValue();
      if ((it=m.find("TkLastChamber"))!=m.end()) tkLastChamber_=it->second.getIntValue();
      if ((it=m.find("OffTimePrescale"))!=m.end()) offTimePrescale_=it->second.getIntValue();
      if ((it=m.find("Seuil"))!=m.end()) theSeuil_=it->second.getIntValue();
      if ((it=m.find("Interactif"))!=m.end()) draw_=it->second.getBoolValue();
      if ((it=m.find("SkipEvents"))!=m.end()) theSkip_=it->second.getIntValue();
      if ((it=m.find("zLastAmas"))!=m.end()) zLastAmas_=it->second.getDoubleValue();
      if ((it=m.find("MonitoringPath"))!=m.end()) theMonitoringPath_=it->second.getStringValue();
      if ((it=m.find("MonitoringPeriod"))!=m.end()) theMonitoringPeriod_=it->second.getIntValue();
      if ((it=m.find("BIFId"))!=m.end()) theBIFId_=it->second.getIntValue();
      if ((it=m.find("PlansUsed"))!=m.end()) theNplansUsed_=it->second.getIntValue();

      //   printf("Apres preset %d\n",theBIFId_);getchar();
      DEBUG_PRINT("Interactif %d \n",draw_);

      //getchar();

    }
  catch (std::string s)
    {
      std::cout<<__PRETTY_FUNCTION__<<" error "<<s<<std::endl;
    }
	
}
bool DQMShowerAnalyzer::decodeTrigger(LCCollection* rhcol, double tcut)
{
  // if (rhcol->getNumberOfElements()==0) return true;

  // Find Trigger information
  IntVec vTrigger;IMPL::RawCalorimeterHitImpl* hit;
  unsigned int difid=0;
  // Find the first read DIF id for this trigger

  if (rhcol->getNumberOfElements()!=0)
    {
		
      try {
	hit = (IMPL::RawCalorimeterHitImpl*) rhcol->getElementAt(0);
      }
      catch (std::exception e)
	{
	  std::cout<<"No hits "<<std::endl;
	  return false;
	}
      if (hit!=0) 
	difid = hit->getCellID0()&0xFF;
    }

  if (difid==0) return false;

  //Find the parameters
  std::stringstream pname("");
  pname <<"DIF"<<difid<<"_Triggers";

  rhcol->getParameters().getIntVals(pname.str(),vTrigger);

  if (vTrigger.size()==0) return false; 
  //for (int i=0;i<vTrigger.size();i++)
  ///  std::cout<<vTrigger[i]<<std::endl;

  // Decode Large Bunch Crossing
  unsigned long long Shift=16777216ULL;//to shift the value from the 24 first bits

  unsigned long long  lbc=0;
  unsigned long long  lbci=0;
  uint32_t  lb5=vTrigger[4] ;
  uint32_t  lb4=vTrigger[3] ;

  lbc = lb5*Shift+ lb4;

  theDTC_=vTrigger[0];
  theGTC_=vTrigger[1];
  theBCID_=lbc;

  //lbc =lb4*Shift+lb5;
  double tTrigger_= lbc*(DCBufferReader::getDAQ_BC_Period()*1E-6);
  DEBUG_PRINT("Time stamp ==========================>: %d %d %llu %f\n",lb4,lb5,lbc,tTrigger_);
  // Fill Trigger info
  // DEBUG_PRINT("creqtion de htspill \n");

  TH1* htspill=theMonitoring_->getRealHistogram1D("SpillDif","spilLDif",500,0.,100.)->get<TH1F>();
  //DEBUG_PRINT("apres creqtion de htspill \n");

  // Calculate tiem differences since the last trigger
  double tdif = tTrigger_-externalTriggerTime_;
  //std::cout<<lbc<<" "<<externalTriggerTime_<<" "<<tdif<<std::endl;
#ifdef DEBUG
  if (tdif>50 || tdif <-1E-3)
    {
      cout<<tdif << " strange time  "<<externalTriggerTime_<<endl;
      //streamlog_out(DEBUG)<<lb4<<endl;
      //streamlog_out(DEBUG)<<lb5<<endl;
      //streamlog_out(DEBUG)<<lbc<<endl;
		
    }
  //streamlog_out(DEBUG)<<lbc <<" "<<tdif<<" # hits "<<rhcol->getNumberOfElements()<<std::endl;
#endif
  isNewSpill_=(tdif>tcut);
  if (tdif>tcut) 
    {
      lastSpill_=tTrigger_;
      std::cout<<"New Spill "<<tdif<<"===========================================================>"<<npi_<<std::endl; 
      npi_=0;
      htspill->Fill(tdif);
    }
  externalTriggerTime_=tTrigger_;
  //  for (unsigned int i=0;i<vTrigger.size();i++) streamlog_out(MESSAGE)<<i<<" "<<vTrigger[i]<<std::endl;




  // Drop the first event of the Spill
  //  streamlog_out(MESSAGE)<<dropFirstSpillEvent_<<std::endl;
  if (tdif>tcut && dropFirstSpillEvent_) return false;
  if ((tTrigger_-lastSpill_)<1. && dropFirstSpillEvent_) 
    {
      DEBUG_PRINT("Event dropped %f %f \n",tTrigger_,lastSpill_);
      return false;
    }
  // TH1* htdiff= rootHandler_->GetTH1("TimeDif");
  // if (htdiff==NULL)
  //     {
  //         htdiff =rootHandler_->BookTH1( "TimeDif",20000,0.,20000.);
  //     }
  //    htdiff->Fill(tdif*1000.);




  return true;
}


void DQMShowerAnalyzer::findTimeSeeds( IMPL::LCCollectionVec* rhcol, int32_t nhit_min,std::vector<uint32_t>& candidate)
{
  map<uint32_t,uint32_t> tcount;
  map<uint32_t,int32_t> tedge;

  // Tcount is the time histo
  for (uint32_t i=0;i<rhcol->getNumberOfElements();i++)
    {
      IMPL::RawCalorimeterHitImpl* hit = (IMPL::RawCalorimeterHitImpl*) rhcol->getElementAt(i);
      if (hit==0) continue;
      uint32_t bc = hit->getTimeStamp();
      map<uint32_t,uint32_t>::iterator it=tcount.find(bc);
      if (it!=tcount.end()) 
	it->second=it->second+1;
      else
	{
	  std::pair<uint32_t,uint32_t> p(bc,1);
	  tcount.insert(p);
	}
    }
	
  std::vector<uint32_t> seed;
  seed.clear();
	
  //d::cout<<"Size =>"<<tcount.size()<<std::endl;
  // Tedge is convolute with +1 -1 +1 apply to tcount[i-1],tcount[i],tcount[i+1]
  for (map<uint32_t,uint32_t>::iterator it=tcount.begin();it!=tcount.end();it++)
    {
      //std::cout<<it->first<<" "<<it->second<<std::endl;
		
      map<uint32_t,uint32_t>::iterator ita=tcount.find(it->first+1);
      map<uint32_t,uint32_t>::iterator itb=tcount.find(it->first-1);
      int32_t c=-1*it->second;
      if (ita!=tcount.end()) c+=ita->second;
      if (itb!=tcount.end()) c+=itb->second;
      std::pair<uint32_t,int32_t> p(it->first,c);
      tedge.insert(p);
		
    }
  //d::cout<<"Size Edge =>"<<tedge.size()<<std::endl;
  // Now ask for a minimal number of hits
  uint32_t nshti=0;
  for (map<uint32_t,int32_t>::iterator it=tedge.begin();it!=tedge.end();)
    {
      //std::cout<<it->first<<"====>"<<it->second<<" count="<<tcount[it->first]<<std::endl;
      if (it->second<-1*(nhit_min-2))
	{
			
	  //std::cout<<it->first<<"====>"<<it->second<<" count="<<tcount[it->first]<<std::endl;

	  seed.push_back(it->first);
	  it++;
	}
      else
	tedge.erase(it++);
    }
	
  // for (std::vector<uint32_t>::iterator is=seed.begin();is!=seed.end();is++)
  //   std::cout<<" seed " <<(*is)<<" count "<<tcount[(*is)]<<std::endl      ;
  // Merge adjacent seeds
  candidate.clear();
  for (uint32_t i=0;i<seed.size();)
    {
      if ((i+1)<=(seed.size()-1))
	{
	  if (seed[i+1]-seed[i]<=5)
	    {
	      //candidate.push_back(int((seed[i+1]+seed[i])/2));
	      uint32_t max_c=0;
	      uint32_t max_it=0;
	      uint32_t imin=seed[i];
	      uint32_t imax=seed[i+1];
	      if (seed[i+1]>seed[i])
		{
		}
	      for (uint32_t it=imin;it<=imax;it++)
		{
		  if (tcount.find(it)==tcount.end()) continue;
		  if (tcount[it]>max_c) {max_c=tcount[it];max_it=it;}
		}
	      if (max_it!=0)
		candidate.push_back(max_it);
	      else
		candidate.push_back(seed[i]);
	      i+=2;
	    }
	  else
	    {
	      candidate.push_back(seed[i]);
	      i++;
	    }
	}
      else
	{
	  candidate.push_back(seed[i]);
	  i++;
	}

		
    }
  //td::cout<<candidate.size()<<" good showers "<< tedge.size()<<std::endl;
  std::sort(candidate.begin(),candidate.end(),std::greater<uint32_t>());

  /*  
      for (std::vector<uint32_t>::iterator is=candidate.begin();is!=candidate.end();is++)
    std::cout<<(*is)<<" ---> "<<tcount[(*is)]<<std::endl;
  */
  return ;
}

void DQMShowerAnalyzer::processSeed(IMPL::LCCollectionVec* rhcol,uint32_t seed)
{
  ShowerParams ish;
  //std::vector<RecoHit*> vrh;

  currentTime_=seed;
  
  theAbsoluteTime_=theBCID_-currentTime_;
  if (theBCIDSpill_==0) theBCIDSpill_=theAbsoluteTime_;
  if (theAbsoluteTime_-theBCIDSpill_>5./2E-7) theBCIDSpill_=theAbsoluteTime_;
 
    DEBUG_PRINT("GTC %d DTC %d BCID %llu Current Time %llu Time SPill %f Distance %f \n",theGTC_,theDTC_,theBCID_,currentTime_,theBCIDSpill_*2E-7,(theAbsoluteTime_-theBCIDSpill_)*2E-7);
  // DEBUG_PRINT("Building voulume for %d \n",seed);
  //uint32_t nhits=buildVolume(rhcol,seed);
  //  DEBUG_ DEBUG_PRINT("1");
  std::map<uint32_t,std::vector<IMPL::RawCalorimeterHitImpl*> >::iterator iseed=reader_->getPhysicsEventMap().find(seed);
   if (iseed==reader_->getPhysicsEventMap().end()) 
   {
      INFO_PRINT("Impossible \n");
      return ;
   }
   // More than 15 hits required
   if (iseed->second.size()<15) return;
   // Fill 3d MAP with hits
  uint32_t nhits=buildVolume(seed);

  if (nhits<15) return;

  //DEBUG_PRINT("Edge detection for %d \n",seed);
  // Construct the volume around each hit
  buildEdges();
  theNplans_=0;
  std::bitset<48> asics[255];
  // DEBUG_PRINT("2");
  for (int ii=0;ii<255;ii++) asics[ii].reset();
  // Fill the Hit vector to calculate PCA axis
  theHitVector_.clear();
  std::vector<RecoHit*> tkhits;
  tkhits.clear();
  for (uint32_t k=0;k<60;k++)
    { bool found=false;
      for (uint32_t i=0;i<96;i++)
	for (uint32_t j=0;j<96;j++)
	  if (theImage_.getValue(k,i,j)>0) 
	    {
	      // DEBUG_PRINT("%d%d %d %d \n",i,j,hitVolume_[k][i][j].getFlag(RecoHit::EDGE),hitVolume_[k][i][j].getFlag(RecoHit::ISOLATED));
	      //if ((hitVolume_[k][i][j].getFlag(RecoHit::EDGE)==1||hitVolume_[k][i][j].getFlag(RecoHit::ISOLATED)==1) )
	      if ((hitVolume_[k][i][j].getFlag(RecoHit::ISOLATED)!=1) )
		tkhits.push_back(&hitVolume_[k][i][j]);
	      theHitVector_.push_back(&hitVolume_[k][i][j]);
	      asics[hitVolume_[k][i][j].dif()].set(hitVolume_[k][i][j].getAsic()-1);
	      found=true;}
      if (found) theNplans_++;
    }
  //INFO_PRINT("On a trouve                                    %d hits                    %d plans -> %d \n",theHitVector_.size(),theNplans_,minChambersInTime_);
  // DEBUG_PRINT("3");
  if (theHitVector_.size()<15) return;
  // Check the number of plans hits
  if (theNplans_<minChambersInTime_) return;  


  // Cerenkov Tag
  theCerenkovTag_=this->CerenkovTagger(theBIFID_,seed);

  //printf("TAG=======================> %d \n",tag);
  TH1* hnoctag= rootHandler_->GetTH1("NoCTag");
  TH1* hctag1= rootHandler_->GetTH1("CTag1");
  TH1* hctag2= rootHandler_->GetTH1("CTag2");
  TH1* hctag3= rootHandler_->GetTH1("CTag3");
  TH1* hctag3notk= rootHandler_->GetTH1("CTag3Notk");
  if (hnoctag==NULL)
    {
      hnoctag =rootHandler_->BookTH1( "NoCTag",1000,0.,3000.);
      hctag1 =rootHandler_->BookTH1( "CTag1",1000,0.,3000.);
      hctag2 =rootHandler_->BookTH1( "CTag2",1000,0.,3000.);
      hctag3 =rootHandler_->BookTH1( "CTag3",1000,0.,3000.);
      hctag3notk =rootHandler_->BookTH1( "CTag3Notk",1000,0.,3000.);

    }

  theMonitoring_->cd("/Cerenkov");
  TH1* hnoctag=theMonitoring_->getRealHistogram1D("NoCTag","Nb hit BIF no Tag ",100,0.1,300.1)->get<TH1F>(); 
  TH1* hctag1=theMonitoring_->getRealHistogram1D("CTag1","Nb hit BIF Tag 1",100,0.1,300.1)->get<TH1F>(); 
  TH1* hctag2= theMonitoring_->getRealHistogram1D("CTag2","Nb hit BIF Tag 2 ",100,0.1,300.1)->get<TH1F>(); 
  TH1* hctag3=theMonitoring_->getRealHistogram1D("CTag3 ","Nb Hit BIF Tag 3",100,0.1,300.1)->get<TH1F>(); 
  if (theCerenkovTag_==0)
    hnoctag->Fill(theHitVector_.size());
  if (theCerenkovTag_==1)
    hctag1->Fill(theHitVector_.size());
  if (theCerenkovTag_==2)
    hctag2->Fill(theHitVector_.size());
  if (theCerenkovTag_==3)
    hctag3->Fill(theHitVector_.size());

 
  // Calculate PCA with all hits but ISOLATED

  Shower::computePrincipalComponents(tkhits,(double*) &ish);

  // remove shower looking like tracks
  if (sqrt((ish.lambda[0])/ish.lambda[2])<thePCARatioCut_) 
    {theNbTracks_++;return;}

  // Remove single shower outside calo
  if (ish.xm[0]<5) return;

  // Now calculate entry point
  double* v=ish.l2;
  double p[3];
  for (uint32_t i=0;i<3;i++) p[i]=v[i]/ish.lambda[2];
  double* x=ish.xm;    
  RecoCandTk t;
  t.ax_ =p[0]/p[2];
  t.ay_ =p[1]/p[2];
  t.bx_=x[0]-t.ax_*x[2];
  t.by_=x[1]-t.ay_*x[2];
  // remove edge showers
  if (t.bx_<5 || t.bx_>95) return;
  if (t.by_<5 || t.by_>95) return;



  // Make a shower analysis
  uint32_t nshower=buildClusters(theHitVector_);
  //INFO_PRINT("After buildCluster for %d \n",seed);
  if (isPion_)
    {
      printf("Pion seed %d \n",seed);
    }
  else
    if (isShower_)
       printf("Shower seed %d \n",seed);
  theEvent_.tracklength=theComputerTrack_->Length()*1.;

  uint32_t counts[3][5];
  uint32_t ir;
  memset(counts,0,15*sizeof(uint32_t));
  for (std::vector<RecoHit*>::iterator ih=theHitVector_.begin();ih!=theHitVector_.end();ih++)
    {
      if ((*ih)->getFlag(RecoHit::THR0)) ir=0;
      if ((*ih)->getFlag(RecoHit::THR1)) ir=1;
      if ((*ih)->getFlag(RecoHit::THR2)) ir=2;
      counts[ir][0]++;
      if ((*ih)->getFlag(RecoHit::MIP)!=0) counts[ir][1]++;
      if ((*ih)->getFlag(RecoHit::CORE)!=0) counts[ir][2]++;
      if ((*ih)->getFlag(RecoHit::EDGE)!=0) counts[ir][3]++;
      if ((*ih)->getFlag(RecoHit::ISOLATED)!=0) counts[ir][4]++;
    }
  /*
  theEvent_.m0=counts[0][1];
  theEvent_.c0=counts[0][2];
  theEvent_.e0=counts[0][3];
  theEvent_.i0=counts[0][4];
  theEvent_.m1=counts[1][1];
  theEvent_.c1=counts[1][2];
  theEvent_.e1=counts[1][3];
  theEvent_.i1=counts[1][4];
  theEvent_.m2=counts[2][1];
  theEvent_.c2=counts[2][2];
  theEvent_.e2=counts[2][3];
  theEvent_.i2=counts[2][4];

  for (uint32_t ir=0;ir<3;ir++)
    {
      for (uint32_t ic=0;ic<5;ic++)
	 DEBUG_PRINT("%d ",counts[ir][ic]);
       DEBUG_PRINT("\n");
    }
  */
  this->ShowerBuilder(theHitVector_);
  //if (nshower>0) getchar();
  // if (nshower>0) 
  // 	{
  // 	   DEBUG_PRINT("Angles %f %f \n",t.ax_,t.ay_);
  // 	  INFO_PRINT(" %d Hits Composantes principales %f %f %f => %.2f  and %.2f  \n",vrh.size(),ish.lambda[0],ish.lambda[1],ish.lambda[2], sqrt((ish.lambda[0]+ish.lambda[1])/ish.lambda[2]), sqrt((ish.lambda[0])/ish.lambda[2]));
  // 	  getchar();
  // 	}
  return;


  /*

  ShowerParams ish;
  ShowerParams isha;



  currentTime_=seed;
  
  theAbsoluteTime_=theBCID_-currentTime_;
  if (theBCIDSpill_==0) theBCIDSpill_=theAbsoluteTime_;
  if (theAbsoluteTime_-theBCIDSpill_>15/2E-7) theBCIDSpill_=theAbsoluteTime_;
 
    DEBUG_PRINT("GTC %d DTC %d BCID %llu Current Time %llu Time SPill %f Distance %f \n",theGTC_,theDTC_,theBCID_,currentTime_,theBCIDSpill_*2E-7,(theAbsoluteTime_-theBCIDSpill_)*2E-7);

  int nhits=0;
  theNplans_=0;
  std::bitset<60> chhit(0);
  std::map<uint32_t,std::vector<IMPL::RawCalorimeterHitImpl*> >::iterator iseed=reader_->getPhysicsEventMap().find(seed);
   if (iseed==reader_->getPhysicsEventMap().end()) 
   {
      INFO_PRINT("Impossible \n");
      return ;
   }
   if (theBIFId_!=0)
     {
       theCerenkovTag_=this->CerenkovTagger(theBIFId_,seed);
     }
   //printf("%d %d \n",seed,theCerenkovTag_);
   if (theHitVector_.size()!=0)
     for (std::vector<RecoHit*>::iterator ih=theHitVector_.begin();ih!=theHitVector_.end();ih++) delete (*ih);

   theHitVector_.clear();
   for (std::vector<IMPL::RawCalorimeterHitImpl*>::iterator ihit=iseed->second.begin();ihit!=iseed->second.end();ihit++)
    {
      IMPL::RawCalorimeterHitImpl* hit =(*ihit);
      unsigned int difid = hit->getCellID0()&0xFF;
      if (difid<1 || difid>255) continue;
      int asicid = (hit->getCellID0()&0xFF00)>>8;
      int channel= (hit->getCellID0()&0x3F0000)>>16;
      bool thr[3];
      //      DEBUG_PRINT("%x \n",hit->getCellID0());
      int ithr= hit->getAmplitude()&0x3;
      if (ithr==0)
	{
	  std::cout<<difid<<" had:"<<asicid<<":"<<channel<<":"<<ithr<<std::endl;
	  continue;
	}
      std::map<unsigned int,DifGeom>::iterator idg = reader_->getDifMap().find(difid);
      DifGeom& difgeom = idg->second;
      int x=0,y=0;
      uint32_t chid = idg->second.getChamberId();

       if ((difid==3|| chid==10) && theCerenkovTag_==0)
	{
	  printf("C'est quoi ce bordel %d \n",seed);
	  for (std::vector<DIFPtr*>::iterator it = reader_->getDIFList().begin();it!=reader_->getDIFList().end();it++)
	    {
	      DIFPtr* d = (*it);
	      if (d->getID()!=difid) continue;
      // Loop on frames
      
	      for (uint32_t i=0;i<d->getNumberOfFrames();i++)
		{
  	  
		  printf("\t Cerenkov %d \n",d->getFrameTimeToTrigger(i));
		}
	    }
	  getchar();
	}
      uint32_t hrtype=2;
            DifGeom::PadConvert(asicid,channel,x,y,hrtype);
      uint32_t I=difgeom.toGlobalX(x);
      if (I<1 || I>96) continue;
      uint32_t J=difgeom.toGlobalY(y);
      if (J<1 || J>96) continue;
      if (chid<1 || chid>60) continue;
      chhit.set(chid,1);
      RecoHit* h=  new RecoHit();
      //planes.set(chid-1,true);
      std::map<unsigned int,ChamberGeom>::iterator icg = reader_->getChamberMap().find( chid);
      ChamberGeom& chgeom = icg->second;
      //printf("Hit beeing filled %d %d %d\n",chid-1,I-1,J-1);
      chgeom.setZ(reader_->getPosition(chid).getZ0());

      h->initialise(difgeom,chgeom,hit,hrtype);
      theHitVector_.push_back(h);

      nhits++;
    }


  if (nhits==0) return;
  theNplans_=chhit.count();
  // if (nhits<30) return;
  //return;
  //DEBUG_PRINT("Edge detection for %d \n",seed);
  
  //printf("TAG=======================> %d \n",tag);
  if (theBIFId_!=0)
    {
      theMonitoring_->cd("/Cerenkov");
      TH1* hnoctag=theMonitoring_->getRealHistogram1D("NoCTag","Nb hit BIF no Tag ",100,0.1,300.1)->get<TH1F>(); 
      TH1* hctag1=theMonitoring_->getRealHistogram1D("CTag1","Nb hit BIF Tag 1",100,0.1,300.1)->get<TH1F>(); 
      TH1* hctag2= theMonitoring_->getRealHistogram1D("CTag2","Nb hit BIF Tag 2 ",100,0.1,300.1)->get<TH1F>(); 
      TH1* hctag3=theMonitoring_->getRealHistogram1D("CTag3 ","Nb Hit BIF Tag 3",100,0.1,300.1)->get<TH1F>(); 
      TH1* hpatnotag= theMonitoring_->getRealHistogram1D("PatternNoTag","Plane pattern no BIF Tag ",100,0.1,100.1)->get<TH1F>();
      TH1* hpattag= theMonitoring_->getRealHistogram1D("PatternTag","Plane Pattern BIF tag ",100,0.1,100.1)->get<TH1F>();

  
      if (theCerenkovTag_==0)
	{
	  //std::cout<<chhit<<"->no tag"<<std::endl;
	  hnoctag->Fill(theHitVector_.size());
	  for(int i=0;i<64;i++)
	    {
	      if (chhit[i]!=0) hpatnotag->Fill(i*1.);
	    }
	  hpatnotag->Fill(100.);
	}
      else
	{
	  for(int i=0;i<64;i++)
	    {
	      if (chhit[i]!=0) hpattag->Fill(i*1.);
	    }
	  hpattag->Fill(100.);
	}
      if (theCerenkovTag_==1)
	hctag1->Fill(theHitVector_.size());
      if (theCerenkovTag_==2)
	hctag2->Fill(theHitVector_.size());
      if (theCerenkovTag_==3)
	hctag3->Fill(theHitVector_.size());
    }
 
  Shower::computePrincipalComponents(theHitVector_,(double*) &isha);
  this->buildTracks(theHitVector_,"/TrackNoCut");
  if (theBIFId_!=0)
    {
      if (theCerenkovTag_!=0)
	{
	  
	  if (theCerenkovTag_==1 || theCerenkovTag_==3)
	    this->buildTracks(theHitVector_,"/TrackPM1");
	  if (theCerenkovTag_==2|| theCerenkovTag_==3)
	    this->buildTracks(theHitVector_,"/TrackPM2");
	  if (theCerenkovTag_==3)
	    this->buildTracks(theHitVector_,"/TrackPM_BOTH");
	  
	  //INFO_PRINT(" Mean event parameter  %d %f %f %f => %f %d TKS \n",theHitVector_.size(),isha.lambda[0],isha.lambda[1],isha.lambda[2],sqrt((isha.lambda[0]+isha.lambda[1])/isha.lambda[2]),theComputerTrack_->getTracks().size()); 
	  
	  //  DEBUG_PRINT("6\n");
	  if (theComputerTrack_->getTracks().size()>0) theCerenkovTag_+=4;
	  
	  return;
	  
	}
      else
	if (theCerenkovTag_==0)
	  {
	    this->buildTracks(theHitVector_,"/TrackNOPM");
	    return;
	  }
    } 
  return;
  */
}
uint32_t DQMShowerAnalyzer::buildVolume(uint32_t seed)
{
 std::map<uint32_t,std::vector<IMPL::RawCalorimeterHitImpl*> >::iterator iseed=reader_->getPhysicsEventMap().find(seed);
 if (iseed==reader_->getPhysicsEventMap().end()) 
   {
      DEBUG_PRINT("Impossible \n");
   return 0;
 }
  //memset(hitVolume_,0,60*96*96*sizeof(RecoHit)); Aie!!!

  //return 0;
 std::bitset<61> planes(0);
 theImage_.clear();
 theImageWeight_.clear();
  

  uint32_t ncount=0;
  
  if (iseed->second.size()>4096) return 0;
  for (std::vector<IMPL::RawCalorimeterHitImpl*>::iterator ihit=iseed->second.begin();ihit!=iseed->second.end();ihit++)
    {
      IMPL::RawCalorimeterHitImpl* hit =(*ihit);
      unsigned int difid = hit->getCellID0()&0xFF;
      if (difid<1 || difid>255) continue;
      int asicid = (hit->getCellID0()&0xFF00)>>8;
      int channel= (hit->getCellID0()&0x3F0000)>>16;
      bool thr[3];
      //      DEBUG_PRINT("%x \n",hit->getCellID0());
      int ithr= hit->getAmplitude()&0x3;
      if (ithr==0)
	{
	  std::cout<<difid<<" had:"<<asicid<<":"<<channel<<":"<<ithr<<std::endl;
	  continue;
	}
      std::map<unsigned int,DifGeom>::iterator idg = reader_->getDifMap().find(difid);
      DifGeom& difgeom = idg->second;
      int x=0,y=0;
      uint32_t chid = idg->second.getChamberId();
      uint32_t hrtype=2;
      if (chid>50) hrtype=11;
      DifGeom::PadConvert(asicid,channel,x,y,hrtype);
      uint32_t I=difgeom.toGlobalX(x);
      if (I<1 || I>96) continue;
      uint32_t J=difgeom.toGlobalY(y);
      if (J<1 || J>96) continue;
      if (chid<1 || chid>60) continue;
		
      //planes.set(chid-1,true);
      theImage_.setValue(chid-1,I-1,J-1,1);
      theImageWeight_.setValue(chid-1,I-1,J-1,(1<<ithr-1));
      std::map<unsigned int,ChamberGeom>::iterator icg = reader_->getChamberMap().find( chid);
      ChamberGeom& chgeom = icg->second;
      //printf("Hit beeing filled %d %d %d\n",chid-1,I-1,J-1);
      chgeom.setZ(reader_->getPosition(chid).getZ0());

      hitVolume_[chid-1][I-1][J-1].initialise(difgeom,chgeom,hit,hrtype);

      ncount++;
    }
  DEBUG_PRINT("Total number of Hit in buildVolume %d  %d \n",ncount,seed);
  return ncount;
}

void DQMShowerAnalyzer::processEvent()
{


  if (reader_->getEvent()==0) return;
  
  evt_=reader_->getEvent();
  //theSkip_=380;
  if (evt_->getEventNumber()<=theSkip_) return;
  printf("Processing %d - %d \n",evt_->getRunNumber(),evt_->getEventNumber());

 
  nAnalyzed_++;
  IMPL::LCCollectionVec* rhcol=NULL;
  bool rhcoltransient=false;
  try {
    rhcol=(IMPL::LCCollectionVec*) evt_->getCollection(collectionName_);
    rebuild_=false;
  }
  catch (...)
    {
      try 
	{
	  evt_->getCollection("RU_XDAQ");
	  rebuild_=true;
	}
      catch (...)
	{
	   DEBUG_PRINT("No raw data or calo hits \n");
	  exit(0);
	}
    }
  if (rebuild_)
    {
      reader_->parseRawEvent();
      DEBUG_PRINT("End of parseraw \n");
      std::vector<uint32_t> seed;
      if (useSynchronised_ )
	{
	  //DEBUG_PRINT("Calling FastFlag2\n");
			
	  reader_->findTimeSeeds(minChambersInTime_,seed);
	  // DEBUG_PRINT("End of FastFlag2 \n");
			
	}
      else
	{
			
	  seed.clear();
	}
      //
      // Use DIF timing
      seed.clear();
      //INFO_PRINT("Calling CreaetRaw %d\n",minChambersInTime_);
      reader_->findDIFSeeds(minDIFsInTime_);
      rhcol=reader_->createRawCalorimeterHits(reader_->getDIFSeeds());
      //rhcol=reader_->createRawCalorimeterHits(seed);
      evt_->addCollection(rhcol,"DHCALRawHits");
      rhcoltransient=false; 

    }
  else
    rhcol=(IMPL::LCCollectionVec*) evt_->getCollection(collectionName_);

  DEBUG_PRINT("End of CreaetRaw %d \n",rhcol->getNumberOfElements());  

  // Limit to 4E6 hits
  if (rhcol==NULL) return;
  if (rhcol->getNumberOfElements()>4E6) return;
  if (rhcol->getNumberOfElements()==0) return;


  //  LCTOOLS::printParameters(rhcol->parameters());
  //DEBUG_PRINT("Time Stamp %d \n",evt_->getTimeStamp());


  //DEBUG_PRINT("Calling decodeTrigger\n");
  // TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if (!decodeTrigger(rhcol,spillSize_) ) { if (rhcoltransient) delete rhcol;return;}
  
  //if (isNewSpill_) return;
  // if (evt_->getEventNumber()%100 ==0)
  //   rootHandler_->writeSQL();
  //    rootHandler_->writeXML(theMonitoringPath_);


  reader_->findTimeSeeds(minChambersInTime_);

  std::vector<uint32_t> vseeds=reader_->getTimeSeeds();

  
   INFO_PRINT("================>  %d  Number of seeds %d \n",evt_->getEventNumber(),(int) vseeds.size());

  if (vseeds.size()==0)  { if (rhcoltransient) delete rhcol;return;}
 
  theNbShowers_=0;
  theNbTracks_=0;
  bool hasPion=false;
  for (uint32_t is=0;is<vseeds.size();is++)
    {

      this->processSeed(rhcol,vseeds[is]);


      
    }

  if ((theBCID_-theLastBCID_)*2E-7>5)
    {
      theBCIDSpill_=theBCID_;
      theIdxSpill_=0;
      memset(theCountSpill_,0,20*sizeof(float));
      memset(theTimeInSpill_,0,20*sizeof(float));
      printf("===============================================================================================> NEW SPILL : %f\n",theBCIDSpill_*2E-7);
    }
  else
    theIdxSpill_+=1;

  //  theCountSpill_[theIdxSpill_%20] =  theNbShowers_+theNbTracks_;
  theCountSpill_[theIdxSpill_%10] =  theNbShowers_;
  theTimeInSpill_[theIdxSpill_%10] = theMonitoring_->getEventIntegratedTime()*2E-7;

  // Integrated 10 last
  float nc=0;
  float tc=0;
  DEBUG_PRINT("showers %uud %d %d ",theIdxSpill_,theNbShowers_,theNbTracks_);
  for (int i=0;i<10;i++)
    {
      nc+=theCountSpill_[i];
      tc+=theTimeInSpill_[i];
      //printf("%f ",theCountSpill_[i]);
    }

  //INFO_PRINT("\n %d Number of showers/tracks %d,%d Event time %f -> Absolute bcid  %f-> Rate %f %f %f\n",evt_->getEventNumber(),theNbShowers_,theNbTracks_,theMonitoring_->getEventIntegratedTime()*2E-7,(theBCID_-theBCIDSpill_)*2E-7,nc,tc,nc/tc);
  theLastRate_=nc/tc;
  theLastBCID_=theBCID_;
  //etchar();

  if (rhcoltransient) delete rhcol;

  return;
}  



#define DBG printf("%d %s\n",__LINE__,__PRETTY_FUNCTION__);


void DQMShowerAnalyzer::buildEdges()
{
  theMonitoring_->cd("/Showers");
  TH1* hweight=theMonitoring_->getRealHistogram1D("HitWeight","Hit local weight",110,-0.1,0.99)->get<TH1F>();
  TH1* hcoreratio=theMonitoring_->getRealHistogram1D("Coreratio","Core / (Edge+Isolated)",600,0.,30.)->get<TH1F>();
  uint32_t nmax=0;
  uint32_t nedge=0,ncore=0,niso=0;
  int32_t ixmin=-6,ixmax=6; // 6 avant
  for (uint32_t k=0;k<theImage_.getXSize();k++)
    {

      for (uint32_t i=0;i<theImage_.getYSize();i++)
	{
	  //if (image2x[k][i]>=-1 ) continue;
			
	  for (uint32_t j=0;j<theImage_.getZSize();j++)

	    {
	      if (theImage_.getValue(k,i,j)>0)
		{
		  int izmin=-2;
		  int izmax=+2;
		  if (k==0) {izmin=0;izmax=4;}
		  if (k==1) {izmin=0;izmax=4;}
		  uint32_t nv=0;
		  std::vector<RecoHit*> vnear_;
		  RecoHit* h0=&hitVolume_[k][i][j];
		  for (int z=izmin;z<=izmax;z++)
		    {
		      if (z+k<0) continue;
		      if (z+k>=theImage_.getXSize()) continue;
		      for (int x=ixmin;x<=ixmax;x++)
			{
			  if (x+i<0) continue;
			  if (x+i>=theImage_.getYSize()) continue;
			for (int y=ixmin;y<=ixmax;y++)
			  {
			    if (y+j<0) continue;
			    if (y+j>=theImage_.getZSize()) continue;
			    if (theImage_.getValue(k+z,i+x,j+y)<=0) continue;
			    RecoHit* h1=&hitVolume_[k+z][i+x][j+y];
			    float x0=h0->X(),y0=h0->Y(),z0=h0->Z(),x1=h1->X(),y1=h1->Y(),z1=h1->Z();
			    float dist=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
			    if (dist<6.) vnear_.push_back(h1); // was 6

			    nv++;
			  }
			}
		    }
		  if (nv>nmax) nmax=nv;
		  hitVolume_[k][i][j].calculateComponents(vnear_);

		  Components* c=(Components*) hitVolume_[k][i][j].Components();
		  double w=0;
		  if (c->l2!=0) w=sqrt((c->l1)/c->l2);
		  hweight->Fill(w);

		  if (c->l2==0)
		    {hitVolume_[k][i][j].setFlag(RecoHit::ISOLATED,true);niso++;}
		  else 
		    if (w<0.3) // 0.3 before  
		      {hitVolume_[k][i][j].setFlag(RecoHit::EDGE,true);nedge++;}
		    else
		      {hitVolume_[k][i][j].setFlag(RecoHit::CORE,true);ncore++;}
		}

	    }
	}
	  
    }
  if ((nedge+niso)>0)
    {
      coreRatio_=ncore*1./(nedge+niso);
      hcoreratio->Fill(ncore*1./(nedge+niso));

      DEBUG_PRINT("Max neigh %d Hit summary Isolated %d Edge %d Core %d => %.2f\n",nmax,niso,nedge,ncore,(nedge*100./(nedge+ncore)));
    }
  return;
}
uint32_t DQMShowerAnalyzer::buildClusters(std::vector<RecoHit*> &vrh)
{
  
  std::vector<RECOCluster*> clusters;
  std::vector<RECOCluster*> realc;
  std::vector<RECOCluster*> intc;
  clusters.clear();
  realc.clear();
  intc.clear();
  HoughCut cuts;
  float *h_x=(float *) malloc(4096*sizeof(float));
  float *h_y= (float *) malloc(4096*sizeof(float));
  float *h_z=(float *) malloc(4096*sizeof(float));
  unsigned int *h_layer=(unsigned int *) malloc(4096*sizeof(unsigned int));
  uint32_t nshower=0;
  uint32_t nh0=0,nh1=0,nh2=0;
  //ComputerTrack ch(&cuts);
  //ch.DefaultCuts();
#define CORE_EDGE_METHOD
#ifdef CORE_EDGE_METHOD
 TH1* hwt = rootHandler_->GetTH1("/Clusters/hitweight");

  if (hwt==NULL)
    {

      hwt =rootHandler_->BookTH1("/Clusters/hitweight",500,0.,1.);
    }
  for (std::vector<RecoHit*>::iterator ih=vrh.begin();ih<vrh.end();ih++)
    {
     if ((*ih)->getAmplitude()==2) nh0++;
      if ((*ih)->getAmplitude()==1) nh1++;
      if ((*ih)->getAmplitude()==3) nh2++;
      Components* ch=(Components*) (*ih)->Components();
      double w=0;
      if (ch->l2!=0) w=sqrt((ch->l1)/ch->l2);
      hwt->Fill(w);
      
      if ((*ih)->getFlag(RecoHit::CORE)==1) continue;
      bool merged=false;
      for (std::vector<RECOCluster*>::iterator ic=realc.begin();ic!=realc.end();ic++)
	{
	  if ((*ih)->chamber()!=(*ic)->chamber()) continue;
	  merged=(*ic)->Append((*(*ih)),2.); // avant 4 et normalement 2
	  if (merged) break;
	}
      if (merged) continue;
      RECOCluster* c= new RECOCluster((*(*ih)));
      realc.push_back(c);
      clusters.push_back(c);
    }
  uint32_t fpi=99,lpi=0;
  for (std::vector<RecoHit*>::iterator ih=vrh.begin();ih<vrh.end();ih++)
    {
      if ((*ih)->getFlag(RecoHit::CORE)!=1) continue;
      bool merged=false;
      for (std::vector<RECOCluster*>::iterator ic=intc.begin();ic!=intc.end();ic++)
	{
	  if ((*ih)->chamber()!=(*ic)->chamber()) continue;
	  merged=(*ic)->Append((*(*ih)),4.); // avant 4 et normalement 2
	  if (merged) break;
	}
      if (merged) continue;
      RECOCluster* c= new RECOCluster((*(*ih)));
      intc.push_back(c);
      clusters.push_back(c);
    }

  // Move small cluster from intc to realc
  for (std::vector<RECOCluster*>::iterator ic=intc.begin();ic!=intc.end();)
    {
      if ((*ic)->getHits()->size()>=5) 
	++ic;
      else
	{
	  realc.push_back((*ic));
	  intc.erase(ic);
	}
    }
  for (std::vector<RECOCluster*>::iterator ic=realc.begin();ic!=realc.end();)
    {
      if ((*ic)->getHits()->size()<5) 
	++ic;
      else
	{
	  intc.push_back((*ic));
	  realc.erase(ic);
	}
    }
  // Now clear clusters >6

   for (std::vector<RECOCluster*>::iterator ic=intc.begin();ic!=intc.end();ic++)
	{
	      uint32_t ch=(*ic)->chamber();
	      if (lpi<ch) lpi=ch;
	      if (fpi>ch) fpi=ch;
	}

#endif
   
   if (false && (realc.size()*100./vrh.size()>40 || realc.size()*100./clusters.size()<95.) )
     {
       
       for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
	 delete (*ic);
       free(h_x);
       free(h_y);
       free(h_z);
       free(h_layer);
       return nshower;
     }
   
   // 2D draw of Hits and real clusters
#define DRAW_HISTOS
#ifdef  DRAW_HISTOS
#define BOOK_HISTOS
#define FILL_HISTOS
#endif
   //#define BOOK_HISTOS

#ifdef BOOK_HISTOS
  
  TH2F* hzx = (TH2F*) rootHandler_->GetTH2("/Clusters/HZXPOS");

  TH2F* hzy = (TH2F*) rootHandler_->GetTH2("/Clusters/HZY");
  TH2F* czx = (TH2F*) rootHandler_->GetTH2("/Clusters/CZX");
  TH2F* czy = (TH2F*) rootHandler_->GetTH2("/Clusters/CZY");
  TH2F* mzx = (TH2F*) rootHandler_->GetTH2("/Clusters/MZX");
  TH2F* mzy = (TH2F*) rootHandler_->GetTH2("/Clusters/MZY");
  TH2F* izx = (TH2F*) rootHandler_->GetTH2("/Clusters/IZX");
  TH2F* izy = (TH2F*) rootHandler_->GetTH2("/Clusters/IZY");
  
  TH1* hest = rootHandler_->GetTH1("/Clusters/MipRatio");
  
  TH1* hest1 = rootHandler_->GetTH1("/Clusters/ShowerNumberOfHit");
  TH1* hest2 = rootHandler_->GetTH1("/Clusters/LowRateNumberOfHit");
  TH1* hest3 = rootHandler_->GetTH1("/Clusters/HighRateNumberOfHit");
  TH1* hlom = rootHandler_->GetTH1("/Clusters/TrackLengthOverMip");
  TProfile* hrate=(TProfile*) rootHandler_->GetTH1("/Clusters/RATE");
  TProfile* hrate0=(TProfile*) rootHandler_->GetTH1("/Clusters/RATE0");
  TProfile* hrate1=(TProfile*) rootHandler_->GetTH1("/Clusters/RATE1");
  TProfile* hrate2=(TProfile*) rootHandler_->GetTH1("/Clusters/RATE2");

  TProfile* hspill=(TProfile*) rootHandler_->GetTH1("/Clusters/MeanSpillRate");
  TH1* hhadron = rootHandler_->GetTH1("/Clusters/hadrons");
  TH1* hpion = rootHandler_->GetTH1("/Clusters/pions");
  TH1* hproton = rootHandler_->GetTH1("/Clusters/protons");
  TH1* helectron = rootHandler_->GetTH1("/Clusters/electrons");
  TH1* hmuon = rootHandler_->GetTH1("/Clusters/muons");

  if (hest1==NULL)
    {
       DEBUG_PRINT("Booking\n");
      hhadron =rootHandler_->BookTH1("/Clusters/hadrons",1000,0.,3000.);
      hpion =rootHandler_->BookTH1("/Clusters/pions",1000,0.,3000.);
      hproton =rootHandler_->BookTH1("/Clusters/protons",1000,0.,3000.);
      helectron =rootHandler_->BookTH1("/Clusters/electrons",1000,0.,3000.);
      hmuon =rootHandler_->BookTH1("/Clusters/muons",1000,0.,1000.);
      hest1 =rootHandler_->BookTH1("/Clusters/ShowerNumberOfHit",1000,0.,3000.);
      hest2 =rootHandler_->BookTH1("/Clusters/LowRateNumberOfHit",1000,0.,3000.);
      hest3 =rootHandler_->BookTH1("/Clusters/HighRateNumberOfHit",1000,0.,3000.);
      hrate =rootHandler_->BookProfile("/Clusters/RATE",50,0.,500.,0.,2000.);
      hrate0 =rootHandler_->BookProfile("/Clusters/RATE0",50,0.,500.,0.,2000.);
      hrate1 =rootHandler_->BookProfile("/Clusters/RATE1",50,0.,500.,0.,2000.);
      hrate2 =rootHandler_->BookProfile("/Clusters/RATE2",50,0.,500.,0.,2000.);
      hspill =rootHandler_->BookProfile("/Clusters/MeanSpillRate",80,0.,20.,0.,2000.);
  
      hzx =(TH2F*)rootHandler_->BookTH2("/Clusters/HZXPOS",120,0.,170.,120,-10.,110.);


      hzy =(TH2F*)rootHandler_->BookTH2("/Clusters/HZY",120,0.,170.,120,-10.,110.);
    
      czy =(TH2F*)rootHandler_->BookTH2("/Clusters/CZY",120,0.,170.,120,-10.,110.);
      czx =(TH2F*)rootHandler_->BookTH2("/Clusters/CZX",120,0.,170.,120,-10.,110.);
      mzy =(TH2F*)rootHandler_->BookTH2("/Clusters/MZY",120,0.,170.,120,-10.,110.);
      mzx =(TH2F*)rootHandler_->BookTH2("/Clusters/MZX",120,0.,170.,120,-10.,110.);
      izy =(TH2F*)rootHandler_->BookTH2("/Clusters/IZY",120,0.,170.,120,-10.,110.);
      izx =(TH2F*)rootHandler_->BookTH2("/Clusters/IZX",120,0.,170.,120,-10.,110.);
      hest = rootHandler_->BookTH1("/Clusters/MipRatio",120,-0.1,1.1);
      hlom = rootHandler_->BookTH1("/Clusters/TrackLengthOverMip",1100,-0.1,2.1);

      
    }
  
  else
    {
      
  

     
      hzx->Reset();
     hzy->Reset();
     czy->Reset();
      czx->Reset();
      izy->Reset();
      izx->Reset();
      mzy->Reset();
      mzx->Reset();
     
    }
   
#endif
#ifdef FILL_HISTOS 
  for (std::vector<RecoHit*>::iterator ih=vrh.begin();ih<vrh.end();ih++)
    {
      hzx->Fill((*ih)->Z(),(*ih)->X());
      hzy->Fill((*ih)->Z(),(*ih)->Y());
    }
#endif

  uint32_t nstub=0;
  for (std::vector<RECOCluster*>::iterator ic=realc.begin();ic!=realc.end();ic++)
    {
#ifdef FILL_HISTOS
      czx->Fill((*ic)->Z(),(*ic)->X());
      czy->Fill((*ic)->Z(),(*ic)->Y());
#endif
      h_x[nstub]=(*ic)->X();
      h_y[nstub]=(*ic)->Y();
      h_z[nstub]=(*ic)->Z();
      h_layer[nstub]=(*ic)->chamber();
      nstub++;
    }
  theComputerTrack_->associate(nstub,h_x,h_y,h_z,h_layer);
 
  uint32_t nmip=0;
   for (unsigned int i=0;i<theComputerTrack_->getCandidates().size();i++)
	{
	  RecoCandTk& tk = theComputerTrack_->getCandidates()[i];
#ifdef CORRECT_GEOM
	  tk.ax_/=1.04125;
	  tk.bx_/=1.04125;
	  tk.ay_/=1.04125;
	  tk.by_/=1.04125;
#endif
	  // DEBUG_PRINT("Track cand tk: %f %f %f %f \n",tk.ax_,tk.bx_,tk.ay_,tk.by_);
// 	  for (std::vector<RecoHit*>::iterator ic=vrh.begin();ic!=vrh.end();ic++)
// 	    {
// 	      float dx=tk.getXext((*ic)->Z())-(*ic)->X();
// 	      float dy=tk.getYext((*ic)->Z())-(*ic)->Y();
// 	      if (sqrt(dx*dx+dy*dy)<2.)
// 		{
// 		  (*ic)->setFlag(RecoHit::MIP,true);
// 		  (*ic)->setFlag(RecoHit::ISOLATED,false);
// #ifdef FILL_HISTOS
// 		  mzx->Fill((*ic)->Z(),(*ic)->X());
// 		  mzy->Fill((*ic)->Z(),(*ic)->Y());
// #endif
// 		}
// 	    }
	  for (std::vector<RECOCluster*>::iterator ic=realc.begin();ic!=realc.end();ic++)
	    {
	      float dx=tk.getXext((*ic)->Z())-(*ic)->X();
	      float dy=tk.getYext((*ic)->Z())-(*ic)->Y();
	      if (sqrt(dx*dx+dy*dy)<2.)
		{
		  // DEBUG_PRINT("Point (%f,%f,%f) dist %f %f -> %f  \n",(*ic)->X(),(*ic)->Y(),(*ic)->Z(),dx,dy,sqrt(dx*dx+dy*dy));
		for (std::vector<RecoHit>::iterator ih=(*ic)->getHits()->begin();ih!=(*ic)->getHits()->end();ih++)
		  {
		    if ((*ih).getFlag(RecoHit::MIP)==0)
		      {
			// DEBUG_PRINT("Point (%f,%f,%f) dist %f %f -> %f  \n",(*ic)->X(),(*ic)->Y(),(*ic)->Z(),dx,dy,sqrt(dx*dx+dy*dy));
			(*ih).setFlag(RecoHit::MIP,true);
			hitVolume_[ih->chamber()-1][ih->I()-1][ih->J()-1].setFlag(RecoHit::MIP,true);
#ifdef FILL_HISTOS
			mzx->Fill(ih->Z(),ih->X());
			mzy->Fill(ih->Z(),ih->Y());
#endif
			nmip++;
		      }
		  }
		}
	    }
	}
   DEBUG_PRINT("==> MIPS hit %d -> %.2f Length %f \n",nmip,nmip*100./vrh.size(),theComputerTrack_->Length()); 
   float lom=theComputerTrack_->Length()*1./vrh.size();
   hlom->Fill(lom);
   bool electron=(theComputerTrack_->Length()<5. || nmip*100./vrh.size()<1);
   bool muon=nmip*100./vrh.size()>60.;
   bool pion = (lom>1E-2 && lom<0.5);
   //if (false && (electron || muon) )
   if (false &&(lom<1E-4 || lom>0.75))
   //if (lom>1E-4) //select electrons
    {
       for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
	 delete (*ic);
       free(h_x);
       free(h_y);
       free(h_z);
       free(h_layer);
       return nshower;
    }
  //#undef DRAW_HISTOS  
  /*
  ch.ComputeOneShot(nstub,h_x,h_y,h_z,h_layer);
  std::vector<RecoCandTk> &v=ch.getCandidates();
  for (std::vector<RecoCandTk>::iterator it=v.begin();it!=v.end();it++)
     DEBUG_PRINT("%f %f %f %f \n",it->ax_,it->bx_,it->ay_,it->by_);
  */
  //drawph(ch.getPh());
   //std::vector<Amas> theAmas_;
  theAmas_.clear();

  for (std::vector<RECOCluster*>::iterator ic=intc.begin();ic!=intc.end();ic++)
    {
      for (std::vector<RecoHit>::iterator ih=(*ic)->getHits()->begin();ih!=(*ic)->getHits()->end();ih++)
	{
#ifdef FILL_HISTOS
	  izx->Fill(ih->Z(),ih->X());
	  izy->Fill(ih->Z(),ih->Y());
#endif
	  RecoHit& h=(*ih);
	  bool appended=false;
	  for (std::vector<Amas>::iterator ia=theAmas_.begin();ia!=theAmas_.end();ia++)
	    {
	      appended=(ia->append(&h,2));
	      if (appended) break;
	    }
	  if (!appended)
	    {
	      Amas a(&h);
	      theAmas_.push_back(a);
	    }
	}
    }

#ifdef DRAW_HISTOS
  bool doPlot=draw_;// &&int(hest1->GetEntries())%30==40;
  doPlot=false;
  uint32_t seed=currentTime_;







  if (doPlot)
    {
  if (TCCluster==NULL)
    {

      TCCluster=new TCanvas("TCCluster","hough1",1000,900);

      TCCluster->Draw();
      TCCluster->Divide(3,2);
      TCCluster->Update();
      TCCluster->Modified();

    }
  hzx->SetLineColor(kRed);
  hzy->SetLineColor(kRed);
  czx->SetLineColor(kBlue);
  czy->SetLineColor(kBlue);
  czx->SetLineColor(kGreen);
  czy->SetLineColor(kGreen);
  mzx->SetMarkerColor(kBlack);
  mzy->SetMarkerColor(kBlack);
  mzx->SetMarkerStyle(21);
  mzy->SetMarkerStyle(21);

  TCCluster->cd(1);
  hzx->Draw("BOX");
  TCCluster->cd(2);
  czx->Draw("BOX");
  TCCluster->cd(3);
  izx->Draw("BOX");
  mzx->Draw("BOXSAME");
  TCCluster->cd(4);
  hzy->Draw("BOX");
  TCCluster->cd(5);
  czy->Draw("BOX");
  TCCluster->cd(6);
  izy->Draw("BOX");
  mzy->Draw("BOXSAME");
    }
#endif
  std::sort(theAmas_.rbegin(),theAmas_.rend());
  uint32_t nag=0;
  TH2F* amx[100];
  TH2F* amy[100];
  Components* c;
  uint32_t ifi=9999,ili=0;
  for (std::vector<Amas>::iterator ia=theAmas_.begin();ia!=theAmas_.end();ia++)
    {
      ia->compute();
      if (ia->getVolume()==0) continue;
      c=(Components*) ia->Components();
      // DEBUG_PRINT("%f %f %d %d \n",c->zmin,c->zmax,int(c->zmin/2.8),int(c->zmax/2.8));
      if (int(c->zmin/2.8)+1<ifi) ifi =int(c->zmin/2.8)+1;
      if (int(c->zmax/2.8)+1>ili) ili =int(c->zmax/2.8)+1;
#ifdef DRAW_HISTOS
      std::stringstream sn;
      sn<<"amasxz"<<nag;
      amx[nag]=new TH2F(sn.str().c_str(),sn.str().c_str(),120,0.,170.,120,-10.,110.);
      for (std::vector<RecoHit*>::iterator iha=ia->getHits().begin(); iha!=ia->getHits().end();iha++)
	{
	  amx[nag]->Fill((*iha)->Z(),(*iha)->X());
	}
      sn<<"y";
      amy[nag]=new TH2F(sn.str().c_str(),sn.str().c_str(),120,0.,170.,120,-10.,110.);
      for (std::vector<RecoHit*>::iterator iha=ia->getHits().begin(); iha!=ia->getHits().end();iha++)
	{
	  amy[nag]->Fill((*iha)->Z(),(*iha)->Y());
	}
      if (doPlot)
	{
      amx[nag]->SetLineColor(nag+1);
      amy[nag]->SetLineColor(nag+1);
      TCCluster->cd(3);
      if (nag==0)
	amx[nag]->Draw("BOX");
      else
	amx[nag]->Draw("BOXSAME");
      TCCluster->cd(6);
      if (nag==0)
	amy[nag]->Draw("BOX");
      else
	amy[nag]->Draw("BOXSAME");
	}
      nag++;
#endif
      // ia->Print();
    }

   DEBUG_PRINT("Amas %d MIP %d hits %d\n",nag,nmip,vrh.size());
  if (nag==0 && theComputerTrack_->Length()==0) 
    {
      for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
	 delete (*ic);
       free(h_x);
       free(h_y);
       free(h_z);
       free(h_layer);
       return nshower;
    }
#ifdef FILL_HISTOS
  isPion_=(lom>1E-2 && lom<0.5) && fpi<15 && theCerenkovTag_!=0 ;
  isProton_=(lom>1E-2 && lom<0.5) && fpi<15 && theCerenkovTag_==0 ;
  isElectron_=lom<1E-2;
  isMuon_=lom>0.5;
  
  isShower_=true;
  // printf("Ck tag = %d %d %d\n",theCerenkovTag_,isPion_,isProton_);
  //  if ((lom>1E-2 && lom<0.5) &&theLastRate_<520. && fpi>3 && fpi<15  )
  if (isProton_)
    {
      printf("Filling hproton\n");
      hproton->Fill(vrh.size());
    }
  if (isPion_)
    {
      hpion->Fill(vrh.size());
    }
  if (isPion_ || isProton_ )
    {
      hhadron->Fill(vrh.size());
      for (std::vector<RecoHit*>::iterator ih=vrh.begin();ih<vrh.end();ih++)
	{
	  bool thr[3];
	  thr[0]=(*ih)->getAmplitude()==2;
	  thr[1]=(*ih)->getAmplitude()==1;
	  thr[2]=(*ih)->getAmplitude()==3;
      
	  std::stringstream namech("");
	  namech<<"/Clusters/Pions/Chamber"<<(*ih)->chamber();

	  TH2* hthr0 = rootHandler_->GetTH2(namech.str()+"/Seuil0");
	  TH2* hthr1 = rootHandler_->GetTH2(namech.str()+"/Seuil1");
	  TH2* hthr2 = rootHandler_->GetTH2(namech.str()+"/Seuil2");
	  if (hthr0==NULL)
	    {
	      hthr0 =rootHandler_->BookTH2( namech.str()+"/Seuil0",96,0.,96.,96,0.,96.);
	      hthr1 =rootHandler_->BookTH2( namech.str()+"/Seuil1",96,0.,96.,96,0.,96.);
	      hthr2 =rootHandler_->BookTH2( namech.str()+"/Seuil2",96,0.,96.,96,0.,96.);
	    }
	  int chamberLocalI=(*ih)->I();
	  int chamberLocalJ=(*ih)->J();
	  if (thr[1]||thr[2]) hthr0->Fill(chamberLocalI*1.,chamberLocalJ*1.);
	  hthr1->Fill(chamberLocalI*1.,chamberLocalJ*1.);
	  if (thr[2]) hthr2->Fill(chamberLocalI*1.,chamberLocalJ*1.);
	}

    }
  if (isElectron_)
    helectron->Fill(vrh.size());
  if (isMuon_)
    hmuon->Fill(vrh.size());
#endif
#ifdef DRAW_HISTOS
  TLine* l[100];
  TLine* l1[100];
  uint32_t nline=0;
  if (doPlot)
    {
  TH1F* hweight=(TH1F*) rootHandler_->GetTH1("showerweight");
  TH1F* hweight2=(TH1F*) rootHandler_->GetTH1("showerweight2");
  
  TCCluster->cd(6);
  //mzy->Draw("BOX");
  //if (hest1!=0) hest1->Draw();
  // if (hweight2!=NULL) hweight2->Draw("text");


  
  
  for (unsigned int i=0;i<theComputerTrack_->getCandidates().size();i++)
	{
	  RecoCandTk& tk = theComputerTrack_->getCandidates()[i];
	  DEBUG_PRINT("%f %d \n",tk.chi2_,tk.getList().size());
	  TCCluster->cd(3);
	  l[i] = new TLine(tk.zmin_,tk.getXext(tk.zmin_),tk.zmax_,tk.getXext(tk.zmax_));
	  l[i]->SetLineColor(i+1);
	  l[i]->SetLineWidth(2);
	  l[i]->Draw("SAME");
	  TCCluster->cd(6);
	  l1[i] = new TLine(tk.zmin_,tk.getYext(tk.zmin_),tk.zmax_,tk.getYext(tk.zmax_));
	  l1[i]->SetLineColor(i+1);
	  l1[i]->SetLineWidth(2);
	  l1[i]->Draw("SAME");
	  nline++;
	}
  TCCluster->cd(3);
  mzx->SetMarkerSize(0.5);
  mzy->SetMarkerSize(0.5);
  mzx->Draw("SAME");
  TCCluster->cd(6);
  mzy->Draw("SAME");

    }
  if (hest!=NULL) hest->Fill(nmip*1./vrh.size());
  if (hest1!=NULL) hest1->Fill(1.*vrh.size());
  if (hrate!=NULL) hrate->Fill(theLastRate_,1.*vrh.size());
  if (hrate0!=NULL) hrate0->Fill(theLastRate_,1.*nh0);
  if (hrate1!=NULL) hrate1->Fill(theLastRate_,1.*nh1);
  if (hrate2!=NULL) hrate2->Fill(theLastRate_,1.*nh2);
  if (hspill!=NULL) hspill->Fill((theBCID_-currentTime_-theBCIDSpill_)*2E-7,theLastRate_);
  if (theLastRate_<120) 
    hest2->Fill(1.*vrh.size());
  else
    hest3->Fill(1.*vrh.size());
  theNbShowers_++;
  nshower++;
#endif

  //  hest1->Draw();
  DEBUG_PRINT("Hits %ld  Clusters (%ld,%.2f)  REAL (%d,%.2f,%.2f) LOM %f FP %d %d LP %d %d \n",vrh.size(),clusters.size(),clusters.size()*100./vrh.size(),realc.size(),realc.size()*170/vrh.size(),realc.size()*100./clusters.size(),lom,fpi,ifi,lpi,ili);
#ifdef DRAW_HISTOS
  if (doPlot)
    {
      TCCluster->Modified();
      TCCluster->Update();
      //TCCluster->WaitPrimitive();
      char ci;ci=getchar();putchar(ci); if (ci=='.') exit(0);
    }
#endif
 end:
   for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
     delete (*ic);
#ifdef DRAW_HISTOS
   
   for (unsigned int i=0;i<nline;i++)
	{
	  if (doPlot)
	    {
	      delete l[i];
	      delete l1[i];
	    }
	}
   for (uint32_t i=0;i<nag;i++)
     {
       delete amx[i];
       delete amy[i];

     }
#endif
   free(h_x);
   free(h_y);
   free(h_z);
   free(h_layer);
   return nshower;
}
