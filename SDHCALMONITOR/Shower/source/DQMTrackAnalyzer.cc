#define NPLANS_USED 6
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMRun.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMCoreTool.h"
#include "DQMTrackAnalyzer.h"
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

uint32_t DQMTrackAnalyzer::PMAnalysis(uint32_t bifid)
{

  theMonitoring_->cd("/BIFAnalysis"); 
  TH1* hpattag= theMonitoring_->getRealHistogram1D("PatternTagNoSeed","BIF Pattern no time seed ",100,0.1,100.1)->get<TH1F>();

  


  for (std::vector<DIFPtr*>::iterator itb = reader_->getDIFList().begin();itb!=reader_->getDIFList().end();itb++)
    {
      DIFPtr* d = (*itb);
      if (d->getID()!=bifid) continue;
      // Loop on frames

      for (uint32_t i=0;i<d->getNumberOfFrames();i++)
  	{
  	  //if (abs(d->getFrameTimeToTrigger(i)-seed)<200)
	  //printf(" Frame %d  Cerenkov time %d  \n",i,d->getFrameTimeToTrigger(i));
	  float ti=d->getFrameTimeToTrigger(i)*1.;
	  
	  std::bitset<64> chb(0);
	  for (std::vector<DIFPtr*>::iterator it = reader_->getDIFList().begin();it!=reader_->getDIFList().end();it++)
	    {
	      DIFPtr* dc = (*it);


	      std::map<unsigned int,DifGeom>::iterator idg = reader_->getDifMap().find(dc->getID());
	      DifGeom& difgeom = idg->second;
	      uint32_t chid = idg->second.getChamberId();
	      std::stringstream s;
	      s<<"BIFPOS"<<chid;

	      TH1* hpattag1=theMonitoring_->getRealHistogram1D(s.str(),s.str(),200,-30.,30.)->get<TH1F>();
	      for (uint32_t j=0;j<dc->getNumberOfFrames();j++)
		{
		  
		  float tj=dc->getFrameTimeToTrigger(j)*1.;
		  hpattag1->Fill(ti-tj); 
		  //if (chid==10)
		  // printf(" Slot10 Frame %d time %d %f %f %f\n",j,dc->getFrameTimeToTrigger(j),ti,tj,ti-tj);
		  if ((ti-tj)>=-4.5 && (ti-tj)<=-1.5)
		    {chb[chid]=1;}

		   if ((ti-tj)>=-0.5 && (ti-tj)<=2.5)
		    {chb[50+chid]=1;}
		}

	    }
	  
	  for(int i=0;i<64;i++)
	    {
	      if (chb[i]!=0) hpattag->Fill(i*1.);
	    }
	  hpattag->Fill(100.);

	  //if (chb[60]==0) {printf("Quoi qu on fout la \n ");getchar();}


	}

      


    }

 
  //  getchar();
  return 0;
}


uint32_t DQMTrackAnalyzer::CerenkovTagger(uint32_t difid,uint32_t seed)
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

void DQMTrackAnalyzer::initHistograms()
{
  //  rootHandler_->BookTH1("/Clusters/EST1",100,0.,300.);
}


DQMTrackAnalyzer::DQMTrackAnalyzer(DHCalEventReader* r,DQMModule* m) :trackIndex_(0),nAnalyzed_(0),clockSynchCut_(8), spillSize_(90000),maxHitCount_(500000),
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

void DQMTrackAnalyzer::initialise()
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
void DQMTrackAnalyzer::initJob(){presetParameters();}
void DQMTrackAnalyzer::endJob(){
 
	
}
void DQMTrackAnalyzer::presetParameters()
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
bool DQMTrackAnalyzer::decodeTrigger(LCCollection* rhcol, double tcut)
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


double  DQMTrackAnalyzer::checkTime()
{

  ftime(&theCurrentTime_);  
  double dt=theCurrentTime_.time-theTime_.time+(theCurrentTime_.millitm-theTime_.millitm)*1E-3;
  theTime_.time= theCurrentTime_.time;  theTime_.millitm= theCurrentTime_.millitm;
  return dt;
}
void DQMTrackAnalyzer::findTimeSeeds( IMPL::LCCollectionVec* rhcol, int32_t nhit_min,std::vector<uint32_t>& candidate)
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

void DQMTrackAnalyzer::processSeed(IMPL::LCCollectionVec* rhcol,uint32_t seed)
{

 

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

}

void DQMTrackAnalyzer::processEvent()
{


  checkTime();
  if (reader_->getEvent()==0) return;
  
  evt_=reader_->getEvent();
  //theSkip_=380;
  if (evt_->getEventNumber()<=theSkip_) return;
  printf("Processing %d - %d \n",evt_->getRunNumber(),evt_->getEventNumber());

  if (nAnalyzed_==0)
    {
      //      std::stringstream s;
      //      s<<"./Shower"<<evt_->getRunNumber()<<".root";
      //      this->createTrees(s.str());
    }
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
      reader_->findDIFSeeds(minChambersInTime_);
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

  // Make noise and time analysis
  theMonitoring_->FillTimeAsic(rhcol);

  //  LCTOOLS::printParameters(rhcol->parameters());
  //DEBUG_PRINT("Time Stamp %d \n",evt_->getTimeStamp());


  //DEBUG_PRINT("Calling decodeTrigger\n");
  // TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if (!decodeTrigger(rhcol,spillSize_) ) { if (rhcoltransient) delete rhcol;return;}
  
  //if (isNewSpill_) return;
  // if (evt_->getEventNumber()%100 ==0)
  //   rootHandler_->writeSQL();
  //    rootHandler_->writeXML(theMonitoringPath_);
  PMAnalysis(theBIFId_);

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


uint32_t DQMTrackAnalyzer::buildTracks(std::vector<RecoHit*> &vrh,std::string vdir)
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
  //ComputerTrack ch(&cuts);
  //ch.DefaultCuts();

  for (std::vector<RecoHit*>::iterator ih=vrh.begin();ih<vrh.end();ih++)
    {
      // DEBUG_PRINT("Hit plan = %d %d \n",(*ih)->chamber(),(*ih)->plan());
      //      if ((*ih)->getFlag(RecoHit::CORE)==1) continue;
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
  
  DEBUG_PRINT(" Number of clusters %d REALC %d INTC %d Hits %d \n",clusters.size(),realc.size(),intc.size(),vrh.size());

  uint32_t nstub=0;
  for (std::vector<RECOCluster*>::iterator ic=realc.begin();ic!=realc.end();ic++)
    {

      ChamberPos& cp=reader_->getPosition((*ic)->chamber());
	  // DEBUG_PRINT(" %d (%f,%f,%f) (%f,%f,%f) (%d,%d) \n",
	  //	 cp.getId(),cp.getX0(),cp.getY0(),cp.getZ0(),cp.getX1(),cp.getY1(),cp.getZ1(),cp.getXsize(),cp.getYsize());
      double x,y,z;
      cp.calculateGlobal((*ic)->X(),(*ic)->Y(),x,y,z);

      h_x[nstub]=x;//(*ic)->X();
      h_y[nstub]=y;//(*ic)->Y();
      h_z[nstub]=z;//(*ic)->Z();
      h_layer[nstub]=(*ic)->plan();

      // DEBUG_PRINT("\t %d :  %d %f %f %f \n",nstub,h_layer[nstub],h_x[nstub],h_y[nstub],h_z[nstub]);
      nstub++;
    }
  //  theComputerTrack_->associate(nstub,h_x,h_y,h_z,h_layer);
  theComputerTrack_->telescope(nstub,h_x,h_y,h_z,h_layer,theNplansUsed_);
  //theComputerTrack_->muonFinder(nstub,h_x,h_y,h_z,h_layer);
 

  //  if (theComputerTrack_->getTracks().size()>0) theNbTracks_++;
  uint32_t nmip=0;
   for (unsigned int i=0;i<theComputerTrack_->getTracks().size();i++)
	{
	  TrackInfo& tk = theComputerTrack_->getTracks()[i];

	  //if (tk.size()<minChambersInTime_) continue;
	  //if (fabs(tk.ax())<1.E-2) continue;
	  //if (fabs(tk.ax())<0.5 && fabs(tk.ay())<0.5) theNbTracks_++;
	  //this->draw(tk);
	  //char c;c=getchar();putchar(c); if (c=='.') exit(0);
	  uint32_t fch=int(ceil(tk.zmin()*10))/28+1;
	  uint32_t lch=int(ceil(tk.zmax()*10))/28+1;

	  std::stringstream st;
	  st<<vdir<<"/";
	  theMonitoring_->cd(st.str());
	  
	  TH1* hnp= theMonitoring_->getRealHistogram1D("Npoints","Number Of points per track",51,-0.1,50.9)->get<TH1F>();
	  TH1* hnpl= theMonitoring_->getRealHistogram1D("Nplanes","Number Of points per plane",51,-0.1,50.9)->get<TH1F>();
	  TH1* hax= theMonitoring_->getRealHistogram1D("ax","(x,z) coefficient directeur",200,-5.,5.)->get<TH1F>();
	  TH1* hay= theMonitoring_->getRealHistogram1D("ay","(y,z) coefficient directeur",200,-5.,5.)->get<TH1F>();

	  // Calcul de l'efficacite

	  // Track info
	  
	  hnp->Fill(tk.size()*1.);
	  hax->Fill(tk.ax());
	  hay->Fill(tk.ay());
	  fch=1;lch=theNplansUsed_;
	  for (int ip=fch;ip<=lch;ip++)
	    if (tk.plane(ip)) hnpl->Fill(ip*1.);
	  //	  std::cout<<tk.planes_<<std::endl;
	  //getchar();
	  for (uint32_t ip=fch;ip<=lch;ip++)
	    {
	      
	      TrackInfo tex;
	      
	      tk.exclude_layer(ip,tex);
	      uint32_t npext=tex.size();
	      /*
	      if (npext<minChambersInTime_) continue; // Au moins 4 plans dans l'estrapolation touches 

	      if (ip>1 && !tex.plane(ip-1)) 
		if (ip>2 && !tex.plane(ip-2)) continue;

	      if (ip<lch && !tex.plane(ip+1))
		if (ip<(lch-1) && !tex.plane(ip+2)) continue;
	      */
	      if (npext<3) continue;

	      std::stringstream s;
	      s<<st.str()<<"Plan"<<ip<<"/";

	      float dz0=0.,distz=60.; // 2.8
	      float xext=tex.xext(dz0+(ip-1)*distz);
	      float yext =tex.xext(dz0+(ip-1)*distz);
	       std::map<uint32_t,ChamberPos>& pos= reader_->getPositionMap();
	       double xi=1000,xa=-1000,yi=1000,ya=-1000;
	      for (std::map<uint32_t,ChamberPos>::iterator ich=pos.begin();ich!=pos.end();ich++)
		{
		  if ((*ich).second.getPlan()!=ip) continue;
		   xext=tex.xext((*ich).second.getZ0());
		   yext =tex.yext((*ich).second.getZ0());
		   if ((*ich).second.getX0()<xi) xi= (*ich).second.getX0();
		   if ((*ich).second.getY0()<yi) yi= (*ich).second.getY0();
		   if ((*ich).second.getX0()>xa) xa= (*ich).second.getX0();
		   if ((*ich).second.getY0()>ya) ya= (*ich).second.getY0();
		   if ((*ich).second.getX1()<xi) xi= (*ich).second.getX1();
		   if ((*ich).second.getY1()<yi) yi= (*ich).second.getY1();
		   if ((*ich).second.getX1()>xa) xa= (*ich).second.getX1();
		   if ((*ich).second.getY1()>ya) ya= (*ich).second.getY1();

		   break;
		}
	      int nx=int(xa-xi)+1;
	      int ny=int(ya-yi)+1;
	      theMonitoring_->cd(s.str());
	      TH2* hext=theMonitoring_->getRealHistogram2D("ext","Extrapolation of the track",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH2* hfound=theMonitoring_->getRealHistogram2D("found","Ext. of the track, cluster found",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH2* hfound1=theMonitoring_->getRealHistogram2D("found1","Ext. of the track, cluster (th2) found",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH2* hfound2=theMonitoring_->getRealHistogram2D("found2","Ext. of the track, cluster (th3) found",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH2* hnear=theMonitoring_->getRealHistogram2D("near","Position of cluster found",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH2* hmul=theMonitoring_->getRealHistogram2D("mul","Multiplicity of cluster found",nx,xi,xa,ny,yi,ya)->get<TH2F>();
	      TH1* hdx=theMonitoring_->getRealHistogram1D("dx","X(Ext)-X(Cluster)",400,-4.,4.)->get<TH1F>();
	      TH1* hdy=theMonitoring_->getRealHistogram1D("dy","Y(Ext)-Y(Cluster)",400,-4.,4.)->get<TH1F>();
	     
	      hext->Fill(xext,yext);
	      //bool 
	      float dist=1E9;
	      bool th1=false,th2=false;
	      float dxi,dyi,xn,yn,nhi;
	      for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
		{
		  if ((*ic)->plan()!=ip) continue;
		  ChamberPos& cp=reader_->getPosition((*ic)->chamber());
	  // DEBUG_PRINT(" %d (%f,%f,%f) (%f,%f,%f) (%d,%d) \n",
	  //	 cp.getId(),cp.getX0(),cp.getY0(),cp.getZ0(),cp.getX1(),cp.getY1(),cp.getZ1(),cp.getXsize(),cp.getYsize());
		  double x,y,z;
		  cp.calculateGlobal((*ic)->X(),(*ic)->Y(),x,y,z);
		  xext=tex.xext(z);
		  yext=tex.yext(z);
		  float dx=xext-x;
		  float dy=yext-y;

		  double dap=tex.closestApproach(x,y,z);
		  //  DEBUG_PRINT(" (%f,%f,%f) %f %f \n",x,y,z,dap,sqrt(dx*dx+dy*dy));
		  //getchar();
		  if (dap<dist)
		    {
		      
		      dist=dap;
		      dxi=dx;
		      dyi=dy;
		      nhi=(*ic)->size();
		      xn=x;
		      yn=y;
		      th1=false,th2=false;
		      for (std::vector<RecoHit>::iterator ih=(*ic)->getHits()->begin();ih!=(*ic)->getHits()->end();ih++)
			{
			  if ((*ih).getFlag(RecoHit::THR1)!=0) th1=true;
			  if ((*ih).getFlag(RecoHit::THR2)!=0) th2=true;
			}

		    }

		}
	      // Cut a 1.5 au lieu de 6
	      if (dist<2.5)
		{
		  hdx->Fill(dxi);
		  hdy->Fill(dyi);
		  hmul->Fill(xext,yext,nhi*1.);
		  hfound->Fill(xext,yext);
		  hnear->Fill(xn,yn);
		  if (th1||th2)  hfound1->Fill(xext,yext);
		  if (th2)  hfound2->Fill(xext,yext);
		}
	    }





	}
   DEBUG_PRINT("==> MIPS hit %d -> %.2f\n",nmip,nmip*100./vrh.size()); 
 
 
   for (std::vector<RECOCluster*>::iterator ic=clusters.begin();ic!=clusters.end();ic++)
     delete (*ic);
   free(h_x);
   free(h_y);
   free(h_z);
   free(h_layer);
   return nshower;
}
