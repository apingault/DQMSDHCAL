#include "DQMSDHCALMonitor.h"
#include "dqm4hep/DQMMonitorElement.h"
#include "dqm4hep/DQMRun.h"
#include "dqm4hep/DQMXmlHelper.h"
#include "dqm4hep/DQMModuleApi.h"
#include "dqm4hep/DQMCoreTool.h"

DQMSDHCALMonitor::DQMSDHCALMonitor(DQMModule* m) : theModule_(m),theTrackIndex_(0),theFirstChamber_(1),theLastChamber_(50),
						   theExtrapolationMinimumPoint_(6),theExtrapolationMinimumChi2_(0.001),theExtrapolationDistanceCut_(8),theTrackAngularCut_(0.1),theChamberEdgeCut_(5.),useTk4_(true),theIntegratedTime_(0)
{
  theAsicCount_.clear();
  DHCalEventReader* reader=DHCalEventReader::instance();

  theDifMap_=reader->getDifMap();
  theChamberMap_=reader->getChamberMap();
}
void DQMSDHCALMonitor::clear()
{
  theAsicCount_.clear();
  theTrackIndex_=0;
  theIntegratedTime_=0;
}
void DQMSDHCALMonitor::setFirstChamber(uint32_t i) 
{theFirstChamber_=i;}

void DQMSDHCALMonitor::setLastChamber(uint32_t i)
{theLastChamber_=i;}
void DQMSDHCALMonitor::setExtrapolationMinimumPoint(uint32_t i) 
{theExtrapolationMinimumPoint_=i;}
void DQMSDHCALMonitor::setExtrapolationMinimumChi2(float i) 
{theExtrapolationMinimumChi2_=i;}
void DQMSDHCALMonitor::setChamberEdgeCut( float i) 
{theChamberEdgeCut_=i;} 
void DQMSDHCALMonitor::setUseTk4(bool t)
{useTk4_=t;};
	
void DQMSDHCALMonitor::FillTimeAsic(IMPL::LCCollectionVec* rhcol)
{
  //  std::map<uint32_t,uint32_t> count;
  //count.clear();
  dqm4hep::DQMModuleApi::cd(theModule_);
  StatusCode rc=dqm4hep::DQMModuleApi::cd(theModule_,"/Timing");
  if (rc!=STATUS_CODE_SUCCESS)
    {
      dqm4hep::DQMModuleApi::mkdir(theModule_,"/Timing");
      dqm4hep::DQMModuleApi::cd(theModule_,"/Timing");
    }
  DQMMonitorElement *moccall=NULL,*moccalldif=NULL,*moccallchamber=NULL,*masic2=NULL,*mtdif=NULL;

  if (dqm4hep::DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancy",moccall)!=STATUS_CODE_SUCCESS)
    {
      DQMModuleApi::bookRealHistogram1D(theModule_,moccall,"AsicOccupancy", "Occupancy per Asic",255*48,0.,255*48.); 
      DQMModuleApi::bookRealHistogram1D(theModule_,moccalldif,"AsicOccupancyDIF", "Occupancy per DIF",255,0.,255); 
      DQMModuleApi::bookRealHistogram1D(theModule_,moccallchamber,"AsicOccupancyChamber", "Occupancy per chamber",255,0.,255); 
      DQMModuleApi::bookRealHistogram2D(theModule_,masic2,"DIFAsicCount", "Occupancy ASIC vs DIF",256,0.1,256.1,48,0.1,48.1);
       DQMModuleApi::bookRealHistogram1D(theModule_,mtdif,"TimeDif", " time integrated in DIF",2000,0.,4000000.);
    }
  DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancy",moccall);
  DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancyDIF",moccalldif);
  DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancyChamber",moccallchamber);
  DQMModuleApi::getMonitorElement(theModule_,"DIFAsicCount",masic2);
  DQMModuleApi::getMonitorElement(theModule_,"TimeDif",mtdif);


  TH1* hoccall= moccall->get<TH1F>();
  TH1* hoccalldif= moccalldif->get<TH1F>();
  TH1* htdif= mtdif->get<TH1F>();
  TH1* hoccallchamber= moccallchamber->get<TH1F>();
  TH2* hasic2= masic2->get<TH2F>();
  
  hoccalldif->Reset();
  hoccallchamber->Reset();
  hasic2->Reset();
  double tmin=99999999.;
  double tmax=0.;
  //	printf("%d %s\n",__LINE__,__PRETTY_FUNCTION__);
  this->DIFStudy(rhcol);
  //	printf("%d %s\n",__LINE__,__PRETTY_FUNCTION__);
  //IMPL::LCCollectionVec* rhcol=(IMPL::LCCollectionVec*) evt_->getCollection(collectionName_);
  for (int i=0;i<rhcol->getNumberOfElements();i++)
    {
      IMPL::RawCalorimeterHitImpl* hit = (IMPL::RawCalorimeterHitImpl*) rhcol->getElementAt(i);
      if (hit==0) continue;
      
      // Decode
      unsigned int difid = hit->getCellID0()&0xFF;
      int asicid = (hit->getCellID0()&0xFF00)>>8;
      //hasic2->Fill(difid*1.,asicid*1);
      int channel= (hit->getCellID0()&0x3F0000)>>16;
      unsigned int bc = hit->getTimeStamp();
      if (bc>5E6) continue;
      if (bc<tmin) tmin=bc;
      if (bc>tmax) tmax=bc;
      
      std::map<unsigned int,DifGeom>::iterator idg = theDifMap_.find(difid);
      DifGeom& difgeom = idg->second;
      uint32_t chid = idg->second.getChamberId();
      //if (chid==0) printf("L'id de la chambre st %d %d\n",chid,difid);
      uint32_t key=(chid<<16)|(difid<<8)|asicid;
      std::map<uint32_t,uint32_t>::iterator it=theAsicCount_.find(key);
      if (theAsicCount_.find(key)!=theAsicCount_.end()) 
	it->second=it->second+1;
      else
	{
	  uint32_t n=1;
	  std::pair<uint32_t,uint32_t> p(key,n);
	  theAsicCount_.insert(p);
	}
    }
  //	printf("%d %s\n",__LINE__,__PRETTY_FUNCTION__);
  theIntegratedTime_+=(tmax-tmin);
  theEventIntegratedTime_=(tmax-tmin);
  //std::cout<<tmin<<" "<<tmax<<" => Event "<<theEventIntegratedTime_<<" total " <<theIntegratedTime_<<std::endl;
  for (std::map<uint32_t,uint32_t>::iterator it=theAsicCount_.begin();it!=theAsicCount_.end();it++)
    {
      uint32_t chid =(it->first>>16)&0xFF;
      uint32_t difid =(it->first>>8)&0xFF;
      uint32_t asicid =(it->first)&0xFF;
      
      std::stringstream namec("");
      namec<<"/Noise/Chamber"<<chid<<"/DIF"<<difid;
      this->cd(namec.str());		

      DQMMonitorElement *mocc=NULL,*moccn=NULL;

      if (dqm4hep::DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancy",mocc)!=STATUS_CODE_SUCCESS)
	{
	  DQMModuleApi::bookRealHistogram1D(theModule_,mocc,"AsicOccupancy", "Occupancy",48,0.,48.); 
	  DQMModuleApi::bookRealHistogram1D(theModule_,moccn,"AsicOccupancyNumber", "Occupancy count",48,0.,48.); 
	}
      DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancy",mocc);
      DQMModuleApi::getMonitorElement(theModule_,"AsicOccupancyNumber",moccn);


      TH1* hocc= mocc->get<TH1F>();
      TH1* hoccn= moccn->get<TH1F>();
		
      hoccn->SetBinContent(asicid,it->second);
      hoccall->SetBinContent(difid*48+asicid,it->second/(theIntegratedTime_*DCBufferReader::getDAQ_BC_Period()*1.E-6));
		
      hocc->SetBinContent(asicid,it->second/(theIntegratedTime_*DCBufferReader::getDAQ_BC_Period()*1.E-6));
      float focc=it->second/(theIntegratedTime_*DCBufferReader::getDAQ_BC_Period()*1.E-6);
      if (focc>hoccallchamber->GetBinContent(chid)) hoccallchamber->SetBinContent(chid,focc);
      if (focc>hoccalldif->GetBinContent(difid)) hoccalldif->SetBinContent(difid,focc);


    }
  //	printf("%d %s\n",__LINE__,__PRETTY_FUNCTION__);


  htdif->Fill((tmax-tmin)*1.);
}
void DQMSDHCALMonitor::cd(std::string name)
{
  //printf("Changing directory to %s\n",name.c_str());
  dqm4hep::DQMModuleApi::cd(theModule_);
  StatusCode rc=dqm4hep::DQMModuleApi::cd(theModule_,name);
  if (rc!=STATUS_CODE_SUCCESS)
    {
      dqm4hep::DQMModuleApi::mkdir(theModule_,name);
      dqm4hep::DQMModuleApi::cd(theModule_,name);
    }

}
DQMMonitorElement* DQMSDHCALMonitor::getRealHistogram1D(std::string name,std::string title,uint32_t nx,float xi,float xa)
{
  DQMMonitorElement* m=NULL;
  if (dqm4hep::DQMModuleApi::getMonitorElement(theModule_,name,m)!=STATUS_CODE_SUCCESS)
    {
      DQMModuleApi::bookRealHistogram1D(theModule_,m,name,title,nx,xi,xa);

    }

  return m;
}
DQMMonitorElement* DQMSDHCALMonitor::getRealHistogram2D(std::string name,std::string title,uint32_t nx,float xi,float xa,float ny,float yi,float ya)
{
  DQMMonitorElement* m=NULL;
  if (dqm4hep::DQMModuleApi::getMonitorElement(theModule_,name,m)!=STATUS_CODE_SUCCESS)
    {
      DQMModuleApi::bookRealHistogram2D(theModule_,m,name,title,nx,xi,xa,ny,yi,ya);

    }

  return m;
}
void DQMSDHCALMonitor::DIFStudy( IMPL::LCCollectionVec* rhcol)
{
    
  for (int i=0;i<rhcol->getNumberOfElements();i++)
    {
      IMPL::RawCalorimeterHitImpl* hit = (IMPL::RawCalorimeterHitImpl*) rhcol->getElementAt(i);
      if (hit==0) continue;

      unsigned int difid = hit->getCellID0()&0xFF;
      std::map<unsigned int,DifGeom>::iterator idg = theDifMap_.find(difid);
      DifGeom& difgeom = idg->second;
      uint32_t chid = idg->second.getChamberId();
      int asicid = (hit->getCellID0()&0xFF00)>>8;
      int channel= (hit->getCellID0()&0x3F0000)>>16;
      //streamlog_out(MESSAGE)<<"ch-"<<channel<<std::endl;
      unsigned int bc = hit->getTimeStamp();

      bool thr[3];
      //      DEBUG_PRINT("%x \n",hit->getCellID0());
      int ithr= hit->getAmplitude()&0x3;
      if (ithr<=0 || ithr>3)
	{
	  std::cout<<difid<<" had:"<<asicid<<":"<<channel<<":"<<bc<<":"<<ithr<<std::endl;
	  return;
	}
      thr[0] = (ithr == 1);
      thr[1] = (ithr == 2);
      thr[2] = (ithr == 3);
      int asic=asicid;int x,y;
      if (difid>1000) asic=(asic-1)%4+1; // Small chamber
      if (chid<49)
	DifGeom::PadConvert(asic,channel,x,y,2);
      else
	DifGeom::PadConvert(asic,channel,x,y,11);
      int difLocalI=int(x);
      int difLocalJ=int(y);
      int chamberLocalI=difgeom.toGlobalX(difLocalI);
      int chamberLocalJ=difgeom.toGlobalY(difLocalJ);

      //INFO_PRINT("%d %d %d %d %d %d \n",x,y,difLocalI,difLocalJ,chamberLocalI,chamberLocalJ);
      std::stringstream namec("");
      namec<<"/Noise/Chamber"<<chid<<"/DIF"<<difid;

      this->cd(namec.str());

      TH1* hhits0 = this->getRealHistogram1D("Hits0","Hits per pad and asic second seuil",48*64,0.1,48*64+0.1)->get<TH1F>();
      TH1* hhits1 = this->getRealHistogram1D("Hits1","Hits per pad and asic premier seuil",48*64,0.1,48*64+0.1)->get<TH1F>();
      TH1* hhits2 = this->getRealHistogram1D("Hits2","Hits per pad and asic troisieme seuil",48*64,0.1,48*64+0.1)->get<TH1F>();
      TH1* hetd = this->getRealHistogram1D("EventTime","Frame time",10000,0.,1.5E7)->get<TH1F>();
      TH1* hetdz = this->getRealHistogram1D("EventTimeZoom","Frame time (zoomed)",10000,0.,10000)->get<TH1F>();

      if (thr[0]||thr[2]) hhits0->SetBinContent((asic-1)*64+channel+1,hhits0->GetBinContent((asic-1)*64+channel+1)+1);
      if (thr[1]||thr[0]||thr[2]) hhits1->SetBinContent((asic-1)*64+channel+1,hhits1->GetBinContent((asic-1)*64+channel+1)+1);
      if (thr[2]) hhits2->SetBinContent((asic-1)*64+channel+1,hhits2->GetBinContent((asic-1)*64+channel+1)+1);
      hetd->Fill(bc*1.);

      hetdz->Fill(bc*1.);

      std::stringstream namech("");
      namech<<"/Noise/Chamber"<<chid;
      this->cd(namech.str());
      TH2* hthr0 =this->getRealHistogram2D("Seuil0","Hit Map Deuxieme Seuil",96,0.,96.,96,0.,96.)->get<TH2F>();
      TH2* hthr1 =this->getRealHistogram2D("Seuil1","Hit Map Premier Seuil",96,0.,96.,96,0.,96.)->get<TH2F>();
      TH2* hthr2 =this->getRealHistogram2D("Seuil2","Hit Map Troisieme Seuil",96,0.,96.,96,0.,96.)->get<TH2F>();
      
      if (thr[0]||thr[2]) hthr0->Fill(chamberLocalI*1.,chamberLocalJ*1.);
      if (thr[1]||thr[2]||thr[0]) hthr1->Fill(chamberLocalI*1.,chamberLocalJ*1.);
      if (thr[2]) hthr2->Fill(chamberLocalI*1.,chamberLocalJ*1.);

    }
}

void DQMSDHCALMonitor::trackHistos(std::vector<RecoCandTk> &tracks,std::vector<RecoPoint> &points,std::string tkdir)
{
#ifdef UNUSED
  /*
    TH1* htimedif = rootHandler_->GetTH1("TimeToTrigger");
    if (htimedif==NULL)
    {
    htimedif=rootHandler_->BookTH1( "TimeToTrigger",5000,0.,5000.);

    }
    htimedif->Fill(currentTime_*1.);
    //std::string tkdir="/Tracking";
    if (currentTime_>=(uint32_t)clockSynchCut_) tkdir="/OtherTracking";
  */
  for (unsigned int i=0;i<points.size();i++)
    {	 

      //if (points[i].getCluster().getHits()->size()>1) continue;
      std::stringstream namec("");


      namec<<tkdir+"/Plan"<<points[i].getChamberId();
      if (points[i].isUsed()) 
	namec<<"/OnTrack";
      else
	namec<<"/OffTrack";
      TH2* hpos = rootHandler_->GetTH2(namec.str()+"/XYPos");	   
      TH2* hcpos = rootHandler_->GetTH2(namec.str()+"/XYClusterPos");	   
      TH2* hposhit = rootHandler_->GetTH2(namec.str()+"/XYPosHit");	   
      TH1* hposmul = rootHandler_->GetTH1(namec.str()+"/Multiplicity");	   
      if (hpos==NULL)
	{
	  hpos=rootHandler_->BookTH2( namec.str()+"/XYPos",115,-10.1,110.1,115,-10.1,110.1);
	  hcpos=rootHandler_->BookTH2( namec.str()+"/XYClusterPos",100,0.,96.,100,0.,96.);
	  hposhit=rootHandler_->BookTH2( namec.str()+"/XYPosHit",96,0.,96.,96,0.,96.);
	  hposmul=rootHandler_->BookTH1( namec.str()+"/Multiplicity",50,0.,50.);
	}
      hpos->Fill(points[i].X(),points[i].Y());
      hcpos->Fill(points[i].getCluster().X(),points[i].getCluster().Y());
      hposmul->Fill(points[i].getCluster().getHits()->size()*1.);
      for (std::vector<RecoHit>::iterator ih=points[i].getCluster().getHits()->begin();ih!=points[i].getCluster().getHits()->end();ih++)
	hposhit->Fill(ih->X()*1.,ih->Y()*1.);	
    }

  TH1* hngood = rootHandler_->GetTH1(tkdir+"/NumberOfTracks");
  if (hngood==0)
    {
      hngood = rootHandler_->BookTH1(tkdir+"/NumberOfTracks",21,-0.1,20.9);
    }
  hngood->Fill(tracks.size()*1.);
  if (tracks.size()==0) return;
  for (unsigned int itk=0;itk<tracks.size();itk++)
    {
      theTrackIndex_++;
      RecoCandTk& tk = tracks[itk];
      //DEBUG_PRINT("Tk=%d Time=%d %f %f %f %f %f \n",theTrackIndex_,currentTime_,tk.ax_,tk.bx_,tk.ay_,tk.by_,tk.chi2_);
      TH1* htchi2 = rootHandler_->GetTH1(tkdir+"/Chi2");
      TH1* htpchi2 = rootHandler_->GetTH1(tkdir+"/ProbChi2");
      TH1* htnpoint = rootHandler_->GetTH1(tkdir+"/NumberOfPoints");
      TH1* htax = rootHandler_->GetTH1(tkdir+"/Ax");
      TH1* htay = rootHandler_->GetTH1(tkdir+"/Ay");
      TH1* htxh = rootHandler_->GetTH1(tkdir+"/xh");
      TH1* htyh = rootHandler_->GetTH1(tkdir+"/yh");

      if (htchi2==0)
	{
	  htchi2 = rootHandler_->BookTH1(tkdir+"/Chi2",500,0.,100.);
	  htpchi2 = rootHandler_->BookTH1(tkdir+"/ProbChi2",1000,0.,1.);
	  htnpoint = rootHandler_->BookTH1(tkdir+"/NumberOfPoints",60,0.,60.);
	  htax = rootHandler_->BookTH1(tkdir+"/Ax",200,-10.,10.);
	  htay = rootHandler_->BookTH1(tkdir+"/Ay",200,-10.,10.);
	  htxh = rootHandler_->BookTH1(tkdir+"/xh",6000,0.,6000.);
	  htyh = rootHandler_->BookTH1(tkdir+"/yh",6000,0.,6000.);

	}

      htchi2->Fill(tk.chi2_/(2*tk.getList().size()-4.));
      htpchi2->Fill(tk.prChi2_);
      //DEBUG_PRINT("track %f %d %f %f \n",tk.chi2_,2*tk.getList().size()-4,TMath::Prob(tk.chi2_,2*tk.getList().size()-4),tkChi2Cut_);
      //getchar();
      htnpoint->Fill(tk.getList().size()*1.);

      htax->Fill(tk.ax_);

      htay->Fill(tk.ay_);


      for (unsigned int ich=0;ich<61;ich++)
	{
	  uint32_t bintk=((theTrackIndex_-1)%100)*60+ich+1;
	  htxh->SetBinContent(bintk,0);
	  htyh->SetBinContent(bintk,0);
	}

      //       if (tracks.size()>1)
      //DEBUG_PRINT("\t %d good hits found \n",tk.getList().size());
      for (unsigned int i =0;i<tk.getList().size();i++)
	{
	  // if (tracks.size()>1)
	  //   DEBUG_PRINT("\t \t %f %f %f \n",tk.getList()[i].X(),tk.getList()[i].Y(),tk.getList()[i].Z());
	  //(tk.getList())[i].Print();
	  std::stringstream namec("");
	  namec<<tkdir+"/Plan"<<(tk.getList())[i]->getChamberId();
	  TH1* hposx = rootHandler_->GetTH1(namec.str()+"/XPos");	   
	  TH1* hposy = rootHandler_->GetTH1(namec.str()+"/YPos");	   
	  TH1* hpullx = rootHandler_->GetTH1(namec.str()+"/XDist");	   
	  TH1* hpully = rootHandler_->GetTH1(namec.str()+"/YDist");	   
	  TH1* hmult = rootHandler_->GetTH1(namec.str()+"/Multiplicity");	   

	  if (hposx==0)
	    {
				
	      hposx =rootHandler_->BookTH1( namec.str()+"/XPos",115,-10.,110.);
	      hposy =rootHandler_->BookTH1( namec.str()+"/YPos",115,-10.,110.);
	      hpullx =rootHandler_->BookTH1( namec.str()+"/XDist",200,-5.,5.);
	      hpully =rootHandler_->BookTH1( namec.str()+"/YDist",200,-5.,5.);
	      hmult =rootHandler_->BookTH1( namec.str()+"/Multiplicity",50,0.,50.);
				
	    }
	  uint32_t bintk=((theTrackIndex_-1)%100)*60+(tk.getList())[i]->getChamberId();
	  htxh->SetBinContent(bintk,(tk.getList())[i]->X());
	  htyh->SetBinContent(bintk,(tk.getList())[i]->Y());

	  hposx->Fill((tk.getList())[i]->X());
	  hposy->Fill((tk.getList())[i]->Y());
	  hpullx->Fill((tk.getList())[i]->X() -tk.getXext((tk.getList())[i]->Z()) );
	  hpully->Fill((tk.getList())[i]->Y() -tk.getYext((tk.getList())[i]->Z()) );
	  hmult->Fill(tk.getList()[i]->getCluster().getHits()->size()*1.);
	}


      if (fabs(tk.ax_)>theTrackAngularCut_ || fabs(tk.ay_)>theTrackAngularCut_) continue;
      if (tracks.size()!=1) continue;

      std::bitset<255> intrack;
      bool synch= false;
      unsigned int s_shift=0;
      if (synch) 
	intrack.set(s_shift,true);
      else
	{
	  s_shift =100;
	  intrack.set(s_shift,true);

	}
      for (std::map<unsigned int,ChamberGeom>::iterator ip=theChamberMap_.begin();ip!=theChamberMap_.end();ip++)
	{
	  // std::cout<<ip->first<<" "<<ip->second.getZ()<<std::endl;
	  //getchar();
	  int chid = ip->first;
	  int32_t interval=3;
	  if (chid<=theFirstChamber_+1) interval=5;
	  if (chid>=theLastChamber_-1) interval=5;

	  int32_t tkFirstEx=((chid-interval)>theFirstChamber_)?(chid-interval):theFirstChamber_;
	  int32_t tkLastEx=((chid+interval)<theLastChamber_)?(chid+interval):theLastChamber_;
	  //std::cout<<chid<<" "<<tkFirstEx<<" "<<tkLastEx<<std::endl;

	  RecoCandTk tk0;
	  for (unsigned int j =0;j<tk.getList().size();j++)
	    {
	      //std::cout<<fabs(tk.getList()[j].Z()-zplane[ip])<<std::endl;
	      if (tk.getList()[j]->getChamberId()==chid) continue;
				
	      tk0.addPoint(*tk.getList()[j]);
	    }
	  //std::cout<<" extra "<<chid<<" " <<tk0.getList().size()<<std::endl;
	  float xext=-999999,yext=-999999;
	  if (useTk4_ && chid>=theFirstChamber_ && chid<=theLastChamber_)
	    {
	      RecoCandTk tk4;
	      for (unsigned int j =0;j<tk.getList().size();j++)
		{
		  //std::cout<<fabs(tk.getList()[j].Z()-zplane[ip])<<std::endl;
		  if (tk.getList()[j]->getChamberId()==chid) continue;
		  if (tk.getList()[j]->getChamberId()<tkFirstEx) continue;
		  if (tk.getList()[j]->getChamberId()>tkLastEx) continue;
		  tk4.addPoint(*tk.getList()[j]);

		  // if (chid==theLastChamber_)
		  //   if (tk.getList()[j]->getChamberId()==(chid-3)) tk4.addPoint(*tk.getList()[j]);
		  // if (chid>=theFirstChamber_+2)
		  //   if (tk.getList()[j]->getChamberId()==(chid-2)) tk4.addPoint(*tk.getList()[j]);
		  // if (chid>=theFirstChamber_+1)
		  // if (tk.getList()[j]->getChamberId()==(chid-1)) tk4.addPoint(*tk.getList()[j]);

		  // if (chid<=theLastChamber_-1)
		  //   if (tk.getList()[j]->getChamberId()==(chid+1)) tk4.addPoint(*tk.getList()[j]);
		  // if (chid<=theLastChamber_-2)
		  //   if (tk.getList()[j]->getChamberId()==(chid+2)) tk4.addPoint(*tk.getList()[j]);
		  // if (chid==theFirstChamber_)
		  //   if (tk.getList()[j]->getChamberId()==(chid+3)) tk4.addPoint(*tk.getList()[j]);
					

		}
	      tk4.regression();
	      if (tk4.getList().size()>=3 && TMath::Prob(tk4.chi2_,2*tk4.getList().size()-4)>theExtrapolationMinimumChi2_)
		{
		  xext=tk4.getXext(ip->second.getZ());
		  yext=tk4.getYext(ip->second.getZ());
		}
	      else 
		continue;

	    }
	  if (xext==-999999 || yext==-999999)
	    {
	      tk0.regression();
	      if (tk0.getList().size()<(uint32_t) theExtrapolationMinimumPoint_ || TMath::Prob(tk0.chi2_,2*tk0.getList().size()-4)<theExtrapolationMinimumChi2_) continue;
	      xext=tk0.getXext(ip->second.getZ());
	      yext=tk0.getYext(ip->second.getZ());
	    }
	  double xchext = ip->second.toLocalX(xext);
	  double ychext = ip->second.toLocalY(yext);
	  double zch=ip->second.getZ();
	  //std::cout <<xext<<" " <<yext<<" "<<xchext<<" " <<ychext<< " "<<zch<<" "<<ip->second.getZ()<<std::endl;
	  ip->second.calculateLocal(xext,yext,0,xchext,ychext,zch);
	  //std::cout <<xext<<" " <<yext<<" "<<xchext<<" " <<ychext<< " "<<zch<<" "<<ip->second.getZ()<<std::endl;
	  //	  getchar();
	  if (xchext<theChamberEdgeCut_ || xchext >96-theChamberEdgeCut_) continue;
	  if (ychext<theChamberEdgeCut_ || ychext >96-theChamberEdgeCut_) continue;
	  if (s_shift+ip->first<255)
	    intrack.set(s_shift+ip->first,true);
	  std::stringstream namec("");
	  namec<<tkdir+"/Plan"<<ip->first;

	  TH2* hextpos = rootHandler_->GetTH2( namec.str()+"/LocalExtrapolationMap");
	  TH2* hfoundpos = rootHandler_->GetTH2( namec.str()+"/LocalFoundMap");
	  TH2* hnearpos = rootHandler_->GetTH2( namec.str()+"/LocalNearestMap");
	  TH2* hmispos = rootHandler_->GetTH2( namec.str()+"/LocalMissedMap");
	  if (hextpos == NULL)
	    {
	      hextpos =rootHandler_->BookTH2( namec.str()+"/LocalExtrapolationMap",96,0.,96.,96,0.,96.);
	      hfoundpos =rootHandler_->BookTH2( namec.str()+"/LocalFoundMap",96,0.,96.,96,0.,96.);
	      hnearpos =rootHandler_->BookTH2( namec.str()+"/LocalNearestMap",96,0.,96.,96,0.,96.);
	      hmispos =rootHandler_->BookTH2( namec.str()+"/LocalMissedMap",96,0.,96.,96,0.,96.);

	    }
	  hextpos->Fill(xchext,ychext);
	  //std::cout<<xext<<" "<<yext<<std::endl;
	  unsigned int imin=999999;double distmin=9999999;
	  for (unsigned int irp=0;irp<points.size();irp++)
	    {
	      if (points[irp].getChamberId()!=chid) continue;
	      double dist =sqrt((points[irp].X()-xext)*(points[irp].X()-xext)+(points[irp].Y()-yext)*(points[irp].Y()-yext));
	      //std::cout<<chid<<" "<<dist<<std::endl;
	      if (dist<distmin)
		{
		  distmin=dist;
		  imin=irp;
		}
	    }
	  if (imin == 999999)
	    {
	      hmispos->Fill(xchext,ychext);
	      //		std::cout<<" chamber "<<chid<<" not found ("<<xchext<<","<<ychext<<") "<<tk0.ax_<<":"<<tk0.ay_<<":"<<tk0.chi2_<<std::endl;
	    }
	  // else
	  //   std::cout<<chid<<" is found "<<std::endl;
	  if (distmin<theExtrapolationDistanceCut_*1000)
	    {
	      std::stringstream namec("");
	      namec<<tkdir+"/Plan"<<ip->first;
	      TH1* hresol = rootHandler_->GetTH1(namec.str()+"/Resolution");	  
	      TH1* hresx = rootHandler_->GetTH1(namec.str()+"/XPull");	  
	      TH1* hresy = rootHandler_->GetTH1(namec.str()+"/YPull");	  
	      if (hresol==0)
		{
		  hresol =rootHandler_->BookTH1(namec.str()+"/Resolution",200,0.,50.);
		  hresx =rootHandler_->BookTH1( namec.str()+"/XPull",2000,-150.,150.);
		  hresy =rootHandler_->BookTH1( namec.str()+"/YPull",2000,-150.,150.);
		}
	      hresol->Fill(distmin);
	      hresx->Fill(points[imin].X()-xext);
	      hresy->Fill(points[imin].Y()-yext);
				

	      if (distmin<theExtrapolationDistanceCut_)
		{
		  if (s_shift+ip->first+60<255)
		    intrack.set(s_shift+ip->first+60,true);
		  double xchnear = ip->second.toLocalX(points[imin].X());
		  double ychnear = ip->second.toLocalY(points[imin].Y());
		  double zch;
		  ip->second.calculateLocal(points[imin].X(),points[imin].Y(),0,xchnear,ychnear,zch);

		  hnearpos->Fill(xchnear,ychnear);
		  hfoundpos->Fill(xchext,ychext);
					
		}
	      else
		if (1<0)
		  DEBUG_PRINT("Dist = %f X =%f Y =%f \n",distmin,points[imin].X()-xext,points[imin].Y()-yext);
	    }
	}
      TH1* hintrack= rootHandler_->GetTH1("PlanInTrack");
      if (hintrack==NULL)
	{
	  hintrack =rootHandler_->BookTH1( "PlanInTrack",255,-0.1,254.9);
	}
      for (unsigned int ib=0;ib<255;ib++)
	if (intrack[ib]!=0) hintrack->Fill(ib*1.);
		
    }
#endif
}

