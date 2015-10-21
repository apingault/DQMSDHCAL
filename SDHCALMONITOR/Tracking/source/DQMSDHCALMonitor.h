#ifndef _SDHCALMONITOR_H

#define _SDHCALMONITOR_H
#include <limits.h>
#include <vector>
#include "DHCalEventReader.h"
#include "dqm4hep/module/DQMModule.h"
using namespace dqm4hep;

#include <iostream>

#include "IMPL/LCTOOLS.h"
#include "EVENT/RawCalorimeterHit.h" 
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"
#include "UTIL/CellIDDecoder.h"
#include "DifGeom.h"
#include "ChamberGeom.h"


#include <vector>
#include <map>
#include "UtilDefs.h"
#include "Array3D.h"
#include "RecoHit.h"
#include "RECOCluster.h"
#include "RecoPoint.h"
#include "RecoCandTk.h"
#include "HC.h"
#include "HTImage.h"
#include "Shower.h"
#include "Amas.h"

using namespace dqm4hep;

class DQMSDHCALMonitor 
{
public:
	DQMSDHCALMonitor(DQMModule* h);
	void clear();
	void FillTimeAsic(IMPL::LCCollectionVec* rhcol);
	void DIFStudy( IMPL::LCCollectionVec* rhcol);	
	void trackHistos(std::vector<RecoCandTk> &tracks,std::vector<RecoPoint> &points,std::string tkdir="/OtherTracking");
	void cd(std::string name);
	DQMMonitorElement* getRealHistogram1D(std::string name,std::string title,uint32_t nx,float xi,float xa);
	DQMMonitorElement* getRealHistogram2D(std::string name,std::string title,uint32_t nx,float xi,float xa,float ny,float yi,float ya);

	void setFirstChamber(uint32_t i);

	void setLastChamber(uint32_t i);
	void setExtrapolationMinimumPoint(uint32_t i);
	void setExtrapolationMinimumChi2(float i);
	void setChamberEdgeCut( float i);
	void setUseTk4(bool t);
	int getEventIntegratedTime(){return  theEventIntegratedTime_;}
private:

	uint32_t theTrackIndex_,theFirstChamber_,theLastChamber_,theExtrapolationMinimumPoint_;
	float theExtrapolationMinimumChi2_,theChamberEdgeCut_ ,theTrackAngularCut_,theExtrapolationDistanceCut_;
	bool useTk4_;
	DQMModule* theModule_;
	// Reader
	//DHCalEventReader* reader_;

	std::map<uint32_t,uint32_t> theAsicCount_;
	int theIntegratedTime_;	
	int theEventIntegratedTime_;	
	std::map<unsigned int,DifGeom> theDifMap_;
	std::map<unsigned int,ChamberGeom> theChamberMap_;
};
#endif
