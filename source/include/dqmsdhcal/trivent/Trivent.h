/// \file Trivent.h
/*
*
* Trivent.h header template automatically generated by a class generator
* Creation date : lun. ao�t 3 2015
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
* @author Yacine Haddad, Arnaud Steen, Remi Ete, Antoine Pingault
* @copyright CNRS , IPNL, UGent
*/


#ifndef TRIVENT_H
#define TRIVENT_H
#define  HISTOGRAM_PARSER true

// -- dqm4hep headers
#include "dqm4hep/core/DQM4HEP.h"

namespace EVENT { class LCEvent; }
namespace dqm4hep { class TiXmlElement; }

namespace dqm_sdhcal
{

/** Trivent class
*/
class Trivent
{
public:
    /** Constructor
     */
    Trivent();

    /** Destructor
     */
    ~Trivent();

	/** Process Trivent on the lcio event.
	 *  Create a CalorimeterHit collection from a RawCalorimeterHit Collection from StreamOut
	 */
    dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

    void init();

    void    processEvent( LCEvent * evtP );

    dqm4hep::StatusCode readGeometry(const std::string &fileName);
    dqm4hep::StatusCode readDifGeometry(TiXmlElement *pElement);
    dqm4hep::StatusCode readChamberGeometry(TiXmlElement *pElement);

    int ijkToKey(int i, int j, int k);
    int findAsicKey(int i, int j, int k);
    unsigned int getCellDifId(int cellId);
    unsigned int getCellAsicId(int cellId);
    unsigned int getCellChanId(int cellId);

    int getMaxTime();
    std::vector<int> getTimeSpectrum();
    std::vector<dqm4hep::dqm_uint> getPadIndex(uint difId, uint asicId, uint chanId);
    void eventBuilder(LCCollection* colEvent,int timePeak, int previousTimePeak);
    bool peakOrNot(std::vector<int> timeSpectrum, int iTime ,int threshold);
    void end();


    void setInputCollectionName(const std::string &inputCollectionName);

    void setOutputCollectionName(const std::string &outputCollectionName);

    void setOutputFileName(const std::string &outputFileName);

    void setOutputNoiseFileName(const std::string &noiseFileName);

    //  GeometryXMLFile (Default = setup_geometry.xml)
    void setGeometryXMLFile(const std::string &geomXMLFile);

    //  Cut in number of layer (Default = 10)
    void setLayerCut(const int &layerCut);

    //  Noise Cut in time spectrum (Default = 10)
    void setNoiseCut(const int &noiseCut);

    //  TimeWindow in TimeBin (Default = 2)
    void setTimeWindow(const int &timeWindow);


    // LayerGap dimension in cm (Default = .9)
    void setLayerGap(const double &layerGap);

    // Cut on number of hit max (Electronic Noise)  (Default = 100000)
    void setElecNoiseCut(const int &elecNoiseCut);

    // Cut on time from previous event in TimeBin (Default = 0)
    void setTime2PreviousEventCut(const int &time2PreviousEventCut);

    // Activate/Deativate GainCorrectionMode (Default = False)
    void setGainCorrectionMode(const bool &gainCorrectionMode);

    // Cerenkov window for event recontstruction (Default = 20)
    void setCerenkovWindow(const int &cerenkovWindow);

    // Cerenkov signal Lenght, in TimeBin (Default = 1)
    void setCerenkovLength(const int &cerenkovLength);

    // Cerenkov DifId (Default = 3)
    void setCerenkovDifId(const int &cerenkovDifId);


    //    void setCellSize(const int &cellSize);

private:
    // xml test
    std::map<std::string,std::string> m_parameters;

    std::vector<EVENT::RawCalorimeterHit*> m_triggerRawHit;
    std::vector<EVENT::RawCalorimeterHit*> m_cerenkovHits;
    std::vector<EVENT::LCEvent*> evts; // Vector of evts reconstructed in the Trigger

    std::string _logRootName;
    std::string _mappingFile;
    std::vector<std::string> m_hcalCollections;

    std::string m_inputCollectionName;
    std::string m_outputCollectionName;
    std::string m_outputFileName;
    std::string m_noiseFileName;
    std::string m_geomXMLFile;

    int m_layerCut;
    int m_noiseCut;
    int m_timeWindow;
    double m_layerGap;
    int m_elecNoiseCut;
    int m_time2PreviousEventCut;
    bool m_gainCorrectionMode;
    int m_cerenkovWindow;
    int m_cerenkovLength;
    int m_cerenkovDifId;

    float m_beamEnergy;
    int m_trigCount;
    int m_maxTime;
    int m_evtNbr;
    int m_previousEvtNbr;
    int m_rejectedEvt;
    int m_selectedEvt;
    std::vector<dqm4hep::dqm_uint> m_index;
    uintVec _zCut;

    int m_bcid1;
    int m_bcid2;
    bool m_cerenkovFlag[3];
    int m_cerenkovCount[3];
    int m_cerenkovCountTotal[3];

    std::map<int, LayerID> m_difMapping;
    std::map<int, double> m_chamberPositions; //chamber , position
};

}

#endif  //  Trivent_H
