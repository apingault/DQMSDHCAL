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

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- dqm sdhcal headers
#include "Mapping.h"

namespace EVENT { class LCEvent; class LCCollection; class RawCalorimeterHit; }
namespace dqm4hep { class TiXmlElement; class TiXmlHandle; }

namespace dqm_sdhcal
{

class TriventListener
{
public:
	/** Call back function to process a Trivent
	 *  reconstructed event as a noisy event
	 */
	virtual dqm4hep::StatusCode processNoisyEvent(EVENT::LCEvent *pLCEvent) = 0;

	/** Call back function to process a Trivent
	 *  reconstructed event as a physical event
	 */
	virtual dqm4hep::StatusCode processPhysicalEvent(EVENT::LCEvent *pLCEvent) = 0;
};

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

    /** Initialize Trivent.
     *  Must be called after setting all the needed parameters and
     *  must be called once.
     */
    dqm4hep::StatusCode init();

    /** Read settings from a xml handle
     */
    dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle &xmlHandle);

	/** Process Trivent on the lcio event.
	 *  Create a CalorimeterHit collection from a RawCalorimeterHit Collection from StreamOut
	 */
    dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

    /** Add a Trivent listener
     */
    void addListener(TriventListener *pListener);

    /** Add a Trivent listener
     */
    void removeListener(TriventListener *pListener);

private:
    void clear();
    dqm4hep::StatusCode readGeometry(const std::string &fileName);
    dqm4hep::StatusCode readDifGeometry(dqm4hep::TiXmlElement *pElement);
    dqm4hep::StatusCode readChamberGeometry(dqm4hep::TiXmlElement *pElement);
    int ijkToKey(int i, int j, int k);
    int findAsicKey(int i, int j, int k);
    unsigned int getCellDifId(int cellId);
    unsigned int getCellAsicId(int cellId);
    unsigned int getCellChanId(int cellId);
    int getMaxTime();
    std::vector<int> getTimeSpectrum();
    std::vector<dqm4hep::dqm_uint> getPadIndex(unsigned int difId, unsigned int asicId, unsigned int chanId);
    dqm4hep::StatusCode eventBuilder(EVENT::LCCollection* colEvent,int timePeak, int previousTimePeak);
    bool peakOrNot(std::vector<int> timeSpectrum, int iTime, int threshold);

private:

    // trivent toolkit
    std::vector<EVENT::RawCalorimeterHit*> m_triggerRawHit;

    std::string m_inputCollectionName;
    std::string m_outputCollectionName;
    std::string m_geomXMLFile;
	float m_cellSizeU;
	float m_cellSizeV;
	float m_layerThickness;

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
    bool m_treatCherenkov;

    float m_beamEnergy;
    int m_maxTime;
    int m_evtNbr;
    int m_previousEvtNbr;
    std::set<unsigned int>    m_foundLayerList;

    int m_bcid1;
    int m_bcid2;
    bool m_cerenkovFlag[3];
    int m_cerenkovCount[3];
    int m_cerenkovCountTotal[3];

    std::map<int, LayerID> m_difMapping;
    std::map<int, double> m_chamberPositions; //chamber , position

    std::set<TriventListener *>     m_listeners;
};

}

#endif  //  Trivent_H