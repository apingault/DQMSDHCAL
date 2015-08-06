/// \file Streamout.cc
/*
*
* Streamout.cc source template automatically generated by a class generator
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

// -- dqm sdhcal headers
#include <dqmsdhcal/trivent/Trivent.h>

// -- lcio headers
#include <EVENT/LCIO.h>
#include <EVENT/LCCollection.h>
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// -- std headers
#include <limits.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <algorithm>

// -- dqm4hep headers
#include "dqm4hep/core/tinyxml.h"

namespace dqm_sdhcal
{

Trivent::Trivent():
    m_inputCollectionName("DHCALRawHits"),
    m_outputCollectionName("SDHCAL_HIT"),
    m_layerCut(10),
    m_noiseCut(10),
    m_timeWindow(2),
    m_geomXMLFile("setup_geometry.xml"),
    m_layerGap(.9),
    m_elecNoiseCut(100000),
    m_time2PreviousEventCut(0),
    m_gainCorrectionMode(false),
    m_cerenkovWindow(20),
    m_cerenkovLength(1),
    m_cerenkovDifId(3),
	m_cellSizeU(10.408f),
	m_cellSizeV(10.408f),
	m_layerThickness(26.131)
{
	/* nop */
}

//-------------------------------------------------------------------------------------------------

Trivent::~Trivent()
{
	clear();
}

//-------------------------------------------------------------------------------------------------

void Trivent::clear()
{
	for(std::vector<EVENT::LCEvent*>::iterator iter = m_reconstructedEvents.begin(), endIter = m_reconstructedEvents.end() ;
			endIter != iter ; ++iter)
	{
		delete *iter;
	}

	for(std::vector<EVENT::LCEvent*>::iterator iter = m_noiseEvents.begin(), endIter = m_noiseEvents.end() ;
			endIter != iter ; ++iter)
	{
		delete *iter;
	}

	m_reconstructedEvents.clear();
	m_noiseEvents.clear();
}

//-------------------------------------------------------------------------------------------------

const std::vector<EVENT::LCEvent*> &Trivent::getReconstructedEvents() const
{
	return m_reconstructedEvents;
}

//-------------------------------------------------------------------------------------------------

const std::vector<EVENT::LCEvent*> &Trivent::getNoiseEvents() const
{
	return m_noiseEvents;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::readGeometry(const std::string &fileName)
{
	dqm4hep::TiXmlDocument document(fileName);

	if(!document.LoadFile())
		return dqm4hep::STATUS_CODE_FAILURE;

	dqm4hep::TiXmlHandle documentHandle(&document);
	dqm4hep::TiXmlElement* pRootElement = documentHandle.FirstChildElement().Element();

	if(NULL == pRootElement)
		return dqm4hep::STATUS_CODE_FAILURE;

	// root element handler
	dqm4hep::TiXmlHandle rootHandle(pRootElement);

	bool difGeomFound = false;
	bool chamberGeomFound = false;

	for(dqm4hep::TiXmlElement *pParameterElement = rootHandle.FirstChild("parameter").Element(); NULL != pParameterElement;
			pParameterElement = pParameterElement->NextSiblingElement("parameter"))
        {
		// parameter name
		const char *pParameterNameStr = pParameterElement->Attribute("name");

		if(NULL == pParameterNameStr)
			return dqm4hep::STATUS_CODE_NOT_FOUND;

		std::string parameterName = pParameterNameStr;

		if(parameterName == "DifGeom")
		{
			RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->readDifGeometry(pParameterElement));
			difGeomFound = true;
		}
		else if(parameterName == "ChamberGeom")
		{
			RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->readChamberGeometry(pParameterElement));
			chamberGeomFound = true;
        }
		else
        {
			std::cout << "Unknown parameter element, name : " << parameterName << std::endl;
			continue;
		}
	}

	if(!difGeomFound || !chamberGeomFound)
		return dqm4hep::STATUS_CODE_FAILURE;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::readDifGeometry(dqm4hep::TiXmlElement *pElement)
{
	if(NULL == pElement)
		return dqm4hep::STATUS_CODE_FAILURE;

	// get the element contents
	const char *pValueStr = pElement->GetText();

	if(NULL == pValueStr)
		return dqm4hep::STATUS_CODE_NOT_FOUND;

	// clear dif mapping
	m_difMapping.clear();

	std::string value = pValueStr;
    std::vector<std::string> lines;

	// split the contents for each lines
	dqm4hep::DQM4HEP::tokenize(value, lines, "\n");

	for(unsigned int i=0 ; i<lines.size() ; i++)
	{
		std::string line = lines.at(i);
		std::vector<std::string> tokens;

		LayerID layerId;
                    int difId;

		// split the line with commas
		dqm4hep::DQM4HEP::tokenize(line, tokens, ",");

		// fill LayerID object
		dqm4hep::DQM4HEP::stringToType(tokens.at(0), difId);
		dqm4hep::DQM4HEP::stringToType(tokens.at(1), layerId.K);
		dqm4hep::DQM4HEP::stringToType(tokens.at(2), layerId.DifX);
		dqm4hep::DQM4HEP::stringToType(tokens.at(3), layerId.DifY);
		dqm4hep::DQM4HEP::stringToType(tokens.at(4), layerId.IncX);
		dqm4hep::DQM4HEP::stringToType(tokens.at(5), layerId.IncY);

		// add the dif entry
		m_difMapping[difId] = layerId;
    }

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::readChamberGeometry(dqm4hep::TiXmlElement *pElement)
{
	if(NULL == pElement)
		return dqm4hep::STATUS_CODE_FAILURE;

	// get the element contents
	const char *pValueStr = pElement->GetText();

	if(NULL == pValueStr)
		return dqm4hep::STATUS_CODE_NOT_FOUND;

	// clear the chamber positions
	m_chamberPositions.clear();

	std::string value = pValueStr;
	std::vector<std::string> lines;

	// split the contents by lines
	dqm4hep::DQM4HEP::tokenize(value, lines, "\n");

	for(unsigned int i=0; i<lines.size(); i++)
	{
		std::string line = lines.at(i);
		std::vector<std::string> tokens;

		double position;
		int difId;

		// split the lines with commas
		dqm4hep::DQM4HEP::tokenize(line, tokens, ",");

		dqm4hep::DQM4HEP::stringToType(tokens.at(0), difId);
		dqm4hep::DQM4HEP::stringToType(tokens.at(0), position);

		// add a chamber entry
		m_chamberPositions[difId] = position;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

// ============ decode the cell ids =============
// Dif 1 => cellID0= 00983297 => DifID=1 / AsicID=1 / ChanID=15

// Applique le masque 1111 1111
unsigned int Trivent::getCellDifId(int cellId)
{
    return cellId & 0xFF;
}
//-------------------------------------------------------------------------------------------------

//  Applique le masque 1111 1111 0000 0000 puis tronque les 8 derniers bits
unsigned int Trivent::getCellAsicId(int cellId)
{
    return (cellId & 0xFF00)>>8;
}

//-------------------------------------------------------------------------------------------------

//  Applique le masque 1111 0000 0000 0000 0000 puis tronque les 16 derniers bits
unsigned int Trivent::getCellChanId(int cellId)
{
    return (cellId & 0x3F0000)>>16;
}

// ============ ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============ ============
// Exemple sur Dif 1:
// cellId0 = 00983297 -> Binaire =  1111 0000 0001 0000 0001
// binaire & 0xFF = 0000 0001 => 2^0=1
// binaire & 0xFF00 = 0000 0001 0000 0000 >>8 = 0000 0001 => 2^0=1
// binaire & 0x3F0000 =  1111 0000 0000 0000 0000 >>16 = 1111 => (2^3)+(2^2)+(2^1)+(2^0)=15
// ============ ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============ ============

//-------------------------------------------------------------------------------------------------

std::vector<dqm4hep::dqm_uint> Trivent::getPadIndex(unsigned int difId, unsigned int asicId, unsigned int chanId)
{
	std::vector<dqm4hep::dqm_uint> index(3, 0);

    int difY = static_cast<int>(m_difMapping.find(difId)->second.DifY);
    int difZ = static_cast<int>(m_difMapping.find(difId)->second.K);
    index[0] = ( 1 + MapILargeHR2[chanId] + AsicShiftI[asicId] );

    index[1] = ( 32 - (MapJLargeHR2[chanId]+AsicShiftJ[asicId]) ) + difY;
    index[2] = abs(difZ);

    return index;
}

//-------------------------------------------------------------------------------------------------

int Trivent::getMaxTime()
{
    int maxTime = 0;

	for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawHit = m_triggerRawHit.begin() ;
			rawHit != m_triggerRawHit.end() ; rawHit++)
	{
		int time = static_cast<int>((*rawHit)->getTimeStamp());

		if(time >= 0)
			m_maxTime = std::max(m_maxTime, time);
     }

	return maxTime;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> Trivent::getTimeSpectrum()
{
	int maxTime = this->getMaxTime();
    std::vector<int> timeSpectrum(maxTime + 1);

	for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawHit = m_triggerRawHit.begin() ;
			rawHit != m_triggerRawHit.end() ; rawHit++)
	{
		int time = int((*rawHit)->getTimeStamp());

		if(time >= 0)
			timeSpectrum.at(time)++;
    }

    return timeSpectrum;
}

//-------------------------------------------------------------------------------------------------

bool Trivent::peakOrNot(std::vector<int> timeSpectrum , int iTime, int threshold)
{
	return (timeSpectrum.at(iTime) >= threshold
		  && timeSpectrum.at(iTime) > timeSpectrum.at(iTime+1));
}

//-------------------------------------------------------------------------------------------------

int Trivent::ijkToKey(int i, int j, int k)
{
	return 100*100*k+100*j+i;
}

//-------------------------------------------------------------------------------------------------

int Trivent::findAsicKey(int i, int j, int k)
{
    if(i>96 || i<0 || j>96 || j<0)
    	return -1;

    return k*1000+(((j-1)/8)*12 + (i-1)/8);
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::eventBuilder(EVENT::LCCollection *pLCCollection, int timePeak, int previousTimePeak)
{
    _zCut.clear();

    pLCCollection->setFlag(pLCCollection->getFlag()|( 1 << EVENT::LCIO::RCHBIT_LONG));
    pLCCollection->setFlag(pLCCollection->getFlag()|( 1 << EVENT::LCIO::RCHBIT_TIME));

    UTIL::CellIDEncoder<IMPL::CalorimeterHitImpl> cellIDDecoder( "M:3,S-1:3,I:9,J:9,K-1:6", pLCCollection) ;

    std::map<int, int> asicMap;
	std::vector<int> hitKeys;
	int previousTime = 0;
	int cerenkovTime = 0;

	for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawHit = m_triggerRawHit.begin() ;
			rawHit != m_triggerRawHit.end() ; rawHit++)
	{
		int time = (*rawHit)->getTimeStamp();
		int difId = getCellDifId((*rawHit)->getCellID0());

		// CTag has an offset of ~5 in December 2014, 7-17 in May/April 2015
		// Need to modify the timeWin in the following loop to analyse it
		if (difId == m_cerenkovDifId)
			m_timeWindow = m_cerenkovWindow;
		else
			m_timeWindow = 2;

		//
		if(fabs(time-timePeak) <= m_timeWindow &&
				(time > previousTimePeak + m_timeWindow ))
		{
			int asicId = getCellAsicId((*rawHit)->getCellID0());
			int channelId = getCellChanId((*rawHit)->getCellID0());
			std::vector<dqm4hep::dqm_uint> padIndex = getPadIndex(difId, asicId, channelId);

			dqm4hep::dqm_uint I = padIndex[0];
			dqm4hep::dqm_uint J = padIndex[1];
			dqm4hep::dqm_uint K = padIndex[2];

			// Find and remove square events
			int asicKey = findAsicKey(I,J,K);

			if(asicMap[asicKey])
				asicMap[asicKey]++;
			else
				asicMap[asicKey] = 1;

			if(asicMap[asicKey] == 64)
			{
				_zCut.clear();
				hitKeys.clear();
				asicMap.clear();

				return dqm4hep::STATUS_CODE_SUCCESS;
			}

			int aHitKey = ijkToKey(I,J,K);

			float pos[3];
			pos[0] = I*m_cellSizeU;
			pos[1] = J*m_cellSizeV;
			pos[2] = K*m_layerThickness;

			// If No more CTag Reset it to 0
			if ( cerenkovTime - previousTime > 1)
			{
				m_cerenkovFlag[0] = m_cerenkovFlag[1] = m_cerenkovFlag[2] = 0;
				m_cerenkovCount[0] = m_cerenkovCount[1] = m_cerenkovCount[2] = 0;
			}

			// Treat CTag
			if (difId == m_cerenkovDifId)
			{
				cerenkovTime = time;
				unsigned short cerenkovAmplitude = (*rawHit)->getAmplitude();

				// In december 2014 only one Cerenkov was used and its signal spanned on 4-5times slot
				// Need to be treated separately to not count the CTag multiple time
				// Besides in Decmeber cerenkovAmplitude for the CTag is equal to 5
				if (cerenkovAmplitude == 5)
				{
					m_cerenkovCount[2]++;
					m_cerenkovCountTotal[2]++;

					if (m_cerenkovFlag[2] == 0)
						m_cerenkovFlag[2] = 1;
				}
				else if (cerenkovAmplitude<3) // In may/April 2015 2 Cerenkov. cerenkovAmplitude = 1 = CTag1
									//                                                             = 2 = CTag2
									//                                                             = 3 = CTag1 + CTag2
				{
					m_cerenkovCount[cerenkovAmplitude-1]++;
					m_cerenkovCountTotal[cerenkovAmplitude-1]++;

					if (m_cerenkovFlag[cerenkovAmplitude-1] == 0)
						m_cerenkovFlag[cerenkovAmplitude-1] = 1;
				}
			}

			// last chamber id
			std::map<int, double>::iterator it = m_chamberPositions.end();
			it--;
			int lastChamberId = it->first;

			if(K < 0 || K > lastChamberId)
			{
				streamlog_out( ERROR ) << " difId  == " << difId
										 << " asicId == " << asicId
										 << " channelId == " << channelId
										 << " I == " << I
										 << " J == " << J
										 << " K == " << K
										 << std::endl;
				continue;
			}

			IMPL::CalorimeterHitImpl *pCaloHit = new IMPL::CalorimeterHitImpl();
			pCaloHit->setTime(float((*rawHit)->getTimeStamp()));

			if(float((*rawHit)->getAmplitude()&3)>2.5)
				pCaloHit->setEnergy(float((*rawHit)->getAmplitude()&3));         // 3rd treshold
			else if(float((*rawHit)->getAmplitude()&3)>1.5)
				pCaloHit->setEnergy(float((*rawHit)->getAmplitude()&3)-1);  // 2nd treshold -1 to shift color to green
			else
				pCaloHit->setEnergy(float((*rawHit)->getAmplitude()&3)+1);                                             // 1st treshold +1 to shift color to blue

			//avoid two hits in the same cell
			std::vector<int>::iterator findIter = std::find(hitKeys.begin(), hitKeys.end(), aHitKey);

			if(findIter != hitKeys.end())
			{
				IMPL::CalorimeterHitImpl* hit =
						dynamic_cast<IMPL::CalorimeterHitImpl*>(pLCCollection->getElementAt(std::distance(hitKeys.begin(), findIter)));

				float hitTime = hit->getTime();

				if( fabs(timePeak - hitTime) > fabs( timePeak - time ))
					hit->setEnergy(pCaloHit->getEnergy());

				continue;
			}

			// set the cell id
			cellIDDecoder["I"] = I;
			cellIDDecoder["J"] = J;
			cellIDDecoder["K-1"] = K-1;
			cellIDDecoder["M"] = 0;
			cellIDDecoder["S-1"] = 3;

			cellIDDecoder.setCellID(pCaloHit);

			if(std::find(_zCut.begin(), _zCut.end(), K)==_zCut.end())
				_zCut.push_back(K);

			pCaloHit->setPosition(pos);
			pLCCollection->addElement(pCaloHit);
			hitKeys.push_back(aHitKey);

			m_previousEvtNbr = m_evtNbr;
		}

		previousTime = time;
	}

	hitKeys.clear();

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::init()
{
    m_cerenkovFlag[0] = m_cerenkovFlag[1] = m_cerenkovFlag[2] = 0;
    m_cerenkovCount[0] = m_cerenkovCount[1] = m_cerenkovCount[2] = 0;
    m_cerenkovCountTotal[0] = m_cerenkovCountTotal[1] = m_cerenkovCountTotal[2] = 0;

    RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, readGeometry(m_geomXMLFile));

    m_evtNbr = 0;// event number
    m_previousEvtNbr = 0;

    return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode Trivent::processEvent(EVENT::LCEvent *pLCEvent)
{
    if (NULL == pLCEvent)
        return dqm4hep::STATUS_CODE_INVALID_PTR;

    m_evtNbr = pLCEvent->getEventNumber();

	// Grab the input collection
	EVENT::LCCollection * pLCCollection = NULL;

	try
	{
		pLCCollection = pLCEvent->getCollection(m_inputCollectionName);
	}
	catch (lcio::DataNotAvailableException &exception)
	{
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	if(NULL == pLCCollection)
		return dqm4hep::STATUS_CODE_FAILURE;

	if(pLCCollection->getTypeName() != EVENT::LCIO::RAWCALORIMETERHIT)
		return dqm4hep::STATUS_CODE_INVALID_PARAMETER;

	int numberOfHits = pLCCollection->getNumberOfElements(); // hit number

	// If numberOfHit too large do not process the event
	if(numberOfHits > m_elecNoiseCut)
	{
		streamlog_out( MESSAGE ) << "TRIGGER SKIPED ... NoiseCut" << std::endl;
		return dqm4hep::STATUS_CODE_SUCCESS;
	}

	// set raw hits
	m_triggerRawHit.clear();
	std::vector<int> vTrigger;

	// clear previous events
	clear();

	for (int iHit=0 ; iHit<numberOfHits; iHit++) // loop over the hits
	{
		EVENT::RawCalorimeterHit *pRawHit = dynamic_cast<EVENT::RawCalorimeterHit*>( pLCCollection->getElementAt(iHit));

		// just in case ...
		if(NULL == pRawHit)
			continue;

		// Get the difId
		unsigned int difId = pRawHit->getCellID0()&0xFF;

		// Extract abolute bcid information
		if(iHit == 0)
		{
			if (difId == 0)
				continue;

			std::stringstream pname("");
			pname << "DIF" << difId << "_Triggers";

			pLCCollection->getParameters().getIntVals(pname.str(),vTrigger);

			if (vTrigger.size() != 0)
			{
				m_bcid1 = vTrigger[4]; // absoluteBCID/(0xFFFFFF+1))&0xFFFFFF;
				m_bcid2 = vTrigger[3]; // absoluteBCID&0xFFFFFF;

				// Shift the value from the 24 first bits
				unsigned long long Shift = 16777216ULL;
				unsigned long long theBCID_ = m_bcid1*Shift + m_bcid2;

				streamlog_out( DEBUG1 ) << "trigger time : " << theBCID_ << std::endl;
			}
		}

		m_triggerRawHit.push_back(pRawHit);
	}

	std::vector<int> timeSpectrum = getTimeSpectrum();

	//---------------------------------------------------------------
	//! Find the candidate event
	//!
	int iBin=0; // The bin number
	int previousBin = 0;

	while(iBin < (m_maxTime+1))
	{
		if(timeSpectrum[iBin] >= m_noiseCut                &&
				timeSpectrum[iBin] >= timeSpectrum[iBin+2] &&
				timeSpectrum[iBin] >= timeSpectrum[iBin+1] &&
				timeSpectrum[iBin] >= timeSpectrum[iBin-1] &&
				timeSpectrum[iBin] >= timeSpectrum[iBin-2]
				)
		{
			IMPL::LCEventImpl *pEvt = new IMPL::LCEventImpl() ;     // create the event

			//---------- set event paramters ------
			const std::string parname_trigger = "trigger";
			const std::string parname_energy  = "beamEnergy";
			const std::string parname_bcid1 = "bcid1";
			const std::string parname_bcid2 = "bcid2";
			const std::string parname_cerenkov1 = "cerenkov1"; // First Cerenkov in April/May 2015
			const std::string parname_cerenkov2 = "cerenkov2"; // Second Cerenkov in April/May 2015
			const std::string parname_cerenkov3 = "cerenkov3"; // Both Cerenkov in April/May 2015 + Cerenkov in December 2014

			pEvt->parameters().setValue(parname_trigger, pLCEvent->getEventNumber());
			pEvt->parameters().setValue(parname_energy, m_beamEnergy);
			pEvt->parameters().setValue(parname_bcid1, m_bcid1);
			pEvt->parameters().setValue(parname_bcid2, m_bcid2);

			pEvt->setRunNumber( pLCEvent->getRunNumber()) ;

			//-------------------------------------
			IMPL::LCCollectionVec* pCalorimeterHitCollection = new IMPL::LCCollectionVec(EVENT::LCIO::CALORIMETERHIT);

			m_evtNbr++;
			m_cerenkovCount[0] = m_cerenkovCount[1] = m_cerenkovCount[2] = 0;

			Trivent::eventBuilder(pCalorimeterHitCollection, iBin, previousBin);
			pEvt->parameters().setValue(parname_cerenkov1, m_cerenkovCount[0]); // Value determined in the eventBuilder
			pEvt->parameters().setValue(parname_cerenkov2, m_cerenkovCount[1]); // Value determined in the eventBuilder
			pEvt->parameters().setValue(parname_cerenkov3, m_cerenkovCount[2]); // Value determined in the eventBuilder
			// ->Need to be after the EventBuilder function

			pEvt->setEventNumber(m_evtNbr) ;

			// Add the collection to the event
			pEvt->addCollection(pCalorimeterHitCollection, m_outputCollectionName);

			if( (int)_zCut.size() > m_layerCut &&                                  // the min layer numb cut
					abs( iBin - previousBin) > m_time2PreviousEventCut)// time2prev event  cut
			{
				m_reconstructedEvents.push_back(pEvt);
			}
			else
			{
				m_noiseEvents.push_back(pEvt);
			}

			previousBin = iBin;
			iBin = iBin + m_timeWindow;
		}
		else
		{
			iBin++;
		}
	}

    return dqm4hep::STATUS_CODE_SUCCESS;
}

//==============================================================

void Trivent::setInputCollectionName(const std::string &inputCollectionName)
{
    m_inputCollectionName = inputCollectionName;
}

void Trivent::setOutputCollectionName(const std::string &outputCollectionName)
{
    m_outputCollectionName = outputCollectionName;
}

void Trivent::setGeometryXMLFile(const std::string &geomXMLFile)
{
    m_geomXMLFile = geomXMLFile;
}

void Trivent::setLayerCut(const int &layerCut)
{
    m_layerCut = layerCut;
}

void Trivent::setNoiseCut(const int &noiseCut)
{
    m_noiseCut = noiseCut;
}

void Trivent::setTimeWindow(const int &timeWindow)
{
    m_timeWindow = timeWindow;
}

void Trivent::setLayerGap(const double &layerGap)
{
    m_layerGap = layerGap;
}

void Trivent::setElecNoiseCut(const int &elecNoiseCut)
{
    m_elecNoiseCut = elecNoiseCut;
}

void Trivent::setTime2PreviousEventCut(const int &time2PreviousEventCut)
{
    m_time2PreviousEventCut = time2PreviousEventCut;
}

void Trivent::setGainCorrectionMode(const bool &gainCorrectionMode)
{
    m_gainCorrectionMode = gainCorrectionMode;
}

void Trivent::setCerenkovWindow(const int &cerenkovWindow)
{
    m_cerenkovWindow = cerenkovWindow;
}

void Trivent::setCerenkovLength(const int &cerenkovLength)
{
    m_cerenkovLength = cerenkovLength;
}

void Trivent::setCerenkovDifId(const int &cerenkovDifId)
{
    m_cerenkovDifId = cerenkovDifId;
}

void Trivent::setCellSizes(const float &cellSizeU, const float &cellSizeV)
{
	m_cellSizeU = cellSizeU;
	m_cellSizeV = cellSizeV;
}

}
