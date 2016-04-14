  /// \file ShmProcessors.h
/*
 *
 * ShmProcessors.h header template automatically generated by a class generator
 * Creation date : jeu. avr. 14 2016
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
 * @author Remi Ete
 * @copyright CNRS , IPNL
 */


#ifndef SHMPROCESSORS_H
#define SHMPROCESSORS_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/evb/DQMShmProcessor.h"

// -- levbdim headers
#include "buffer.hh"

namespace IMPL { class LCFlagImpl; class RawCalorimeterHitImpl; }

namespace dqm_sdhcal
{

class DIFPtr;

/** EventInfoShmProcessor class
 *  Fill the EVENT::LCEvent parameters
 */
class EventInfoShmProcessor : public dqm4hep::DQMShmProcessor
{
public:
	/** Constructor
	 */
	EventInfoShmProcessor();

	/** Destructor
	 */
	~EventInfoShmProcessor();

	/** Call back function on start of run
	 */
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);

	/** Call back function on end of run
	 */
	dqm4hep::StatusCode endOfRun(const dqm4hep::DQMRun *const pRun);

	/** Called when an event is reconstructed.
	 *  The key is a unique identifier for the event.
	 *  The buffer list is the reconstructed list of buffers for the target key
	 *  of the different sources
	 */
	dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *pEvent, uint32_t key, const std::vector<levbdim::buffer*> &bufferList);

	/** Read settings from xml handle
	 */
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);

private:
	unsigned int                      m_eventNumber;
	unsigned int                      m_runNumber;
	std::string                       m_detectorName;
	std::string                       m_creationTimeParameterName;
};

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

/** SDHCALDifHelper class
 */
class SDHCALDifHelper
{
public:
	/** Create a sdhcal dif pointer from the raw buffer.
	 */
	static DIFPtr *createDIFPtr(levbdim::buffer *pBuffer, unsigned int xdaqShift);

	/** Create a raw calorimeter hit from the dif ptr at a target frame and channel
	 */
	static IMPL::RawCalorimeterHitImpl *createRawCalorimeterHit(DIFPtr *pDifPtr, uint32_t frame, uint32_t channel);

	/** Whether the pad has been fired in the dif
	 */
	static bool isEmptyPad(DIFPtr *pDifPtr, uint32_t frame, uint32_t channel);

	/** Append the dif info at the end of vector
	 */
	static void fillDifTriggerInfo(DIFPtr *pDifPtr, dqm4hep::IntVector &vector);

	/** Fill the lc flag for a sdhcal raw calorimeter hit collection
	 */
	static void setLCFlag(IMPL::LCFlagImpl &lcFlag);
};

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

/** SDHCALShmProcessor class
 *  Produce a collection of EVENT::RawCalorimeterHit from the dif raw buffers
 */
class SDHCALShmProcessor : public dqm4hep::DQMShmProcessor
{
	typedef std::set<unsigned int> UIntSet;
public:
	/** Constructor
	 */
	SDHCALShmProcessor();

	/** Destructor
	 */
	~SDHCALShmProcessor();

	/** Call back function on start of run
	 */
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);

	/** Call back function on end of run
	 */
	dqm4hep::StatusCode endOfRun(const dqm4hep::DQMRun *const pRun);

	/** Called when an event is reconstructed.
	 *  The key is a unique identifier for the event.
	 *  The buffer list is the reconstructed list of buffers for the target key
	 *  of the different sources
	 */
	dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *pEvent, uint32_t key, const std::vector<levbdim::buffer*> &bufferList);

	/** Read settings from xml handle
	 */
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);

private:
	bool                                   m_skipFullAsics;
	bool                                   m_dropFirstRU;
	unsigned int                           m_xdaqShift;
	unsigned int                           m_detectorId;
	unsigned int                           m_noiseLimit;
	std::string                            m_outputCollectionName;
	UIntSet                                m_difMaskList;
};

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

/** CherenkovShmProcessor class
 */
class CherenkovShmProcessor : public dqm4hep::DQMShmProcessor
{
public:
	/** Constructor
	 */
	CherenkovShmProcessor();

	/** Destructor
	 */
	~CherenkovShmProcessor();

	/** Call back function on start of run
	 */
	dqm4hep::StatusCode startOfRun(dqm4hep::DQMRun *const pRun);

	/** Call back function on end of run
	 */
	dqm4hep::StatusCode endOfRun(const dqm4hep::DQMRun *const pRun);

	/** Called when an event is reconstructed.
	 *  The key is a unique identifier for the event.
	 *  The buffer list is the reconstructed list of buffers for the target key
	 *  of the different sources
	 */
	dqm4hep::StatusCode processEvent(dqm4hep::DQMEvent *pEvent, uint32_t key, const std::vector<levbdim::buffer*> &bufferList);

	/** Read settings from xml handle
	 */
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle xmlHandle);

private:
	unsigned int                           m_detectorId;
	unsigned int                           m_xdaqShift;
	unsigned int                           m_cherenkovDifId;
	int                                    m_cherenkovTimeShift;
	std::string                            m_outputCollectionName;
};

} 

#endif  //  SHMPROCESSORS_H
