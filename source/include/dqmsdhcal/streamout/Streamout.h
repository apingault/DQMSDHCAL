  /// \file Streamout.h
/*
 *
 * Streamout.h header template automatically generated by a class generator
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
 * @author Laurent Mirabito, Remi Ete
 * @copyright CNRS , IPNL
 */


#ifndef STREAMOUT_H
#define STREAMOUT_H

#include "dqm4hep/DQM4HEP.h"

namespace EVENT { class LCEvent; }

namespace dqm_sdhcal
{

/** Streamout class
 */
class Streamout 
{
public:
	/** Constructor
	 */
	Streamout();

	/** Destructor
	 */
	~Streamout();

	/** Process streamout on the lcio event.
	 *  Create a RawCalorimeterHit collection from a LCGenericObject
	 *  collection of sdhcal raw dif buffers
	 */
	dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

	/** Set the RU shift
	 */
	void setRuShift(int ruShift);

	/** Set the collection name used as input for streamout.
	 *  The input collection must be a LCGenericObject collection
	 */
	void setInputCollectionName(const std::string &collectionName);

	/** Set the collection name used as output after streamout processing.
	 *  The out put collection is a RawCalorimeterHit collection
	 */
	void setOutputCollectionName(const std::string &collectionName);

	/** Set whether the first RU in the collection has to be dropped
	 */
	void setDropFirstRU(bool drop);

	/** Set the XDAQ shift (dif ptr)
	 */
	void setXDaqShift(unsigned int shift);

	/** Set whether to skip full asics (default is true)
	 */
	void setSkipFullAsic(bool skip);

private:

	int                m_ruShift;
	unsigned int      m_xdaqShift;
	bool               m_dropFirstRU;
	bool               m_skipFullAsics;
	std::string        m_inputCollectionName;
	std::string        m_outputCollectionName;
}; 

} 

#endif  //  STREAMOUT_H
