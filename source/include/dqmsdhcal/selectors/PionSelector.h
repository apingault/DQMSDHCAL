  /// \file PionSelector.h
/*
 *
 * PionSelector.h header template automatically generated by a class generator
 * Creation date : mar. sept. 1 2015
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


#ifndef PIONSELECTOR_H
#define PIONSELECTOR_H

// -- dqm4hep headers
#include "dqm4hep/DQM4HEP.h"

// -- lcio headers
#include "EVENT/LCEvent.h"
#include "EVENT/CalorimeterHit.h"

namespace dqm4hep { class TiXmlHandle; }

namespace dqm_sdhcal
{

/** PionSelector class
 */ 
class PionSelector 
{
public:
	/** Constructor
	 */
	PionSelector();

	/** Destructor
	 */
	~PionSelector();

	/** Read the settings needed to configure the pion event selector
	 *  from a xml handle
	 */
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle &value);

	/** Process an event and fills the properties to decide
	 *  whether the event is pion event
	 */
	dqm4hep::StatusCode processEvent(EVENT::LCEvent *pLCEvent);

	/** Whether the last processed event was a pion event
	 */
	bool isPion() const;

private:
	/** Convert i j k to a unique key
	 */
	int IJKToKey( const int i , const int j , const int k );

	/** Convert the key to i j and k
	 */
	std::vector<int> KeyToIJK( const int &key );

	/** Static function for std::sort() to sort the calo hit collection by increasing layer
	 */
	static bool sortByLayer( EVENT::CalorimeterHit *caloHit1 , EVENT::CalorimeterHit *caloHit2 );

	/** Compute the fractal dimension
	 */
	double getFractalDimension();

	/** Helper for function for fractal dimension computation
	 */
	int nbOfHitsInCube( int cubeSize );

    /** Determine whether the event is a single particle shower event
     */
    bool isSingleParticle();

private:
	std::string m_decoderString;
	std::string m_sdhcalCollectionName;
	std::vector<std::string> m_ijkEncoding;

	std::vector<EVENT::CalorimeterHit*> m_caloHitCollection;

	bool m_isPion;

	int m_nbOfLayers;
	int m_nbOfCellSizeX;
	int m_nbOfCellSizeY;
	std::vector<int> m_nlayers;
	std::vector<int> m_nHit;
	std::vector<double> m_cog;
	double m_radius;
	double m_fractalDimension;
	int m_showerStartingLayer;
	int m_nHitsInEdge;
	int m_nHolesAfterStartingPoint;
	int m_nHitsInCentralCells;


	// processor cut parameters. See ctor for info
	int               m_nHitCut;
	double           m_nHitOverNLayerCut;
	double           m_radiusOverCog2Cut;
	int              m_showerStartingLayerCut;
	int              m_nTouchedLayersCut;
	double          m_fractalTimesCentralCellsCut;
	int              m_nHolesCut;
	double          m_nHitEdgePercentCut;
	int              m_electronStartingLayerCut;
	double           m_largeRMSCut;
	double           m_barycenterPositionCut;
	double           m_cosThetaCut;
	int               m_neutralFirstLayerCut;
}; 

} 

#endif  //  PIONSELECTOR_H
