  /// \file PionSelector.cc
/*
 *
 * PionSelector.cc source template automatically generated by a class generator
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

// -- dqmsdhcal headers
#include "dqmsdhcal/selectors/PionSelector.h"

// -- lcio headers
#include "UTIL/CellIDDecoder.h"

// -- std headers
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Linear3DFit.hh"

namespace dqm_sdhcal
{

PionSelector::PionSelector() :
		m_isPion(false),
		m_decoderString("M:3,S-1:3,I:9,J:9,K-1:6"),
		m_sdhcalCollectionName("SDHCAL_HIT"),
		m_nbOfLayers(48),
		m_nbOfCellSizeX(96),
		m_nbOfCellSizeY(96),
		m_radius(0.f),
		m_fractalDimension(0.f),
		m_showerStartingLayer(-1),
		m_nHitsInEdge(0),
		m_nHolesAfterStartingPoint(0),
		m_nHitsInCentralCells(0),
		m_nHitCut(100),
		m_nHitOverNLayerCut(3),
		m_radiusOverCog2Cut(0.4),
		m_showerStartingLayerCut(0),
		m_nTouchedLayersCut(0),
		m_fractalTimesCentralCellsCut(0),
		m_nHolesCut(0),
		m_nHitEdgePercentCut(0.5),
		m_electronStartingLayerCut(5),
		m_largeRMSCut(4.0),
		m_barycenterPositionCut(4.0),
		m_cosThetaCut(0.95),
		m_neutralFirstLayerCut(5)
{
	m_ijkEncoding.push_back("I");
	m_ijkEncoding.push_back("J");
	m_ijkEncoding.push_back("K-1");
}

//-------------------------------------------------------------------------------------------------

PionSelector::~PionSelector() 
{

}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode PionSelector::readSettings(const Json::Value &value)
{
	try
	{
		m_nHitCut = value.get("nHitCut", m_nHitCut).asInt();
		m_nHitOverNLayerCut = value.get("NHitOverNLayerCut", m_nHitOverNLayerCut).asDouble();
		m_radiusOverCog2Cut = value.get("RadiusOverCog2Cut", m_radiusOverCog2Cut).asDouble();
		m_showerStartingLayerCut = value.get("ShowerStartingLayerCut", m_showerStartingLayerCut).asInt();
		m_nTouchedLayersCut = value.get("NTouchedLayersCut", m_nTouchedLayersCut).asInt();
		m_fractalTimesCentralCellsCut = value.get("FractalTimesCentralCellsCut", m_fractalTimesCentralCellsCut).asDouble();
		m_nHolesCut = value.get("NHolesCut", m_nHolesCut).asInt();
		m_nHitEdgePercentCut = value.get("NHitEdgePercentCut", m_nHitEdgePercentCut).asDouble();
		m_electronStartingLayerCut = value.get("ElectronStartingLayerCut", m_electronStartingLayerCut).asInt();
		m_largeRMSCut = value.get("m_largeRMSCut", m_largeRMSCut).asDouble();
		m_barycenterPositionCut = value.get("m_barycenterPositionCut", m_barycenterPositionCut).asDouble();
		m_cosThetaCut = value.get("m_cosThetaCut", m_cosThetaCut).asDouble();
		m_neutralFirstLayerCut = value.get("m_neutralFirstLayerCut", m_neutralFirstLayerCut).asInt();
	}
	catch(const std::runtime_error &exception)
	{
		streamlog_out(ERROR) << "PionSelector::readSettings(v): " << exception.what() << std::endl;
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

int PionSelector::IJKToKey( const int i , const int j , const int k )
{
	return 100*100*k+100*j+i;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> PionSelector::KeyToIJK( const int &key )
{
	std::vector<int> vec;
	vec.push_back( key%100 );
	vec.push_back( key/100%100 );
	vec.push_back( key/10000 );
	return vec;
}

//-------------------------------------------------------------------------------------------------

bool PionSelector::sortByLayer( EVENT::CalorimeterHit *caloHit1 , EVENT::CalorimeterHit *caloHit2 )
{
	return ( caloHit1->getPosition()[2] < caloHit2->getPosition()[2] );
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode PionSelector::processEvent(EVENT::LCEvent *pLCEvent)
{
	m_nlayers = std::vector<int>( m_nbOfLayers , 0 );
	m_radius = 0.0;
	m_showerStartingLayer = -1;
	m_nHitsInEdge = 0;
	m_nHolesAfterStartingPoint = 0;
	m_fractalDimension = 0.0;
	m_nHitsInCentralCells = 0;

	m_nHit.clear();
	m_cog.clear();
	m_caloHitCollection.clear();

	m_isPion = false;

	int Nhit1 = 0;
	int Nhit2 = 0;
	int Nhit3 = 0;

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double weight = 0.0;
	double sumweight = 0.0;

	std::vector<ThreeVector> positions;

	LCCollection *pLCCollection = 0;

	try
	{
		pLCCollection = pLCEvent->getCollection( m_sdhcalCollectionName );
	}
	catch( DataNotAvailableException &e )
	{
		streamlog_out(ERROR) << "PionSelector::processEvent(evt): Collection '" + m_sdhcalCollectionName + "' not found !" << std::endl;
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	UTIL::CellIDDecoder<EVENT::CalorimeterHit> cellIdDecoder(pLCCollection);

	for( unsigned int i=0 ; i<pLCCollection->getNumberOfElements() ; i++ )
	{
		if( pLCCollection->getTypeName() != EVENT::LCIO::CALORIMETERHIT )
			return dqm4hep::STATUS_CODE_INVALID_PARAMETER;

		EVENT::CalorimeterHit *pCaloHit = static_cast<EVENT::CalorimeterHit*> ( pLCCollection->getElementAt(i) );
		int K = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(2).c_str() ];

		if(K > m_nbOfLayers - 1)
			continue;

		m_caloHitCollection.push_back(pCaloHit);

		positions.push_back( ThreeVector(
				pCaloHit->getPosition()[0],
				pCaloHit->getPosition()[1],
				pCaloHit->getPosition()[2] ) );
	}

	UTIL::CellIDDecoder<EVENT::CalorimeterHit>::setDefaultEncoding(m_decoderString);
	std::sort(m_caloHitCollection.begin(), m_caloHitCollection.end() , PionSelector::sortByLayer );

	// first cut on NHit
	if( m_caloHitCollection.size() < m_nHitCut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
	{
		EVENT::CalorimeterHit *pCaloHit = m_caloHitCollection.at(h);

		float fThr = pCaloHit->getEnergy();
		int I = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(0).c_str() ];
		int J = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(1).c_str() ];
		int K = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(2).c_str() ];

		if( I < 5 || I > 95 || J < 5 || J > 95 )
			m_nHitsInEdge ++;

		// the touched layers
		m_nlayers.at(K)++;

		// number of hits for each threshold + weight for cog computation
		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() ) {
			weight = 10.0;
			Nhit1++;
		}
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() ) {
			weight = 5.0;
			Nhit2++;
		}
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() ) {
			weight = 1.0;
			Nhit3++;
		}

		x += weight*I;
		y += weight*J;
		z += weight*K;
		sumweight += weight;
	}

	// get the number of touched layers in the event
	int nbOfTouchedLayers = 0;
	for( unsigned int l=0 ; l<m_nlayers.size() ; l++ )
		if( m_nlayers.at(l) != 0 )
			nbOfTouchedLayers++;

	m_nHit.push_back( m_caloHitCollection.size() );
	m_nHit.push_back( Nhit1 );
	m_nHit.push_back( Nhit2 );
	m_nHit.push_back( Nhit3 );

	// debug check
	assert( (Nhit1+Nhit2+Nhit3) == m_nHit.at(0) );

	m_cog.push_back( x/sumweight );
	m_cog.push_back( y/sumweight );
	m_cog.push_back( z/sumweight );

	m_fractalDimension = this->getFractalDimension();

	// calculate the shower starting layer
	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
	{
		EVENT::CalorimeterHit *pCaloHit = m_caloHitCollection.at(h);

		float fThr = pCaloHit->getEnergy();
		int I = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(0).c_str() ];
		int J = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(1).c_str() ];
		int K = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(2).c_str() ];

		if( fabs( I - m_cog.at(0) ) > 5 || fabs( J - m_cog.at(1) ) > 5 )
			continue;

		int count = 1;
		std::vector<int> count2(3,0);

		for( unsigned int h2=0 ; h2<m_caloHitCollection.size() ; h2++ )
		{
			EVENT::CalorimeterHit *pCaloHit2 = m_caloHitCollection.at(h2);

			int I2 = cellIdDecoder(pCaloHit2)[ m_ijkEncoding.at(0).c_str() ];
			int J2 = cellIdDecoder(pCaloHit2)[ m_ijkEncoding.at(1).c_str() ];
			int K2 = cellIdDecoder(pCaloHit2)[ m_ijkEncoding.at(2).c_str() ];

			if( h == h2 )
				continue;

			if( fabs( I2 - m_cog.at(0) ) < 5 && fabs( J2 - m_cog.at(1) ) < 5 )
				count++;

			if( fabs( I2 - m_cog.at(0) ) > 5 || fabs( J2 - m_cog.at(1) ) > 5 )
				continue;

			if( K2 == K+1 ) count2.at(0) ++;
			if( K2 == K+2 ) count2.at(1) ++;
			if( K2 == K+3 ) count2.at(2) ++;
		}

		if( count <= 4 )
			continue;

		if( count2.at(0) >= 4 && count2.at(1) >= 4 && count2.at(2) >= 4 ) {

			m_showerStartingLayer = K;
			break;
		}
	}

	sumweight = 0.0;
	int count = 0;

	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ ) {

		EVENT::CalorimeterHit *pCaloHit = m_caloHitCollection.at(h);

		int I = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(0).c_str() ];
		int J = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(1).c_str() ];
		int K = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(2).c_str() ];

		if( K >= m_showerStartingLayer ) {

			m_radius += (m_cog.at(0) - I)*(m_cog.at(0) - I) + (m_cog.at(1) - J)*(m_cog.at(1) - J);
			sumweight ++;
		}

		if( ( I - m_cog.at(0) ) < 3 && ( J - m_cog.at(1) ) < 3 )
			count ++;
	}

	m_nHitsInCentralCells = count;

	if( sumweight != 0 )
		std::sqrt( m_radius /= sumweight );


	for( int Kiter = m_showerStartingLayer+1 ; Kiter < m_showerStartingLayer+8  ; Kiter++ )
	{

		int count = 0;

		for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
		{
			EVENT::CalorimeterHit *pCaloHit = m_caloHitCollection.at(h);

			float fThr = pCaloHit->getEnergy();
			int I = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(0).c_str() ];
			int J = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(1).c_str() ];
			int K = cellIdDecoder(pCaloHit)[ m_ijkEncoding.at(2).c_str() ];

			if( K == Kiter && ( I - m_cog.at(0) ) < 10 && ( J - m_cog.at(1) ) < 10 ) {
				count++;
				break;
			}
		}
		if( count == 0 )
			m_nHolesAfterStartingPoint ++;
	}

	bool isNeutral = true;

	for(unsigned int k=0 ; k<m_neutralFirstLayerCut ; k++)
	{
		if(m_nlayers.at(k) != 0)
		{
			isNeutral = false;
			break;
		}
	}

	if(isNeutral)
		return dqm4hep::STATUS_CODE_SUCCESS;

	if( (double)m_nHit.at(0)/(double)nbOfTouchedLayers < m_nHitOverNLayerCut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	if( m_radius / m_cog.at(2) < m_radiusOverCog2Cut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	if(m_showerStartingLayer < m_showerStartingLayerCut)
		return dqm4hep::STATUS_CODE_SUCCESS;

	if(m_showerStartingLayer < m_electronStartingLayerCut && nbOfTouchedLayers < m_nTouchedLayersCut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	if( m_nHolesAfterStartingPoint > m_nHolesCut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	if( m_nHitsInEdge / m_nHit.at(0) > m_nHitEdgePercentCut )
		return dqm4hep::STATUS_CODE_SUCCESS;

	// cog too fare from center, then cut
	if(m_cog.at(0) < m_barycenterPositionCut || m_nbOfCellSizeX - m_cog.at(0) < m_barycenterPositionCut
	|| m_cog.at(1) < m_barycenterPositionCut || m_nbOfCellSizeY - m_cog.at(1) < m_barycenterPositionCut)
		return dqm4hep::STATUS_CODE_SUCCESS;

	std::vector<int> weights(positions.size(), 1);
	Linear3DFit *pFitter = new Linear3DFit(positions, weights);
	pFitter->Fit();

	ThreeVector px(-1, 0, pFitter->GetFitParameters()[1]);
	ThreeVector py(0, -1, pFitter->GetFitParameters()[3]);
	float cosTheta = px.Cross(py).CosTheta();

	delete pFitter;

	if(cosTheta < m_cosThetaCut)
		return dqm4hep::STATUS_CODE_SUCCESS;

	if(!this->isSingleParticle())
		return dqm4hep::STATUS_CODE_SUCCESS;

	m_isPion = true;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

double PionSelector::getFractalDimension()
{
	int vec[] = {2,3,4,6,8,12,16};
	float f3D = 0;

	for(int i=0; i<7; i++)
	{
		int Ncube = nbOfHitsInCube(vec[i]);

		if(Ncube >= m_nHit.at(0))
			return 0;

		f3D += std::log(float(m_nHit.at(0))/Ncube)/std::log(vec[i]);
	}
	return f3D/7.0;
}

//-------------------------------------------------------------------------------------------------

int PionSelector::nbOfHitsInCube( int cubeSize )
{
	std::vector<int> keys;
	int ncube = 0;
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder( m_decoderString );

	for( unsigned int h=0 ; h<m_caloHitCollection.size() ; h++ )
	{
		EVENT::CalorimeterHit *pCaloHit = m_caloHitCollection.at(h);

		int newI = idDecoder(pCaloHit)[m_ijkEncoding.at(0).c_str()]/(cubeSize+1);
		int newJ = idDecoder(pCaloHit)[m_ijkEncoding.at(1).c_str()]/(cubeSize+1);
		int newK = idDecoder(pCaloHit)[m_ijkEncoding.at(2).c_str()]/(cubeSize+1);

		int key = IJKToKey( newI , newJ , newK );

		if( std::find( keys.begin() , keys.end() , key ) != keys.end() )
			continue;

		ncube++;
		keys.push_back( key );
	}
	return ncube;
}

//-------------------------------------------------------------------------------------------------

bool PionSelector::isSingleParticle()
{
	std::vector<double> layerCogX(m_nbOfLayers, 0);
	std::vector<double> layerCogY(m_nbOfLayers, 0);
	std::vector<double> layerRMSX(m_nbOfLayers, 0);
	std::vector<double> layerRMSY(m_nbOfLayers, 0);
	std::vector<double> weightSum(m_nbOfLayers, 0);

	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(m_decoderString);

	int firstLayersCut = 5;
	double largeRMSCut = 5.0;
	int largeRMSCounterCut = 4;
	std::vector<EVENT::CalorimeterHit*> firstCaloHitCollection;

	for(unsigned int i=0 ; i<m_caloHitCollection.size() ; i++)
	{
		EVENT::CalorimeterHit *caloHit = m_caloHitCollection.at(i);
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		if(K < firstLayersCut)
			firstCaloHitCollection.push_back(caloHit);
	}

	// compute the cog for each layer
	for(unsigned int i=0 ; i<firstCaloHitCollection.size() ; i++)
	{
		EVENT::CalorimeterHit *caloHit = firstCaloHitCollection.at(i);

		float fThr = caloHit->getEnergy();

		int I = idDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = idDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		double weight = 0.;

		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() )
			weight = 1.0;
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() )
			weight = 5.0;
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() )
			weight = 10.0;

		layerCogX.at(K) += I*weight;
		layerCogY.at(K) += J*weight;
		weightSum.at(K) += weight;
	}

	// normalize the cog
	for(unsigned int i=0 ; i<layerCogX.size() ; i++)
	{
		if(weightSum.at(i) != 0)
		{
			layerCogX.at(i) /= double(weightSum.at(i));
			layerCogY.at(i) /= double(weightSum.at(i));
		}
	}


	// compute the rms for each layer
	for(unsigned int i=0 ; i<firstCaloHitCollection.size() ; i++)
	{
		EVENT::CalorimeterHit *caloHit = firstCaloHitCollection.at(i);

		float fThr = caloHit->getEnergy();

		int I = idDecoder(caloHit)[m_ijkEncoding.at(0).c_str()];
		int J = idDecoder(caloHit)[m_ijkEncoding.at(1).c_str()];
		int K = idDecoder(caloHit)[m_ijkEncoding.at(2).c_str()];

		double weight = 0.;

		if( (fThr-1.f) < std::numeric_limits<float>::epsilon() )
			weight = 1.0;
		else if( (fThr-2.0)< std::numeric_limits<float>::epsilon() )
			weight = 5.0;
		else if( (fThr-3.0)< std::numeric_limits<float>::epsilon() )
			weight = 10.0;

		layerRMSX.at(K) += (I - layerCogX.at(K))*(I - layerCogX.at(K))*weight;
		layerRMSY.at(K) += (J - layerCogY.at(K))*(J - layerCogY.at(K))*weight;
	}

	int largeRMSCounter = 0;

	// normalize the rms
	for(unsigned int i=0 ; i<layerCogX.size() ; i++)
	{
		if(weightSum.at(i) != 0)
		{
			layerRMSX.at(i) = std::sqrt(layerRMSX.at(i)/double(weightSum.at(i)));
			layerRMSY.at(i) = std::sqrt(layerRMSY.at(i)/double(weightSum.at(i)));

			if(layerRMSX.at(i) > m_largeRMSCut || layerRMSY.at(i) > m_largeRMSCut)
				largeRMSCounter++;
		}
	}

	return largeRMSCounter < largeRMSCounterCut ? true : false;
}

} 
