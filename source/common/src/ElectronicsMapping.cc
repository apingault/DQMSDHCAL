  /// \file ElectronicsMapping.cc
/*
 *
 * ElectronicsMapping.cc source template automatically generated by a class generator
 * Creation date : ven. avr. 8 2016
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


#include "ElectronicsMapping.h"

#include "dqm4hep/DQMPlugin.h"
#include "dqm4hep/DQMXmlHelper.h"

namespace dqm_sdhcal
{

// declare dqm plugin
DQM_PLUGIN_DECL(SDHCALElectronicsMapping, "SDHCALElectronicsMapping")

//      (J Axis)
// 7	|46|63|61|59|58|56|54|52|
// 6	|43|45|62|60|57|55|53|51|
// 5	|42|44|47|48|49|50|41|40|
// 4	|32|33|34|35|36|37|38|39|
// 3	|31|30|29|28|27|26|25|24|  	     TOP VIEW (ASICs SIDE)
// 2	|20|19|16|15|14|21|22|23|
// 1	|18|00|02|04|06|09|11|13|
// 0	|17|01|03|05|07|08|10|12|
//       0  1  2  3  4  5  6  7    (I Axis)  ----->
//				 |	  |
//				 |DIF |
//				 |____|

/// table of channels within an asic
/// channelID = channelTable[ i%96/12 + (j%96/12 )*4 ]
const unsigned short SDHCALElectronicsMapping::m_channelTable[] =
  {17,  1,  3,  5,  7,  8, 10, 12,
   18,  0,  2,  4,  6,  9, 11, 13,
   20, 19, 16, 15, 14, 21, 22, 23,
   31, 30, 29, 28, 27, 26, 25, 24,
   32, 33, 34, 35, 36, 37, 38, 39,
   42, 44, 47, 48, 49, 50, 41, 40,
   43, 45, 62, 60, 57, 55, 53, 51,
   46, 63, 61, 59, 58, 56, 54, 52 };

//-------------------------------------------------------------------------------------------------

// (J Axis)
//	 	|04|05|12|13|
//	 	|03|06|11|14|  ASIC MAPPING
//	 	|02|07|10|15|  WITHIN A DIF
//	 	|01|08|09|16|
//                     (I Axis)  ----->

/// the asic table given a channel index in i and j cells direction
/// asicID = asicTable[ i%32/8 + (j%32/8 )*4 ]
const unsigned short SDHCALElectronicsMapping::m_asicTable[] =
{ 1,  8,  9,  16,
  2,  7,  10, 15,
  3,  6,  11, 14,
  4,  5,  12, 13 };

//-------------------------------------------------------------------------------------------------

/// projection of channel id to i cell id within an asic
const unsigned short SDHCALElectronicsMapping::m_channelToIMapping[] =
  {1,0,1,0,1,0,1,0,0,1,0,1,0,1,2,2,
   2,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,
   4,4,4,4,4,4,4,4,5,5,5,6,5,6,7,5,
   5,5,5,6,7,6,7,6,7,6,7,7,6,7,6,7};

//-------------------------------------------------------------------------------------------------

/// projection of channel id to j cell id within an asic
const unsigned short SDHCALElectronicsMapping::m_channelToJMapping[] =
  {1,1,2,2,3,3,4,4,5,5,6,6,7,7,4,3,
   2,0,0,1,0,5,6,7,7,6,5,4,3,2,1,0,
   0,1,2,3,4,5,6,7,7,6,0,0,1,1,0,2,
   3,4,5,7,7,6,6,5,5,4,4,3,3,2,2,1};

//-------------------------------------------------------------------------------------------------

/// conversion between the asic id to channel shift in i cell direction
const unsigned short SDHCALElectronicsMapping::m_asicToChannelShiftI[]=
  {	0,	0,	0,	0,	0,	8,	8,
	8,	8,	16,	16,	16,	16,	24,
	24,	24,	24,	32,	32,	32,	32,
	40,	40,	40,	40,	48,	48,	48,
	48,	56,	56,	56,	56,	64,	64,
	64,	64,	72,	72,	72,	72,	80,
	80,	80,	80,	88,	88,	88,	88  };

//-------------------------------------------------------------------------------------------------

/// conversion between the asic id to channel shift in j cell direction
const unsigned short SDHCALElectronicsMapping::m_asicToChannelShiftJ[]=
  {	0,	0,	8,	16,	24,	24,	16,
	8,	0,	0,	8,	16,	24,	24,
	16,	8,	0,	0,	8,	16,	24,
	24,	16,	8,	0,	0,	8,	16,
	24,	24,	16,	8,	0,	0,	8,
	16,	24,	24,	16,	8,	0,	0,
	8,	16,	24,	24,	16,	8,	0  };

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

SDHCALElectronicsMapping::SDHCALElectronicsMapping() :
		m_isInitialized(false),
		m_cellReferencePosition(0.f, 0.f, 0.f),
		m_cellSize0(10.408f),
		m_cellSize1(10.408f),
		m_layerThickness(26.73f),
		m_globalDifShiftY(32)
{
	/* nop */
}

//-------------------------------------------------------------------------------------------------

SDHCALElectronicsMapping::~SDHCALElectronicsMapping() 
{
	/* nop */
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::cellToElectronics( const dqm4hep::DQMElectronicsMapping::Cell &cell, dqm4hep::DQMElectronicsMapping::Electronics &electronics )
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	electronics.m_difId = 0;
	electronics.m_asicId = 0;
	electronics.m_channelId = 0;

	unsigned int indexAsicI = (cell.m_iCell%32)/8;
	unsigned int indexAsicJ = (cell.m_jCell%32)/8;

	electronics.m_asicId = m_asicTable[ indexAsicI + 4*indexAsicJ ];

	unsigned int indexChannelI = (cell.m_iCell%96/12);
	unsigned int indexChannelJ = (cell.m_jCell%96/12);

	electronics.m_channelId = m_channelTable[ indexChannelI + 4*indexChannelJ ];

	unsigned int difShiftY = cell.m_jCell/32;

	LayerDifGeoMap::const_iterator iter = m_layerDifGeoMap.find(cell.m_layer);

	if( m_layerDifGeoMap.end() == iter )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Couldn't find layer id " << cell.m_layer << " in SDHCALElectronicsMapping maps" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	DifGeoMap::const_iterator iter2 = iter->second.find(difShiftY);

	if( iter->second.end() == iter2 )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Couldn't find dif with shift " << difShiftY << " in SDHCALElectronicsMapping maps" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	electronics.m_difId = iter2->second.m_difId;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::electronicstoCell( const dqm4hep::DQMElectronicsMapping::Electronics &electronics, dqm4hep::DQMElectronicsMapping::Cell &cell )
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	cell.m_iCell = 0;
	cell.m_jCell = 0;
	cell.m_layer = 0;

	DifGeoMap::const_iterator iter = m_difGeoMap.find(electronics.m_difId);

	if( m_difGeoMap.end() == iter )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Couldn't find dif id " << electronics.m_difId << " in SDHCALElectronicsMapping map" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	// 1 because pad ids start from 1
	cell.m_iCell = 1 + m_channelToIMapping[electronics.m_channelId] + m_asicToChannelShiftI[electronics.m_asicId];
	cell.m_jCell = 32 - m_channelToJMapping[electronics.m_channelId] + m_asicToChannelShiftJ[electronics.m_asicId] + iter->second.m_yShift;
	cell.m_layer = iter->second.m_layer;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::positionToCell(const dqm4hep::DQMCartesianVector &position, dqm4hep::DQMElectronicsMapping::Cell &cell)
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	const float iCellFloat = position.getX() - m_cellReferencePosition.getX() / m_cellSize0;
	const float jCellFloat = position.getY() - m_cellReferencePosition.getY() / m_cellSize1;
	const float layerFloat = position.getZ() - m_cellReferencePosition.getZ() / m_layerThickness;

	if( iCellFloat < 0.f || jCellFloat < 0.f || layerFloat < 0.f )
		return dqm4hep::STATUS_CODE_FAILURE;

	cell.m_iCell = round(iCellFloat);
	cell.m_jCell = round(jCellFloat);
	cell.m_layer = round(layerFloat);

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::cellToPosition(const dqm4hep::DQMElectronicsMapping::Cell &cell, dqm4hep::DQMCartesianVector &position)
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	const float x = cell.m_iCell * m_cellSize0 + m_cellReferencePosition.getX();
	const float y = cell.m_jCell * m_cellSize1 + m_cellReferencePosition.getY();
	const float z = cell.m_layer * m_layerThickness + m_cellReferencePosition.getZ();

	position.setValues(x, y, z);

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::readSettings(const dqm4hep::TiXmlHandle xmlHandle)
{
	if(m_isInitialized)
		return dqm4hep::STATUS_CODE_ALREADY_INITIALIZED;

	bool readFromDB = false;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"ReadFromDB", readFromDB));

	if( readFromDB )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , "Read DIF mapping from DB not yet implemented please set it to false !" );
		return dqm4hep::STATUS_CODE_INVALID_PARAMETER;
	}

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellReferencePosition", m_cellReferencePosition));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellSize0", m_cellSize0));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellSize1", m_cellSize1));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"LayerThickness", m_layerThickness));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"GlobalDifShiftY", m_globalDifShiftY));

	// require 3 difs per layer
	dqm4hep::UIntVector difList;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"DifList", difList, [] (const dqm4hep::UIntVector &vec) { return vec.size() % 3 == 0; }));

	for(unsigned int d=0 ; d<difList.size() ; d++)
	{
		unsigned int layer = d/3;
		unsigned int difId = difList.at(d);
		unsigned int difShiftY = m_globalDifShiftY * d%3;

		if( difId == 0 )
			continue;

		DifGeo difGeo;
		difGeo.m_layer = layer;
		difGeo.m_difId = difId;
		difGeo.m_yShift = difShiftY;

		m_difGeoMap[ difId ] = difGeo;
		m_layerDifGeoMap[ layer ][ difShiftY ] = difGeo;
	}

	m_isInitialized = true;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

} 
