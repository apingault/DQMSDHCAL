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

#include <algorithm>

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

const unsigned short SDHCALElectronicsMapping::m_channelTable[] =
  { 12, 13, 23, 24, 39, 40, 51, 52,
	10, 11, 22, 25, 38, 41, 53, 54,
	8, 9, 21, 26, 37, 50, 55, 56,
	7, 6, 14, 27, 36, 49, 57, 58,
	5, 4, 15, 28, 35, 48, 60, 59,
	3, 2, 16, 29, 34, 47, 62, 61,
	1, 0, 19, 30, 33, 44, 45, 63,
	17, 18, 20, 31, 32, 42, 43, 46 };

//-------------------------------------------------------------------------------------------------

// (J Axis)
//	 	|04|05|12|13|
//	 	|03|06|11|14|  ASIC MAPPING
//	 	|02|07|10|15|  WITHIN A DIF
//	 	|01|08|09|16|
//                     (I Axis)  ----->

/// the asic table given a channel index in i and j cells direction
const unsigned short SDHCALElectronicsMapping::m_asicTable[] =
  { 4,  5,  12, 13, 20, 21, 28, 29, 36, 37, 44, 45,
	3,  6,  11, 14, 19, 22, 27, 30, 35, 38, 43, 46,
	2,  7,  10, 15, 18, 23, 26, 31, 34, 39, 42, 47,
	1,  8,  9,  16, 17, 24, 25, 32, 33, 40, 41, 48 };

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
    m_layerThickness(0),
    m_rotateAxes(0),
		m_moduleLogStr("[SDHCALElectronicsMapping]")
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

	// first of all, find dif id
	unsigned int difShiftY = ((cell.m_jCell-1)/32)*32;

	Geometry::const_iterator iter = m_geometry.find(cell.m_layer);

	if( m_geometry.end() == iter )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - cellToElectronics : Couldn't find layer id " << cell.m_layer << " in SDHCALElectronicsMapping maps" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	DifMapping::const_iterator iter2 = std::find_if( iter->second.m_difList.begin() , iter->second.m_difList.end(),
			[&] (const DifMapping::value_type &value) {
		return value.second.m_shiftY == difShiftY;
	});

	if( iter->second.m_difList.end() == iter2 )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - cellToElectronics : Couldn't find dif with shift " << difShiftY << " in SDHCALElectronicsMapping maps" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	electronics.m_difId = iter2->first;

	unsigned int indexAsicI = (cell.m_iCell-1)/8;
	unsigned int indexAsicJ = ((cell.m_jCell-1)%32)/8;

	if(indexAsicI >= 12 || indexAsicJ >=4 )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger, m_moduleLogStr << " - cellToElectronics : Wrong asic index for cell to electronics conversion ! - indexAsicI: " << indexAsicI << " (expected <12) \t indexAsicJ: "  << indexAsicJ << " (expected <4)" )
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	electronics.m_asicId = m_asicTable[ indexAsicI + 12*indexAsicJ ];

	unsigned int indexChannelI = ((cell.m_iCell-1)%8);
	unsigned int indexChannelJ = ((cell.m_jCell-1)%8);

	if(indexChannelI >= 8 || indexChannelJ >=8 )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger, m_moduleLogStr << " - cellToElectronics : Wrong channel index for cell to electronics conversion ! - indexChannelI: " << indexChannelI << "(expected <8) \t indexChannelJ: "  << indexChannelJ << "(expected <8)" )
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	electronics.m_channelId = m_channelTable[ indexChannelI + 8*indexChannelJ ];

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

	DifMapping::const_iterator iter = m_difMapping.find(electronics.m_difId);

	if( m_difMapping.end() == iter )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - electronicstoCell : Couldn't find dif id " << electronics.m_difId << " in SDHCALElectronicsMapping map" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	// 1 because pad ids start from 1
	cell.m_iCell = 1 + m_channelToIMapping[electronics.m_channelId] + m_asicToChannelShiftI[electronics.m_asicId];
	cell.m_jCell = 32 - (m_channelToJMapping[electronics.m_channelId] + m_asicToChannelShiftJ[electronics.m_asicId]) + iter->second.m_shiftY;
	cell.m_layer = iter->second.m_layerId;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::positionToCell(const dqm4hep::DQMCartesianVector &position, dqm4hep::DQMElectronicsMapping::Cell &cell)
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	// find layer
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, this->findClosestLayer(position - m_cellReferencePosition, cell.m_layer));

	// compute i and j expected values
	float iCellFloat(0.f);
	float jCellFloat(0.f);

	if(m_rotateAxes)
	{
		iCellFloat = (-1.f * position.getY() - m_cellReferencePosition.getY()) / m_cellSize1;
		jCellFloat = (-1.f * position.getX() - m_cellReferencePosition.getX()) / m_cellSize0;
	}
	else
	{
		iCellFloat = (position.getX() - m_cellReferencePosition.getX()) / m_cellSize0;
		jCellFloat = (position.getY() - m_cellReferencePosition.getY()) / m_cellSize1;
	}

	if( iCellFloat < 0.f || jCellFloat < 0.f )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger, m_moduleLogStr << " - positionToCell : Wrong cell index for position to cell conversion ! - iCellFloat: " << iCellFloat << "(expected >=0) \t jCellFloat: "  << jCellFloat << "(expected >=0)" )
		return dqm4hep::STATUS_CODE_FAILURE;
	}

	cell.m_iCell = round(iCellFloat);
	cell.m_jCell = round(jCellFloat);

	if(cell.m_iCell < 1 || cell.m_iCell > 96
	|| cell.m_jCell < 1 || cell.m_jCell > 96)
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger, m_moduleLogStr << " - positionToCell : Wrong cell index for position to cell conversion ! - cell.m_iCell: " << cell.m_iCell << " (expected >=1 <=96) \t cell.m_jCell: "  << cell.m_jCell << " (expected >=1 <=96)" )
		return dqm4hep::STATUS_CODE_FAILURE;
	} 
	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::cellToPosition(const dqm4hep::DQMElectronicsMapping::Cell &cell, dqm4hep::DQMCartesianVector &position)
{
	if( ! m_isInitialized )
		return dqm4hep::STATUS_CODE_NOT_INITIALIZED;

	Geometry::iterator iter = m_geometry.find(cell.m_layer);

	if( m_geometry.end() == iter )
	{
		LOG4CXX_ERROR( dqm4hep::dqmMainLogger , m_moduleLogStr << " - cellToPosition : Couldn't find layer id " << cell.m_layer << " in SDHCALElectronicsMapping maps" );
		return dqm4hep::STATUS_CODE_NOT_FOUND;
	}

	float x(0.f);
	float y(0.f);

	if(m_rotateAxes)
	{
		x = -1.f * cell.m_jCell * m_cellSize1 + m_cellReferencePosition.getX();
		y = -1.f * cell.m_iCell * m_cellSize0 + m_cellReferencePosition.getY();
	}
	else
	{
		x = cell.m_iCell * m_cellSize0 + m_cellReferencePosition.getX();
		y = cell.m_jCell * m_cellSize1 + m_cellReferencePosition.getY();
	}

	const float z = iter->second.m_z0 + m_cellReferencePosition.getZ();

	position.setValues(x, y, z);

	return dqm4hep::STATUS_CODE_SUCCESS;
}

//-------------------------------------------------------------------------------------------------

dqm4hep::StatusCode SDHCALElectronicsMapping::findClosestLayer(const dqm4hep::DQMCartesianVector &position, unsigned int &layer)
{
	float closestZ(std::numeric_limits<float>::max());
	layer = std::numeric_limits<unsigned int>::max();

	for(Geometry::iterator iter = m_geometry.begin(), endIter = m_geometry.end() ;
			endIter != iter ; ++iter)
	{
		if(closestZ > fabs(iter->second.m_z0-position.getZ()))
		{
			closestZ = fabs(iter->second.m_z0-position.getZ());
			layer = iter->first;
		}
	}

	if( std::numeric_limits<unsigned int>::max() == layer )
		return dqm4hep::STATUS_CODE_NOT_FOUND;

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

	std::string geometryFileName = "sdhcalGeometry.xml";
	RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"GeometryFileName", geometryFileName));

	if( readFromDB )
	{
		std::string host, user, password, database, beamTest;

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
				"Host", host));

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
				"User", user));

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
				"Password", password));

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
				"Database", database));

		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
				"BeamTest", beamTest));

		GeometryDBInterface interface;
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, interface.connect( host , user , password , database ));
		RETURN_RESULT_IF(dqm4hep::STATUS_CODE_SUCCESS, !=, interface.dumpGeometry( geometryFileName , beamTest ));
	}

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellReferencePosition", m_cellReferencePosition));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellSize0", m_cellSize0));

	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"CellSize1", m_cellSize1));

	dqm4hep::UIntVector layerMask;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"LayerMask", layerMask));

	dqm4hep::UIntVector difMask;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValues(xmlHandle,
			"DifMask", difMask));

	m_rotateAxes = true;
	RETURN_RESULT_IF_AND_IF(dqm4hep::STATUS_CODE_SUCCESS, dqm4hep::STATUS_CODE_NOT_FOUND, !=, dqm4hep::DQMXmlHelper::readParameterValue(xmlHandle,
			"RotateAxes", m_rotateAxes));

	GeometryXmlIO reader;
	RETURN_RESULT_IF( dqm4hep::STATUS_CODE_SUCCESS, !=, reader.load( geometryFileName , m_geometry, layerMask , difMask ) );

	for(Geometry::iterator iter = m_geometry.begin(), endIter = m_geometry.end() ;
			endIter != iter ; ++iter)
	{
		for(DifMapping::iterator difIter = iter->second.m_difList.begin(), difEndIter = iter->second.m_difList.end() ;
				difEndIter != difIter ; ++difIter)
		{
			m_difMapping[ difIter->first ] = difIter->second;
		}
	}

	m_isInitialized = true;

	return dqm4hep::STATUS_CODE_SUCCESS;
}

} 

