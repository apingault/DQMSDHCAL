  /// \file ECalHelper.h
/*
 *
 * ECalHelper.h header template automatically generated by a class generator
 * Creation date : sam. juin 11 2016
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


#ifndef ECALHELPER_H
#define ECALHELPER_H

#include "dqm4hep/DQM4HEP.h"
#include "dqm4hep/DQMXmlHelper.h"

class TH2;

namespace dqm_sdhcal
{

/** 
 * @brief ECalCalibration class
 */ 
class ECalHelper
{
public:
	/** Constructor
	 */
	ECalHelper();

	/** Destructor
	 */
	~ECalHelper();

	/** Get the number of layers
	 */
	unsigned int getNLayers() const;

	/** Get the ordered layer in the beam direction from the hardware layer
	 */
	dqm4hep::StatusCode getLayer(int hwLayer, unsigned int &layer) const;

	/** Get the hardware layer in the beam direction from the hardware layer
	 */
	dqm4hep::StatusCode getHwLayer(unsigned int layer, int &hwlayer) const;

	/** Get the slab id from the layer id (Not the hardware one. Use conversion for that)
	 */
	dqm4hep::StatusCode getSlab(unsigned int layer, unsigned int &slab) const;

	/** Get pedestal value given :
	 *    - a layer id
	 *    - an asic id
	 *    - a channel id
	 *    - a column id
	 */
	dqm4hep::StatusCode getPedestal(unsigned int layer, unsigned int asic, unsigned int channel,
			unsigned int column, int &pedestal) const;

	/** Get the corrected ADC count given :
	 *    - a layer id
	 *    - an asic id
	 *    - a channel id
	 *    - a column id
	 */
	dqm4hep::StatusCode getCorrectedADCCount(unsigned int layer, unsigned int asic, unsigned int channel,
			unsigned int column, int adcCount, int &correctedAdcCount) const;

	/** Get the calibrated electromagnetic energy given :
	 *    - a layer id
	 *    - an asic id
	 *    - a channel id
	 *    - a column id
	 */
	dqm4hep::StatusCode getCalibratedEnergy(unsigned int layer, unsigned int asic, unsigned int channel,
			unsigned int column, int adcCount, float &calibratedEnergy) const;

	/** Read settings from xml handle
	 */
	dqm4hep::StatusCode readSettings(const dqm4hep::TiXmlHandle &handle);

private:
	typedef std::map<unsigned int, TH2*>                  ColumnToPedestalMap;
	typedef std::map<unsigned int, ColumnToPedestalMap>   SlabPedestalMap;
	typedef std::map<int, unsigned int>                   IntToUIntMap;
	typedef std::map<unsigned int, int>                   UIntToIntMap;

	unsigned int                    m_nColumns;
	unsigned int                    m_correctedPedestalCut;
	unsigned int                    m_mipCalibrationFactor;
	int                             m_pedestalShift;
	float                           m_mipEquivalentEnergy;
	UIntToIntMap                    m_layerToHwLayerMap;
	IntToUIntMap                    m_hwLayerToLayerMap;
	IntToUIntMap                    m_hwLayerToSlabMap;
	UIntToIntMap                    m_slabToHwLayerMap;
	SlabPedestalMap                 m_slabPedestalMap;
};

} 

#endif  //  ECALHELPER_H