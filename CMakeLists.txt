#################################################
# cmake file for building DQMSDHCAL package
# @author Eté Rémi
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# project name
PROJECT( DQMSDHCAL )

SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# --------------- ilcsoft dependencies -----------
LIST( APPEND CMAKE_MODULE_PATH $ILCSOFT )
LIST( APPEND CMAKE_MODULE_PATH $ENV{LCDQM_DIR} )
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# load default settings
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( LCDQM REQUIRED )

  INCLUDE_DIRECTORIES( ${LCDQM_INCLUDE_DIRS} )
  LINK_LIBRARIES( ${LCDQM_LIBRARIES} )
  ADD_DEFINITIONS ( ${LCDQM_DEFINITIONS} )

ADD_SUBDIRECTORY( ./source )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
