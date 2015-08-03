#################################################
# cmake file for building DQMSDHCAL package
# @author Eté Rémi
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# -- project name --
PROJECT( DQMSDHCAL )

SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

#-- dependencies --
LIST( APPEND CMAKE_MODULE_PATH $ENV{DQM4HEP_DIR}/cmake )
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# load default settings from dqm4hep
INCLUDE( default_settings )

FIND_PACKAGE( DQM4HEP REQUIRED )

FOREACH( pkg DQM4HEP )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

ADD_SUBDIRECTORY( ./source )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
