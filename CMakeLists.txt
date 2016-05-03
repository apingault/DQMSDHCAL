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

###############################
# Require C++11
INCLUDE( CheckCXXCompilerFlag )

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

IF( COMPILER_SUPPORTS_CXX11 )
	LIST( APPEND CMAKE_CXX_FLAGS "-std=c++11" )
ELSE()
	MESSAGE( SEND_ERROR "DQMCore requires C++11 support. Please upgrade your compiler !" )
ENDIF()

#-- dependencies --
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# ----- Find DQMCore package -----
FIND_PACKAGE( DQMCore         REQUIRED )
LIST( APPEND CMAKE_MODULE_PATH ${DQMCore_DIR}/cmake )

INCLUDE( default_settings )

# ----- Find DQM4ILC package -----
FIND_PACKAGE( DQM4ILC         REQUIRED )

# ----- Find Trivent package -----
FIND_PACKAGE( Trivent         REQUIRED )

# ----- Find CaloSoftWare package -----
FIND_PACKAGE( CaloSoftWare    REQUIRED )

FOREACH( pkg DQMCore DQM4ILC Trivent CaloSoftWare )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR} )

ADD_SUBDIRECTORY( source )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
