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

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES streamlog REQUIRED )
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( DQMCore   REQUIRED )
FIND_PACKAGE( DQM4ILC   REQUIRED )
FIND_PACKAGE( LOG4CXX   REQUIRED )
FIND_PACKAGE( LCIO      REQUIRED )
FIND_PACKAGE( xdrstream REQUIRED )
FIND_PACKAGE( xdrlcio   REQUIRED )
FIND_PACKAGE( Marlin    REQUIRED )

FOREACH( pkg DQMCore DQM4ILC LCIO streamlog xdrstream xdrlcio LOG4CXX Marlin )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

ADD_SUBDIRECTORY( source )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
