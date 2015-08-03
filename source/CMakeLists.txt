#######################################################
# cmake file for building DQMSDHCAL package
# @author Eté Rémi, Victor Buridon
# Copyright (c) CNRS / IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")

# -------------------------------------------------
# build the global library

FILE (
	GLOB_RECURSE
	SRC_FILES
	src/daq/*
)

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${SRC_FILES} )

INSTALL(
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION lib
)

ADD_EXECUTABLE( start_sdhcal_shm_reader main/start_sdhcal_shm_reader.cc )
TARGET_LINK_LIBRARIES( start_sdhcal_shm_reader ${PROJECT_NAME} )

INSTALL(
    TARGETS start_sdhcal_shm_reader
    RUNTIME DESTINATION bin
)
