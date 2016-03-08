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
	src/*
)

ADD_SHARED_LIBRARY( ${PROJECT_NAME}Common ${SRC_FILES} )

INSTALL(
	TARGETS ${PROJECT_NAME}Common
	LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build shm reader executable

ADD_EXECUTABLE( dqmsdhcal_start_shm_reader main/dqmsdhcal_start_shm_reader.cc )
TARGET_LINK_LIBRARIES( dqmsdhcal_start_shm_reader ${PROJECT_NAME}Common )

INSTALL(
    TARGETS dqmsdhcal_start_shm_reader
    RUNTIME DESTINATION bin
)


