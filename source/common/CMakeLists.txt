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
# build shm publisher executable

ADD_EXECUTABLE( dqmsdhcal_start_shm_publisher main/dqmsdhcal_start_shm_publisher.cc )
TARGET_LINK_LIBRARIES( dqmsdhcal_start_shm_publisher ${PROJECT_NAME}Common )

INSTALL(
    TARGETS dqmsdhcal_start_shm_publisher
    RUNTIME DESTINATION bin
)


# -------------------------------------------------
# build dump db geometry executable

ADD_EXECUTABLE( dqmsdhcal_dump_db_geometry main/dqmsdhcal_dump_db_geometry.cc )
TARGET_LINK_LIBRARIES( dqmsdhcal_dump_db_geometry ${PROJECT_NAME}Common )

INSTALL(
    TARGETS dqmsdhcal_dump_db_geometry
    RUNTIME DESTINATION bin
)


