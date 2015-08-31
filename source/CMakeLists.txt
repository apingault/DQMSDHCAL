#######################################################
# cmake file for building DQMSDHCAL package
# @author Eté Rémi, Victor Buridon
# Copyright (c) CNRS / IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )
INCLUDE_DIRECTORIES( include/dqmsdhcal/asicAnalysisTools )


# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")

# -------------------------------------------------
# build the global library

FILE (
	GLOB_RECURSE
	SRC_FILES
	src/daq/*
	src/streamout/*
	src/trivent/*
	src/asicAnalysisTools/*
	src/services/*
)

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${SRC_FILES} )

INSTALL(
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build shm reader executable

ADD_EXECUTABLE( dqmsdhcal_start_shm_reader main/dqmsdhcal_start_shm_reader.cc )
TARGET_LINK_LIBRARIES( dqmsdhcal_start_shm_reader ${PROJECT_NAME} )

INSTALL(
    TARGETS dqmsdhcal_start_shm_reader
    RUNTIME DESTINATION bin
)

# -------------------------------------------------
# build raw data converter executable

ADD_EXECUTABLE( dqmsdhcal_start_raw_data_converter main/dqmsdhcal_start_raw_data_converter.cc )
TARGET_LINK_LIBRARIES( dqmsdhcal_start_raw_data_converter ${PROJECT_NAME} )

INSTALL(
    TARGETS dqmsdhcal_start_raw_data_converter
    RUNTIME DESTINATION bin
)

# -------------------------------------------------
# build test raw data module library

ADD_SHARED_LIBRARY( TestRawDataModule src/modules/TestRawDataModule.cc )
TARGET_LINK_LIBRARIES( TestRawDataModule ${PROJECT_NAME} )

INSTALL(
    TARGETS TestRawDataModule
    LIBRARY DESTINATION lib
)


# -------------------------------------------------
# build asic analysis module library

ADD_SHARED_LIBRARY( AsicAnalysisModule src/modules/AsicAnalysisModule.cc )
TARGET_LINK_LIBRARIES( AsicAnalysisModule ${PROJECT_NAME} )

INSTALL(
    TARGETS AsicAnalysisModule
    LIBRARY DESTINATION lib
)

