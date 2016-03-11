#######################################################
# cmake file for building DQMSDHCAL package
# @author Eté Rémi, Victor Buridon
# Copyright (c) CNRS / IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( ../common/include )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")

# -------------------------------------------------
# build AsicAnalysisModule library

INCLUDE_DIRECTORIES( AsicAnalysisModule )
ADD_SHARED_LIBRARY( AsicAnalysisModule AsicAnalysisModule/AsicAnalysisModule.cc )
TARGET_LINK_LIBRARIES( AsicAnalysisModule ${PROJECT_NAME}Common )

INSTALL(
    TARGETS AsicAnalysisModule
    LIBRARY DESTINATION lib
)


# -------------------------------------------------
# build SlowControlModule library

INCLUDE_DIRECTORIES( SlowControlModule )
ADD_SHARED_LIBRARY( SlowControlModule SlowControlModule/SlowControlModule.cc )
TARGET_LINK_LIBRARIES( SlowControlModule ${PROJECT_NAME}Common )

INSTALL(
    TARGETS SlowControlModule
    LIBRARY DESTINATION lib
)

# -------------------------------------------------
# build EventDisplayModule library

INCLUDE_DIRECTORIES( EventDisplayModule )
ADD_SHARED_LIBRARY( EventDisplayModule EventDisplayModule/EventDisplayModule.cc )
TARGET_LINK_LIBRARIES( EventDisplayModule ${PROJECT_NAME}Common )

INSTALL(
    TARGETS EventDisplayModule
    LIBRARY DESTINATION lib
)



