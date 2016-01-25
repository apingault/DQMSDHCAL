# CMakeList.txt for ExampleModule
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
# project name
PROJECT( ExampleModule )
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# append DQM4HEPConfig.cmake location to the default cmake module list
LIST( APPEND CMAKE_MODULE_PATH /opt/dqm4hep/cmake )
# load dqm4hep package
FIND_PACKAGE( DQMCore REQUIRED )
INCLUDE( default_settings )
# include the dqm4hep directories
INCLUDE_DIRECTORIES( ${DQMCore_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( /home/mirabito/SDHCAL/SDHCAL_common/include  )
INCLUDE_DIRECTORIES( /home/mirabito/SDHCAL/SDHCAL_Analysis/include  )
INCLUDE_DIRECTORIES(/usr/local/include )
INCLUDE_DIRECTORIES(/usr/local/cuda/include )
INCLUDE_DIRECTORIES(/usr/X11R6/include )
INCLUDE_DIRECTORIES(/usr/include/jsoncpp )
INCLUDE_DIRECTORIES(/usr/include/libxml2 )
INCLUDE_DIRECTORIES(/usr/local/include/eigen3 )
INCLUDE_DIRECTORIES(/usr/include/python2.7 )
# link lbDQM4HEP.so (.dylib) to our library
LINK_LIBRARIES( ${DQMCore_LIBRARIES} )
LINK_DIRECTORIES(/opt/dhcal/root/lib)
LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu )
LINK_DIRECTORIES(/usr/lib )
LINK_DIRECTORIES(/usr/local/lib )
LINK_DIRECTORIES(/opt/dhcal/lcio/v02-00/lib )
LINK_DIRECTORIES(/home/mirabito/SDHCAL/opt/dhcal/lib)
LINK_LIBRARIES(mysqlclient)
LINK_LIBRARIES(pthread )
LINK_LIBRARIES(z )
LINK_LIBRARIES(m )
LINK_LIBRARIES(rt )
LINK_LIBRARIES(dl )
LINK_LIBRARIES(SDHCALcommon )
LINK_LIBRARIES(SDHCALAnalysis )
LINK_LIBRARIES(boost_system )
LINK_LIBRARIES(boost_thread )
LINK_LIBRARIES(lcio )
LINK_LIBRARIES(xml2 )
LINK_LIBRARIES(sqlite3 )
LINK_LIBRARIES(X11 )
LINK_LIBRARIES(jsoncpp)
# add specific flags to compilation from dqm4hep soft
ADD_DEFINITIONS ( ${DQMCore_DEFINITIONS} )
# specify where our includes are
INCLUDE_DIRECTORIES( source )
# build our shared library libModuleExample.so (.dylib)
ADD_SHARED_LIBRARY( ExampleModule source/ExampleModule.cc source/DQMTrackAnalyzer.cc source/DQMSDHCALMonitor.cc )
# install the shared library in the lib directory
INSTALL(
TARGETS ExampleModule
LIBRARY DESTINATION lib
)
