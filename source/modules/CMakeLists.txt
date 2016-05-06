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

# Macro definition to add a module library
MACRO( ADD_MODULE  module_name )
	
	INCLUDE_DIRECTORIES( ${module_name} )
	
	ADD_SHARED_LIBRARY( ${module_name} ${module_name}/${module_name}.cc )
	TARGET_LINK_LIBRARIES( ${module_name} ${PROJECT_NAME}Common )
	
	INSTALL(
	    TARGETS ${module_name}
	    LIBRARY DESTINATION lib
	)

ENDMACRO()


ADD_MODULE( BeamAnalysisModule )
ADD_MODULE( RawAnalysisModule )
ADD_MODULE( HitAnalysisModule )
ADD_MODULE( AsicAnalysisModule )
ADD_MODULE( SlowControlModule )
ADD_MODULE( HoughTransformModule )
ADD_MODULE( EventDisplayModule )
ADD_MODULE( ParticleIDModule )

