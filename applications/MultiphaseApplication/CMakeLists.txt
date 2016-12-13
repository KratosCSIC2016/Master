set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiphaseApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
 
## generate variables with the sources
set( KRATOS_MULTIPHASE_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/multiphase_application.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/porous_media_application_variables.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiphase_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/flowpressuretrans_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/multiphasefem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksourcepressure.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/massflow.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksourcenonfem.cpp
) 

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMultiphaseApplication SHARED ${KRATOS_MULTIPHASE_APPLICATION_SOURCES})
target_link_libraries(KratosMultiphaseApplication KratosCore )
set_target_properties(KratosMultiphaseApplication PROPERTIES PREFIX "")
install(TARGETS KratosMultiphaseApplication DESTINATION libs )

#QT5#########################################
target_link_libraries(KratosMultiphaseApplication ${QT_LIBRARIES_1})
target_link_libraries(KratosMultiphaseApplication ${QT_LIBRARIES_2})

#qt5_use_modules(KratosMultiphaseApplication Core Xml)
target_link_libraries(KratosMultiphaseApplication Qt5::Core Qt5::Xml)	  
#########################################

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMultiphaseApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiphaseApplication.py" DESTINATION KratosMultiphysics )
 
