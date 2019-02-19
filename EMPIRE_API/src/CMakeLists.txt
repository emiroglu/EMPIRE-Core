#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EMPIRE_API_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
FILE(GLOB EMPIRE_API_SOURCES_src_subdirectories "*")
FOREACH(f ${EMPIRE_API_SOURCES_src_subdirectories})
  IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
    add_subdirectory(${f})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
get_property(EMPIRE_API_INCLUDES GLOBAL PROPERTY EMPIRE_API_INCLUDES) 
get_property(EMPIRE_API_SOURCES GLOBAL PROPERTY EMPIRE_API_SOURCES) 
get_property(EMPIRE_thirdparty_INCLUDES GLOBAL PROPERTY EMPIRE_thirdparty_INCLUDES) 
#------------------------------------------------------------------------------------#
include_directories("${EMPIRE_API_INCLUDES};${EMPIRE_thirdparty_INCLUDES}")
#------------------------------------------------------------------------------------#
# Compiler flags for EMPIRE_API
#------------------------------------------------------------------------------------#
add_library(EMPIRE_API_Shared SHARED ${EMPIRE_API_SOURCES} $<TARGET_OBJECTS:EMPIRE_thirdparty>)
add_library(EMPIRE_API_Static STATIC ${EMPIRE_API_SOURCES} $<TARGET_OBJECTS:EMPIRE_thirdparty>)
#------------------------------------------------------------------------------------#
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
#------------------------------------------------------------------------------------#
set_target_properties(EMPIRE_API_Shared PROPERTIES OUTPUT_NAME EMPIRE_API)
set_target_properties(EMPIRE_API_Static PROPERTIES OUTPUT_NAME EMPIRE_API)
#------------------------------------------------------------------------------------#
# Add dependencies to targets
add_dependencies(EMPIRE_API_Shared EMPIRE_thirdparty)
add_dependencies(EMPIRE_API_Static EMPIRE_thirdparty)
#------------------------------------------------------------------------------------#
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  target_link_libraries(EMPIRE_API_Shared)
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  target_link_libraries(EMPIRE_API_Shared "-ldl")
ELSE()
   MESSAGE(FATAL_ERROR "Your CXX compiler is ${CMAKE_CXX_COMPILER_ID} and your C compiler is ${CMAKE_C_COMPILER_ID} they are not supported by the build system yet")
ENDIF()
target_link_libraries(EMPIRE_API_Static)
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
MESSAGE("C Compiler Flags for EMPIRE_API are:\t${CMAKE_C_FLAGS} \nC++ Compiler Flags for EMPIRE_API are:\t${CMAKE_CXX_FLAGS}")
MESSAGE("Emperor_LIBS: ${Emperor_LIBS}")
MESSAGE("EMPIRE_API_INCLUDES: ${EMPIRE_API_INCLUDES}")
MESSAGE("EMPIRE_API_SOURCES: ${EMPIRE_API_SOURCES}")
