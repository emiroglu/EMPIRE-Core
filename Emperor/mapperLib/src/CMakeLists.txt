#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EMPIRE_MAPPER_LIB_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
# add_subdirectory(include)
FILE(GLOB EMPIRE_MAPPER_LIB_src_subdirectories "*")
FOREACH(f ${EMPIRE_MAPPER_LIB_src_subdirectories})
  IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
    add_subdirectory(${f})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
get_property(EMPIRE_MAPPER_LIB_INCLUDES GLOBAL PROPERTY EMPIRE_MAPPER_LIB_INCLUDES) 
get_property(EMPIRE_MAPPER_LIB_SOURCES GLOBAL PROPERTY EMPIRE_MAPPER_LIB_SOURCES) 
#------------------------------------------------------------------------------------#
get_property(Emperor_INCLUDES GLOBAL PROPERTY Emperor_INCLUDES)
get_property(EMPIRE_thirdparty_INCLUDES GLOBAL PROPERTY EMPIRE_thirdparty_INCLUDES) 
#------------------------------------------------------------------------------------#
include_directories(${MKL_INCLUDE_DIR})
include_directories(${Emperor_INCLUDES})
include_directories(${EMPIRE_thirdparty_INCLUDES})
include_directories(${EMPIRE_MAPPER_LIB_INCLUDES})
#------------------------------------------------------------------------------------#
# Compiler flags for EMPIRE_MAPPER_LIB
#------------------------------------------------------------------------------------#
add_library(EMPIRE_MAPPER_LIB_Static STATIC ${EMPIRE_MAPPER_LIB_SOURCES})
add_library(EMPIRE_MAPPER_LIB_Shared SHARED ${EMPIRE_MAPPER_LIB_SOURCES})
#------------------------------------------------------------------------------------#
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
#------------------------------------------------------------------------------------#
set_target_properties(EMPIRE_MAPPER_LIB_Static PROPERTIES OUTPUT_NAME EMPIRE_MapperLib)
set_target_properties(EMPIRE_MAPPER_LIB_Shared PROPERTIES OUTPUT_NAME EMPIRE_MapperLib)
#------------------------------------------------------------------------------------#
add_dependencies(EMPIRE_MAPPER_LIB_Static EmperorLib EMPIRE_thirdparty)
add_dependencies(EMPIRE_MAPPER_LIB_Shared EmperorLib EMPIRE_thirdparty)
#------------------------------------------------------------------------------------#
target_link_libraries(EMPIRE_MAPPER_LIB_Static EmperorLib ${MKL_LIBRARIES} ${Emperor_LIBS}) #${MKL_INTEL_LP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} #${Emperor_LIBS} #${MKL_LIBRARIES}
target_link_libraries(EMPIRE_MAPPER_LIB_Shared EmperorLib ${MKL_LIBRARIES} ${Emperor_LIBS}) #${MKL_INTEL_LP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} #${Emperor_LIBS} #${MKL_LIBRARIES}
#------------------------------------------------------------------------------------#
# add_dependencies(EMPIRE_MAPPER_LIB_Static EmperorLib EMPIRE_thirdparty) # ${MKL_LIBRARIES}
# add_dependencies(EMPIRE_MAPPER_LIB_Shared EmperorLib EMPIRE_thirdparty) # ${MKL_LIBRARIES}
#------------------------------------------------------------------------------------#
# target_link_libraries(EMPIRE_MAPPER_LIB_Static EmperorLib) #${MKL_INTEL_LP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} #${Emperor_LIBS} #${MKL_LIBRARIES}
# target_link_libraries(EMPIRE_MAPPER_LIB_Shared EmperorLib) #${MKL_INTEL_LP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} #${Emperor_LIBS} #${MKL_LIBRARIES}
#------------------------------------------------------------------------------------#
# MESSAGE("Emperor_INCLUDES: ${Emperor_INCLUDES}")
# MESSAGE("Emperor_INCLUDES: ${Emperor_SOURCES}")
MESSAGE("C Compiler Flags for EMPIRE_MAPPER_LIB are:\t${CMAKE_C_FLAGS} \nC++ Compiler Flags for EMPIRE_MAPPER_LIB are:\t${CMAKE_CXX_FLAGS}")
# MESSAGE("Emperor_LIBS: ${Emperor_LIBS}")
# 
# MESSAGE("EMPIRE_MAPPER_LIB_INCLUDES: ${EMPIRE_MAPPER_LIB_INCLUDES}")
# MESSAGE("EMPIRE_MAPPER_LIB_SOURCES: ${EMPIRE_MAPPER_LIB_SOURCES}")
# 
# MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
# MESSAGE("MKL_INTEL_LP64_LIBRARY: ${MKL_INTEL_LP64_LIBRARY}")
# MESSAGE("MKL_INTEL_THREAD_LIBRARY: ${MKL_INTEL_THREAD_LIBRARY}")
# MESSAGE("MKL_CORE_LIBRARY: ${MKL_CORE_LIBRARY}")
