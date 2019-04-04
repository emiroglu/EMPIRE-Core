#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
FILE(GLOB Emperor_src_sub_subdirectories "*")
FOREACH(sf ${Emperor_src_sub_subdirectories})
  IF(NOT sf MATCHES "\\.txt" AND IS_DIRECTORY ${sf} )
    add_subdirectory(${sf})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#