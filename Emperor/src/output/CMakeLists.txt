#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
file(GLOB SOURCES Message.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EMPIRE_MAPPER_LIB_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EMPIRE_MAPPER_LIB_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#