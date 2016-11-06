cmake_minimum_required(VERSION 2.8)

####
# project settings
####

project(Cones)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	cones.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(Cones ${HEADERS} ${SOURCES})

target_link_libraries(Cones ${GLAMER_LIBRARIES})


