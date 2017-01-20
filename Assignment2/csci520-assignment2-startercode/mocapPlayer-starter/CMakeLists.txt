cmake_minimum_required (VERSION 3.1)
project (Assignment2)

#########################################################
# SET COMPILER 
#########################################################
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

#########################################################
# INTERPOLATOR FILE VARIABLES
#########################################################
SET(INT_SOURCE
    motion.cpp
    posture.cpp
    skeleton.cpp
    transform.cpp
    vector.cpp
    interpolator.cpp
    quaternion.cpp
    interpolate.cpp
)

SET(INT_HEADERS
    motion.h
    posture.h
    skeleton.h
    transform.h
    vector.h
    interpolator.h
    quaternion.h
)

#########################################################
# ADD EXECUTABLES
#########################################################
add_executable(interpolator ${INT_SOURCE} ${INT_HEADERS})
