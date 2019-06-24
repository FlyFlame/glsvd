# Install script for directory: /scratch/glu19b/code/GKlib2

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/scratch/glu19b/code/GKlib2/build/Linux-x86_64/libGKlib.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/scratch/glu19b/code/GKlib2/gk_arch.h"
    "/scratch/glu19b/code/GKlib2/gk_struct.h"
    "/scratch/glu19b/code/GKlib2/gk_mkutils.h"
    "/scratch/glu19b/code/GKlib2/ms_stdint.h"
    "/scratch/glu19b/code/GKlib2/gk_mkblas.h"
    "/scratch/glu19b/code/GKlib2/gk_mkpqueue2.h"
    "/scratch/glu19b/code/GKlib2/gk_mkpqueue.h"
    "/scratch/glu19b/code/GKlib2/ms_stat.h"
    "/scratch/glu19b/code/GKlib2/GKlib.h"
    "/scratch/glu19b/code/GKlib2/ms_inttypes.h"
    "/scratch/glu19b/code/GKlib2/gk_proto.h"
    "/scratch/glu19b/code/GKlib2/gk_macros.h"
    "/scratch/glu19b/code/GKlib2/gkregex.h"
    "/scratch/glu19b/code/GKlib2/gk_mksort.h"
    "/scratch/glu19b/code/GKlib2/gk_getopt.h"
    "/scratch/glu19b/code/GKlib2/gk_defs.h"
    "/scratch/glu19b/code/GKlib2/gk_types.h"
    "/scratch/glu19b/code/GKlib2/gk_externs.h"
    "/scratch/glu19b/code/GKlib2/gk_mkmemory.h"
    "/scratch/glu19b/code/GKlib2/gk_mkrandom.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/scratch/glu19b/code/GKlib2/build/Linux-x86_64/test/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/scratch/glu19b/code/GKlib2/build/Linux-x86_64/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/scratch/glu19b/code/GKlib2/build/Linux-x86_64/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
